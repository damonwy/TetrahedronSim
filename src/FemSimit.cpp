#define TETLIBRARY
#include <iostream>
#include <tetgen.h>

#include "ChronoTimer.h"
#include "FemSimit.h"
#include "Particle.h"
#include "Program.h"
#include "GLSL.h"
#include "MatrixStack.h"
#include "Spring.h"
#include "stdlib.h"

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <igl/mosek/mosek_quadprog.h>
#include <omp.h>
#include <Eigen/Core>
#include <cmath>        // std::abs

using namespace std;
using namespace Eigen;
typedef Eigen::Triplet<double> T;

#define COROTATED_LINEAR 1
#define STVK 2
#define NEOHOOKEAN 3
#define LINEAR 4

// Demo NeoHookean
#define YOUNG 100
#define POISSON 0.4
#define MOSEK 0
#define MATERIAL 3

FemSimit::FemSimit(double density, const Eigen::Vector2d &damping)
{
	assert(density > 0.0);

	this->damping = damping;
	this->young = YOUNG;
	this->poisson = POISSON;

	// Compute Lame coefficients: [mu, lambda]
	this->mu = young / (2.0 * (1.0 + poisson));
	this->lambda = young * poisson / ((1.0 + poisson)*(1.0 - 2.0 * poisson));

	double r = 0.02; 

	//in_2.load_off("octtorus");
	//in_2.load_off("Y");
	//in_2.load_ply("tetrahedron");
	//in_2.load_ply("icosahedron");
	//in_2.load_ply("dodecahedron");

	in_2.load_ply("bunny3");
	tetrahedralize("pqz", &in_2, &out);

	// Get tets info
	nVerts = out.numberofpoints;
	nTets = out.numberoftetrahedra;
	nFacets = out.numberoffacets;
	nTriFaces = out.numberoftrifaces;

	// Build system matrices and vectors
	n = 3 * nVerts;

	M.resize(n, n);
	Bm.resize(3 * nTets, 3);

	v.resize(n);

	mass.resize(nVerts);
	mass.setZero();

	f.resize(n);
	X.resize(n);

	W.resize(nTets);
	W.setZero();

	K.resize(n, n);
	K.setZero();

	P.resize(3 * nTets, 3);
	P.setZero();

	I.setIdentity();

	// Compute tet mass and distribute to vertices
	//std::cout << "Start computing mass and Bm..." << endl;

	for (int itet = 0; itet < nTets; itet++) {
		Vector3d xa, xb, xc, xd;
		int a = out.tetrahedronlist[itet * 4];
		int b = out.tetrahedronlist[itet * 4 + 1];
		int c = out.tetrahedronlist[itet * 4 + 2];
		int d = out.tetrahedronlist[itet * 4 + 3];

		int ia = 3 * a;
		xa << out.pointlist[ia], out.pointlist[ia + 1], out.pointlist[ia + 2];

		int ib = 3 * b;
		xb << out.pointlist[ib], out.pointlist[ib + 1], out.pointlist[ib + 2];

		int ic = 3 * c;
		xc << out.pointlist[ic], out.pointlist[ic + 1], out.pointlist[ic + 2];

		int id = 3 * d;
		xd << out.pointlist[id], out.pointlist[id + 1], out.pointlist[id + 2];

		// Compute volume and mass of each tet
		double vol = abs((xa - xd).transpose()*(xb - xd).cross(xc - xd)) / 6.0;
		double m = vol * density;
		W(itet) = vol;

		// Distribute 1/4 mass to each vertices
		mass(a) += m / 4.0;
		mass(b) += m / 4.0;
		mass(c) += m / 4.0;
		mass(d) += m / 4.0;

		// Precompute Bm, constant throughout the simulation
		MatrixXd Dm(3, 3);
		Dm.col(0) = xa - xd;
		Dm.col(1) = xb - xd;
		Dm.col(2) = xc - xd;
		Bm.block(itet * 3, 0, 3, 3) = Dm.inverse();
	}

	//std::cout << "Finish computing mass and Bm!" << endl;

	//std::cout << "Start creating particles..." << endl;
	for (int i = 0; i < nVerts; i++) {
		auto p = make_shared<Particle>();
		particles.push_back(p);
		p->r = r;
		p->x0 << out.pointlist[3 * i], out.pointlist[3 * i + 1], out.pointlist[3 * i + 2];
		p->x = p->x0;
		p->v0 << 0.0, 0.0, 0.0;
		p->v = p->v0;
		p->m = mass(i);
		p->i = i;
		p->fixed = false;
	}
	//std::cout << "Finish creating particles!" << endl;

	// Build vertex buffers
	posBuf.clear();
	norBuf.clear();
	texBuf.clear();
	eleBuf.clear();
	posBuf.resize(nTriFaces * 9);
	norBuf.resize(nTriFaces * 9);
	eleBuf.resize(nTriFaces * 3);
	updatePosNor();

	for (int i = 0; i < nTriFaces; i++) {
		eleBuf[3 * i + 0] = 3 * i;
		eleBuf[3 * i + 1] = 3 * i + 1;
		eleBuf[3 * i + 2] = 3 * i + 2;
	}
}

Matrix3d FemSimit::computeStress(Eigen::Matrix3d dp, int material_type, double mu, double lambda) {
	Matrix3d E;
	E.setZero();
	double psi;
	Matrix3d P;
	P.setZero();

	switch (material_type) {
	case COROTATED_LINEAR:
	{
		Matrix3d A = dp.adjoint() * dp;
		SelfAdjointEigenSolver<Matrix3d> es(A);
		Matrix3d S = es.operatorSqrt();
		Matrix3d R = dp * S.inverse();

		E = S - I;
		P = 2.0 * mu *(dp - R) + lambda * (R.transpose()*dp - I).trace() * R;
		break;
	}

	case STVK:
	{
		E = 1.0 / 2.0 * (dp.transpose() * dp - I);
		P = dp * (2.0 * mu * E + lambda * E.trace() * I);
		break;
	}

	case NEOHOOKEAN:
	{
		double I1 = (dp.transpose() * dp).trace();
		double I2 = ((dp.transpose() * dp) *  (dp.transpose() * dp)).trace();
		double I3 = (dp.transpose() * dp).determinant();
		double J = sqrt(I3);
		P = mu * (dp - dp.inverse().transpose()) + lambda * log(J)*(dp.inverse().transpose());
		break;
	}
	case LINEAR:
	{
		E = 1.0 / 2.0 * (dp + dp.transpose()) - I;
		P = 2.0 * mu * E + lambda * E.trace() * I;
		break;
	}
	default:
		break;
	}
	return P;
}

Matrix3d FemSimit::computeStressDerivative(Eigen::Matrix3d F, Eigen::Matrix3d dF, int material_type, double mu, double lambda) {
	Matrix3d E, dE, P, dP;
	E.setZero();
	dE.setZero();
	P.setZero();
	dP.setZero();

	switch (material_type) {
	case COROTATED_LINEAR:
	{
		Matrix3d A = F.adjoint() * F;
		SelfAdjointEigenSolver<Matrix3d> es(A);
		Matrix3d S = es.operatorSqrt();
		Matrix3d R = F * S.inverse();
		E = S - I;
		P = 2.0 * mu *(F - R) + lambda * (R.transpose()*F - I).trace() * R;
		break;
	}

	case STVK:
	{
		E = 1.0 / 2.0 * (F.transpose() * F - I);
		dE = 1.0 / 2.0 * (dF.transpose() * F + F.transpose() * dF);
		P = F * (2.0 * mu * E + lambda * E.trace() * I);
		dP = dF * (2.0 * mu * E + lambda * E.trace() * I) + F * (2.0 * mu * dE + lambda * dE.trace() * I);
		break;
	}

	case NEOHOOKEAN:
	{
		MatrixXd FT = F.transpose();
		MatrixXd FIT = F.inverse().transpose();
		double I3 = (FT * F).determinant();
		double J = sqrt(I3);
		P = mu * (F - FIT) + lambda * log(J) * FIT;
		dP = mu * dF + (mu - lambda * log(J)) * FIT * (dF.transpose()) * FIT + lambda * ((F.inverse() * dF)).trace() * FIT;
	break;
	}
	case LINEAR:
	{
		E = 1.0 / 2.0 * (F + F.transpose()) - I;
		dE = 1.0 / 2.0 * (dF + dF.transpose());
		P = 2.0 * mu * E + lambda * E.trace() * I;
		dP = 2.0 * mu * dE + lambda * dE.trace() * I;
		break;
	}
	default:
		break;
	}
	return dP;
}

void FemSimit::step(double h, const Vector3d &grav) {
	M.setZero();
	v.setZero();
	f.setZero();
	X.setZero();

	vector<T> A_;

	for (int i = 0; i < particles.size(); i++) {
		int idx = particles[i]->i;
		double mass = particles[i]->m;
	
		M.block<3, 3>(3 * idx, 3 * idx) = I * mass; 
		v.segment<3>(3 * idx) = particles[i]->v; 
		f.segment<3>(3 * idx) = mass * grav; 
		X.segment<3>(3 * idx) = particles[i]->x;

		for (int j = 0; j < 3; j++) {
			A_.push_back(T(3 * i + j, 3 * i + j, mass + h * damping(0) * mass));
		}
	}
	
	for (int itet = 0; itet < nTets; itet++) {
		int i = out.tetrahedronlist[itet * 4];
		int j = out.tetrahedronlist[itet * 4 + 1];
		int k = out.tetrahedronlist[itet * 4 + 2];
		int l = out.tetrahedronlist[itet * 4 + 3];

		// Compute Ds
		Matrix3d Ds;
		Vector3d pb = particles[i]->x - particles[l]->x;
		Vector3d pc = particles[j]->x - particles[l]->x;
		Vector3d pd = particles[k]->x - particles[l]->x;
		Ds.col(0) = pb;
		Ds.col(1) = pc;
		Ds.col(2) = pd;

		// Compute Deformation Gradient
		Matrix3d F = Ds * (Bm.block(3 * itet, 0, 3, 3));

		// Compute P 
		Matrix3d p = computeStress(F, MATERIAL, mu, lambda); // Input different materials

		P.block(3 * itet, 0, 3, 3) = p;

		// Compute H
		Matrix3d H;
		H.setZero();
		H = -W(itet) * p * (Bm.block(3 * itet, 0, 3, 3)).transpose();

		Vector3d h1 = H.col(0);
		Vector3d h2 = H.col(1);
		Vector3d h3 = H.col(2);

		f.segment<3>((3 * i)) += h1;
		f.segment<3>((3 * j)) += h2;
		f.segment<3>((3 * k)) += h3;
		f.segment<3>((3 * l)) += (-h1 - h2 - h3);
	}

	if (MATERIAL != 1) {
		MatrixXd w(n, n);
		w.setIdentity();
		double co = h * h * damping(1);
		for (int itet = 0; itet < nTets; itet++) {
			int i = out.tetrahedronlist[itet * 4 + 0];
			int j = out.tetrahedronlist[itet * 4 + 1];
			int k = out.tetrahedronlist[itet * 4 + 2];
			int l = out.tetrahedronlist[itet * 4 + 3];

			// Compute Ds
			Matrix3d Ds;
			Vector3d pb = particles[i]->x - particles[l]->x;
			Vector3d pc = particles[j]->x - particles[l]->x;
			Vector3d pd = particles[k]->x - particles[l]->x;
			Ds.col(0) = pb;
			Ds.col(1) = pc;
			Ds.col(2) = pd;

			Matrix3d F = Ds * (Bm.block(3 * itet, 0, 3, 3));

			MatrixXd dFRow(4, 3);

			for (int i = 0; i < 3; i++) {
				dFRow.row(i) = Bm.row(3 * itet + i);
				dFRow(3, i) = -Bm(3 * itet + 0, i) - Bm(3 * itet + 1, i) - Bm(3 * itet + 2, i);
			}

			for (int row = 0; row < 4; row++) {
				MatrixXd Kb(12, 3);
				Kb.setZero();
				for (int kk = 0; kk < 3; kk++) {
					Matrix3d dF;
					dF.setZero();
					dF.row(kk) = dFRow.row(row);
					MatrixXd dP = computeStressDerivative(F, dF, MATERIAL, mu, lambda);
					MatrixXd dH = -W(itet)*dP*(Bm.block(3 * itet, 0, 3, 3)).transpose();

					// fill Kb
					
					for (int ii = 0; ii < 3; ii++) {
						for (int ll = 0; ll < 3; ll++) {
							Kb(ii * 3 + ll, kk) = dH(ll, ii);
						}
						Kb(9+ii, kk) = -dH(ii, 0) - dH(ii, 1) - dH(ii, 2);
					}
				}

				for (int jj = 0; jj < 4; jj++) {
					int aa = (out.tetrahedronlist[itet * 4 + jj]);
					int bb = (out.tetrahedronlist[itet * 4 + row]);

					K.block(3 * aa, 3 * bb, 3, 3) += Kb.block(3 * jj, 0, 3, 3);
					for (int irow = 0; irow < 3; irow++) {
						for (int icol = 0; icol < 3; icol++) {
							A_.push_back(T(3 * aa + irow, 3 * bb + icol, -co*Kb(3 * jj + irow, icol)));
						}
					}
				}

			}
		}
	}
	
	Eigen::VectorXd b;
	b.resize(n);
	b.setZero();
	b = M * v + h * f;

	Eigen::SparseMatrix<double> A(n, n);
	A.setFromTriplets(A_.begin(), A_.end());

	ConjugateGradient< SparseMatrix<double> > cg;
	cg.setMaxIterations(25);
	cg.setTolerance(1e-3);
	cg.compute(A);

	VectorXd v_new = cg.solveWithGuess(b, v);

	for (int i = 0; i < particles.size(); i++) {
		if (particles[i]->i != -1) {
			particles[i]->v = v_new.segment<3>(3 * particles[i]->i);
		}
	}

	for (int i = 0; i < particles.size(); i++) {
		if (particles[i]->i != -1) {
			particles[i]->x += particles[i]->v * h;
		}
	}

	if (MOSEK == 0) {
		for (int i = 0; i < particles.size(); i++) {
			if (particles[i]->x(1) < 0.0 && particles[i]->v(1) < -0.001) {
				particles[i]->v(1) = 0.0;
				particles[i]->x(1) = 0.0;
			}
		}
		updatePosNor();
	}
}

void FemSimit::init() {
	glGenBuffers(1, &posBufID);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size() * sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);

	/*glGenBuffers(1, &texBufID);
	glBindBuffer(GL_ARRAY_BUFFER, texBufID);
	glBufferData(GL_ARRAY_BUFFER, texBuf.size() * sizeof(float), &texBuf[0], GL_STATIC_DRAW);*/

	glGenBuffers(1, &eleBufID);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, eleBuf.size() * sizeof(unsigned int), &eleBuf[0], GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	assert(glGetError() == GL_NO_ERROR);
}

void FemSimit::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> p) const
{
	glUniform3fv(p->getUniform("kdFront"), 1, Vector3f(1.0, 0.0, 0.0).data());
	glUniform3fv(p->getUniform("kdBack"), 1, Vector3f(1.0, 1.0, 0.0).data());
	MV->pushMatrix();

	glUniformMatrix4fv(p->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));

	int h_pos = p->getAttribute("aPos");
	glEnableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size() * sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);
	glVertexAttribPointer(h_pos, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);

	int h_nor = p->getAttribute("aNor");
	glEnableVertexAttribArray(h_nor);
	glBindBuffer(GL_ARRAY_BUFFER, norBufID);
	glBufferData(GL_ARRAY_BUFFER, norBuf.size() * sizeof(float), &norBuf[0], GL_DYNAMIC_DRAW);
	glVertexAttribPointer(h_nor, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);
	glDrawElements(GL_TRIANGLES, 3 * nTriFaces, GL_UNSIGNED_INT, (const void *)(0 * sizeof(unsigned int)));

	glDisableVertexAttribArray(h_nor);
	glDisableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	MV->popMatrix();
}

void FemSimit::tare() {
	for (int k = 0; k < (int)particles.size(); k++) {
		particles[k]->tare();
	}
}

void FemSimit::reset() {
	for (int k = 0; k < (int)particles.size(); k++) {
		particles[k]->reset();
	}
	updatePosNor();
}

void FemSimit::updatePosNor()
{
	// Normal
	//#pragma omp parallel for
	for (int iface = 0; iface < nTriFaces; iface++) {
		Vector3d p1 = particles[out.trifacelist[3 * iface + 0]]->x;
		Vector3d p2 = particles[out.trifacelist[3 * iface + 1]]->x;
		Vector3d p3 = particles[out.trifacelist[3 * iface + 2]]->x;

		//Position
		Vector3d e1 = p2 - p1;
		Vector3d e2 = p3 - p1;
		Vector3d normal = e1.cross(e2);
		normal.normalize();
	//#pragma omp parallel for
		for (int idx = 0; idx < 3; idx++) {
			posBuf[9 * iface + 0 + idx] = p1(idx);
			posBuf[9 * iface + 3 + idx] = p2(idx);
			posBuf[9 * iface + 6 + idx] = p3(idx);
			norBuf[9 * iface + 0 + idx] = normal(idx);
			norBuf[9 * iface + 3 + idx] = normal(idx);
			norBuf[9 * iface + 6 + idx] = normal(idx);
		}
	}
}

FemSimit::~FemSimit()
{
}