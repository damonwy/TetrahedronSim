#define TETLIBRARY
#include <iostream>
#include <tetgen.h>
#include <cmath>        // std::abs

#include "ChronoTimer.h"
#include "FemTet.h"
#include "Particle.h"
#include "Program.h"
#include "GLSL.h"
#include "MatrixStack.h"
#include "Spring.h"
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <igl/mosek/mosek_quadprog.h>
#include <omp.h>
#include "stdlib.h"
#define CRTDBG_MAP_ALLOC

#include <crtdbg.h>
#include <Eigen/Core>
#include "MatrixReplacement.h"

using namespace std;
using namespace Eigen;

typedef Eigen::Triplet<double> T;

#define COROTATED_LINEAR 1
#define STVK 2
#define NEOHOOKEAN 3
#define LINEAR 4

// Demo Linear:
//#define YOUNG 100
//#define POISSON 0.4
//#define MOSEK 0
//#define MATERIAL 4

// Demo STVK
//#define YOUNG 100
//#define POISSON 0.4
//#define MOSEK 0
//#define MATERIAL 2

// Demo Corotated 
//#define YOUNG 100
//#define POISSON 0.4
//#define MOSEK 0
//#define MATERIAL 1

// Demo NeoHookean
#define YOUNG 100
#define POISSON 0.4
#define MOSEK 0
#define MATERIAL 3

FemTet::FemTet(double density, const Eigen::Vector2d &damping)
{
	assert(density > 0.0);

	this->damping = damping;
	this->young = YOUNG;
	this->poisson = POISSON;

	// Compute Lame coefficients: [mu, lambda]
	this->mu = young / (2.0 * (1.0 + poisson));
	this->lambda = young * poisson / ((1.0 + poisson)*(1.0 - 2.0 * poisson));

	double r = 0.02; // Used for collisions

	// Load mesh from files

	//in_2.load_off("octtorus");
	//in_2.load_off("Y");
	//in_2.load_ply("tetrahedron");
	//in_2.load_ply("icosahedron");
	//in_2.load_ply("dodecahedron");
	in_2.load_ply("bunny250");
	tetrahedralize("pqz", &in_2, &out);

	// Get tets info
	nVerts = out.numberofpoints;
	nTets = out.numberoftetrahedra;
	nFacets = out.numberoffacets;
	nTriFaces = out.numberoftrifaces;

	// Build system matrices and vectors
	n = 3 * nVerts;

	M.resize(n, n);
	Bm.resize(3*nTets, 3);

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
		Bm.block(itet*3, 0, 3, 3) = Dm.inverse();
	}

	//std::cout << "Finish computing mass and Bm!" << endl;

	// Create particles
	//std::cout << "Start creating particles..." << endl;
	for (int i = 0; i < nVerts; i++) {
		auto p = make_shared<Particle>();
		particles.push_back(p);
		p->r = r;
		// Init postion
		p->x0 << out.pointlist[3 * i], out.pointlist[3 * i + 1], out.pointlist[3 * i + 2];
		p->x = p->x0;
		// Init velocity
		p->v0 << 0.0, 0.0, 0.0;
		p->v = p->v0;
		// Init mass
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
	//posBuf.resize(nVerts * 3);
	posBuf.resize(nTriFaces * 9);
	//norBuf.resize(nVerts * 3);
	norBuf.resize(nTriFaces * 9);
	eleBuf.resize(nTriFaces * 3);
	updatePosNor();

	for (int i = 0; i < nTriFaces; i++) {
		//eleBuf[3 * i] = out.trifacelist[3 * i];
		//eleBuf[3 * i + 1] = out.trifacelist[3 * i + 1];
		//eleBuf[3 * i + 2] = out.trifacelist[3 * i + 2];
		eleBuf[3 * i + 0] = 3 * i;
		eleBuf[3 * i + 1] = 3 * i + 1;
		eleBuf[3 * i + 2] = 3 * i + 2;

	}
}

Matrix3d FemTet::computeStress(Eigen::Matrix3d dp, int material_type, double mu, double lambda) {
	Matrix3d E;
	E.setZero();
	double psi;
	Matrix3d P;
	P.setZero();

	switch (material_type) {
	case COROTATED_LINEAR:
	{
		// Polar decomposition
		Matrix3d A = dp.adjoint() * dp;
		SelfAdjointEigenSolver<Matrix3d> es(A);
		Matrix3d S = es.operatorSqrt();
		Matrix3d R = dp * S.inverse();
		
		E = S - I;
		//psi = mu * E.norm() * E.norm() + 1.0 / 2.0 * lambda * E.trace() * E.trace();
		//P = R * (2.0 * mu * E + lambda * E.trace() * I);
		P = 2.0 * mu *(dp - R) + lambda * (R.transpose()*dp - I).trace() * R;
		break;
	}

	case STVK:
	{	
		E = 1.0 / 2.0 * (dp.transpose() * dp - I);
		//psi = mu * E.norm()*E.norm() + 1.0 / 2.0 * lambda * E.trace() * E.trace();
		P = dp * (2.0 * mu * E + lambda * E.trace() * I);
		break;
	}

	case NEOHOOKEAN:
	{
		double I1 = (dp.transpose() * dp).trace();
		double I2 = ((dp.transpose() * dp) *  (dp.transpose() * dp)).trace();
		double I3 = (dp.transpose() * dp).determinant();
		double J = sqrt(I3);
		//psi = 1.0 / 2.0 * mu *(I1 - 3.0) - mu * log(J) + 1.0 / 2.0 * lambda * log(J)*log(J);
		P = mu * (dp - dp.inverse().transpose()) + lambda * log(J)*(dp.inverse().transpose());
		break;
	}
	case LINEAR:
	{
		E = 1.0 / 2.0 * (dp + dp.transpose()) - I;
		//psi = mu * E.norm() * E.norm() + 1.0 / 2.0 * lambda * E.trace() * E.trace();
		P = 2.0 * mu * E + lambda * E.trace() * I;
		break;
	}
	default:
		break;
	}
	return P;
}

Matrix3d FemTet::computeStressDerivative(Eigen::Matrix3d F, Eigen::Matrix3d dF, int material_type, double mu, double lambda) {
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
		//double I1 = (F.norm()) * (F.norm());
		//double I2 = ((F.transpose() * F) *  (F.transpose() * F)).trace();
		MatrixXd FT = F.transpose();
		MatrixXd FIT = F.inverse().transpose();
		//double I3 = (F.transpose() * F).determinant();
		double I3 = (FT * F).determinant();
		double J = sqrt(I3);
		P = mu * (F - FIT) + lambda * log(J) * FIT;
		dP = mu * dF + (mu - lambda * log(J)) * FIT * (dF.transpose()) * FIT + lambda * ((F.inverse() * dF)).trace() * FIT;
		
		//P = mu * (F - (F.inverse().transpose())) + lambda * log(J) * (F.inverse().transpose());
		//dP = mu * dF + (mu - lambda * log(J)) * (F.inverse().transpose()) * (dF.transpose()) * (F.inverse().transpose()) + lambda * ((F.inverse() * dF)).trace() * (F.inverse().transpose());
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

MatrixXd FemTet::hooke(double young, double poisson) {
	double v = poisson;
	double v1 = 1 - v;
	double v2 = 1 - 2 * v;
	double s = young / ((1 + v)*(1 - 2 * v));

	MatrixXd E(6, 6);
	E <<
		v1, v, v, 0, 0, 0,
		v, v1, v, 0, 0, 0,
		v, v, v1, 0, 0, 0,
		0, 0, 0, v2, 0, 0,
		0, 0, 0, 0, v2, 0,
		0, 0, 0, 0, 0, v2;
	E = s*E;
	return E;
}

void FemTet::step(double h, const Vector3d &grav) {
	ChronoTimer timer("test", 1);
	timer.tic();

	Eigen::initParallel();
	omp_set_num_threads(4);
	Eigen::setNbThreads(4);
	M.setZero();
	v.setZero();
	f.setZero();
	X.setZero();
	
	vector<T> A_;

#pragma omp parallel for shared(A_)
	for (int i = 0; i < particles.size(); i++) {
		int idx = particles[i]->i;
		double mass = particles[i]->m;
		// filling sparse matrix A_ with mass
		M.block<3, 3>(3 * idx, 3 * idx) = I * mass; // filling M
		v.segment<3>(3 * idx) = particles[i]->v; // filling v 
		f.segment<3>(3 * idx) = mass * grav; // filling f with fg
		X.segment<3>(3 * idx) = particles[i]->x;// filling X with x
		
#pragma omp critical 
		for (int j = 0; j < 3; j++) {
			A_.push_back(T(3 * i + j, 3 * i + j, mass + h * damping(0) * mass));
		}		
	}
	

	//cout << "Finishing filling matrices" << endl;
	//timer.toc();
	//timer.print();
#pragma omp parallel for
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
	
	//cout << "Finishing computing forces" << endl;
	//timer.toc();
	//timer.print();

	// Compute K matrix and fill the sparse matrix A_

	if (MATERIAL != 1) {
		MatrixXd w(n, n);
		w.setIdentity();
		double co = h * h * damping(1);


#pragma omp parallel for 
		for (int icol = 0; icol < n; icol++) {
			VectorXd dx = w.col(icol);
			VectorXd df(n);
			df.setZero();

#pragma omp parallel for shared(A_) 
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

				// Compute dDs
				Matrix3d dDs;
				Vector3d dpb = dx.segment<3>(3 * i) - dx.segment<3>(3 * l);
				Vector3d dpc = dx.segment<3>(3 * j) - dx.segment<3>(3 * l);
				Vector3d dpd = dx.segment<3>(3 * k) - dx.segment<3>(3 * l);
				dDs.col(0) = dpb;
				dDs.col(1) = dpc;
				dDs.col(2) = dpd;

				// Compute F
				Matrix3d F = Ds * (Bm.block(3 * itet, 0, 3, 3));

				// Compute dF
				Matrix3d dF = dDs * (Bm.block(3 * itet, 0, 3, 3));

				// Compute dP
				Matrix3d dP = computeStressDerivative(F, dF, MATERIAL, mu, lambda);

				// Compute dH
				Matrix3d dH = -W(itet) * dP * (Bm.block(3 * itet, 0, 3, 3)).transpose();

				Vector3d dh1 = dH.col(0);
				Vector3d dh2 = dH.col(1);
				Vector3d dh3 = dH.col(2);
				Vector3d dh4 = -dh1 - dh2 - dh3;

				/*df.segment<3>((3 * i)) += dh1;
				df.segment<3>((3 * j)) += dh2;
				df.segment<3>((3 * k)) += dh3;
				df.segment<3>((3 * l)) += dh4;*/

				#pragma omp critical 
				// Filling the A_ matrix
				for (int idx = 0; idx < 3; idx++) {
					A_.push_back(T(3 * i + idx, icol, -co * dh1(idx)));
					A_.push_back(T(3 * j + idx, icol, -co * dh2(idx)));
					A_.push_back(T(3 * k + idx, icol, -co * dh3(idx)));
					A_.push_back(T(3 * l + idx, icol, -co * dh4(idx)));
				}
			}
			// Filling the ith col of K matrix
			//K.col(icol) = df;	
		}
	}

	//cout << "Finishing computing K matrix" << endl;
	//timer.toc();
	//timer.print();
	// Sparse Matrix Version ....
	
	Eigen::VectorXd b;
	b.resize(n);
	b.setZero();
	b = M * v + h * f;

	Eigen::SparseMatrix<double> A(n, n);
	A.setFromTriplets(A_.begin(), A_.end());

	//MatrixReplacement MR;
	//MR.attachMyMatrix(A);
	//
	//// Solve Ax = b using various iterative solver with matrix-free version:
	//
	//Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> cg2;
	//cg2.setTolerance(25);
	//cg2.setTolerance(1e-3);
	//cg2.compute(MR);
	//VectorXd v_new = cg2.solve(b);
	//std::cout << "CG:       #iterations: " << cg2.iterations() << ", estimated error: " << cg2.error() << std::endl;

	ConjugateGradient< SparseMatrix<double> > cg;
	cg.setMaxIterations(25);
	cg.setTolerance(1e-3);
	cg.compute(A);

	VectorXd v_new = cg.solveWithGuess(b, v);
	
	//cout << "Finishing computing velocity" << endl;

	//timer.toc();
	//timer.print();
	
	// Common ....
	/*MatrixXd LHS(n, n);
	LHS.setZero();
	LHS = M + h * damping(0) * M - h * h * damping(1) * K;
	VectorXd RHS(n);
	RHS.setZero();
	RHS = M * v + h * f;

	VectorXd result = (LHS).ldlt().solve(RHS);
	VectorXd v_new = result;*/

#pragma omp parallel for
	// Update velocity
	for (int i = 0; i < particles.size(); i++) {
		if (particles[i]->i != -1) {
			particles[i]->v = v_new.segment<3>(3 * particles[i]->i);
		}
	}
#pragma omp parallel for
	// Update position
	for (int i = 0; i < particles.size(); i++) {
		if (particles[i]->i != -1) {
			particles[i]->x += particles[i]->v * h;
		}
	}

	// Collision with floor without using mosek
	if (MOSEK == 0) {
#pragma omp parallel for
		for (int i = 0; i < particles.size(); i++) {
			if (particles[i]->x(1) < -2.0 && particles[i]->v(1) < -0.001) {
				particles[i]->v(1) = 0.0;
				particles[i]->x(1) = -2.0;
			}
		}
		updatePosNor();
	}

	// Collision with the floor with using mosek
	if (MOSEK == 1) {

		// for constraints matrix
		vector<T> C_;
		int row = 0;

		double y_floor = -1.0;
		Vector3d floor_normal(0.0, 1.0, 0.0);

		for (int i = 0; i < particles.size(); i++) {

			shared_ptr<Particle> p = particles[i];

			Vector3d x_i = p->x;
			Vector3d v_i = p->v;
			double xi = x_i[1];

			if (xi < y_floor && v_i(1) < 0.000001) {
				// filling the constraint matrix
				for (int t = 0; t < 3; t++) {
					C_.push_back(T(row, 3 * i + t, floor_normal(t)));
				}

				row++;
			}
		}

		if (row != 0) {
			cout << "floor!" << endl;
			Eigen::SparseMatrix<double> C(row, n);
			C.setFromTriplets(C_.begin(), C_.end());
			VectorXd lc = VectorXd::Zero(row); // linear inequality
			VectorXd uc = VectorXd::Zero(row); // +Inifity
			VectorXd lx = VectorXd::Zero(n);
			VectorXd ux = VectorXd::Zero(n);
			VectorXd cc = -b;

			for (int i = 0; i < n; i++) {
				lx(i) = -Infinity;
				ux(i) = +Infinity;
			}

			for (int i = 0; i < row; i++) {
				uc(i) = +Infinity;
			}

			VectorXd results;
			results.resize(n);
			igl::mosek::MosekData mosek_data;
			bool r = mosek_quadprog(A, cc, 0, C, lc, uc, lx, ux, mosek_data, results);

			// Update velocity
			for (int i = 0; i < particles.size(); i++) {
				if (particles[i]->i != -1) {
					particles[i]->v = results.segment<3>(3 * particles[i]->i);
				}
			}

			// Update position
			for (int i = 0; i < particles.size(); i++) {
				if (particles[i]->i != -1) {
					particles[i]->x += particles[i]->v * h;
				}
			}
			updatePosNor();
		}
	}
	//cout << "Finishing updating velocity and position " << endl;
	//.toc();
	//timer.print();
}

void FemTet::init() {
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

void FemTet::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> p) const
{

	// Draw mesh
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

void FemTet::tare() {
	for (int k = 0; k < (int)particles.size(); k++) {
		particles[k]->tare();
	}
}

void FemTet::reset() {
	for (int k = 0; k < (int)particles.size(); k++) {
		particles[k]->reset();
	}
	updatePosNor();
}

void FemTet::updatePosNor()
{

	// Position
	/*for (int i = 0; i < (int)particles.size(); i++) {
		Vector3d x = particles[i]->x;
		posBuf[3 * i + 0] = x(0);
		posBuf[3 * i + 1] = x(1);
		posBuf[3 * i + 2] = x(2);
	}*/

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

		//particles[out.trifacelist[3 * iface]]->normal += normal;
		//particles[out.trifacelist[3 * iface + 1]]->normal += normal;
		//particles[out.trifacelist[3 * iface + 2]]->normal += normal;
	}

	/*for (int ipt = 0; ipt < particles.size(); ipt++) {
		Vector3d nor = particles[ipt]->normal;
		nor.normalize();
		norBuf[3 * ipt + 0] = nor(0);
		norBuf[3 * ipt + 1] = nor(1);
		norBuf[3 * ipt + 2] = nor(2);
	}*/
}

FemTet::~FemTet()
{
}
