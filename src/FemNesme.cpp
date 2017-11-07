#define TETLIBRARY

#include <iostream>
#include <tetgen.h>
#include <cmath>        // std::abs

#include "FemNesme.h"
#include "Particle.h"
#include "Program.h"
#include "GLSL.h"
#include "MatrixStack.h"
#include "Spring.h"
#include <omp.h>
#include "stdlib.h"
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

using namespace std;
using namespace Eigen;

FemNesme::FemNesme(
	double density,
	const Eigen::Vector2d &damping)
{
	assert(density > 0.0);

	this->damping = damping;
	this->young = 1e3;
	this->poisson = 0.3;

	// Compute Lame coefficients: [mu, lambda]
	this->mu = young / (2.0 * (1.0 + poisson));
	this->lambda = young * poisson / ((1.0 + poisson)*(1.0 - 2.0 * poisson));

	double r = 0.02; // Used for collisions

					//in_2.load_ply("icosahedron");
					 //in_2.load_ply("dodecahedron");
					 //in_2.load_off("octtorus");
					 //in_2.load_off("N");
					 //in_2.load_ply("bunny");
	in_2.load_ply("tetrahedron");

	tetrahedralize("pqz", &in_2, &out);

	nVerts = out.numberofpoints;
	nTets = out.numberoftetrahedra;
	nFacets = out.numberoffacets;
	nTriFaces = out.numberoftrifaces;

	// Build system matrices and vectors
	n = 3 * nVerts;

	M.resize(n, n);
	v.resize(n);
	K.resize(n, n);
	K.setZero();

	mass.resize(nVerts);
	mass.setZero();

	volume.resize(nTets);
	volume.setZero();

	X_invs.resize(3 * nTets, 3);
	F.resize(n);
	D.resize(6, 6);

	// Compute tet mass and distribute to vertices
	std::cout << "Start computing mass and X_inv..." << endl;

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
		volume(itet) = vol;

		double m = vol * density;

		// Distribute 1/4 mass to each vertices
		mass(a) += m / 4;
		mass(b) += m / 4;
		mass(c) += m / 4;
		mass(d) += m / 4;

		// Precompute X_inv, constant throughout the simulation
		MatrixXd _X_inv(3, 3);
		_X_inv.col(0) = xb - xa;
		_X_inv.col(1) = xc - xa;
		_X_inv.col(2) = xd - xa;
		X_invs.block(itet * 3, 0, 3, 3) = _X_inv.inverse();

	}

	// Create particles
	std::cout << "Start creating particles..." << endl;
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

	D = computeD(lambda, mu);

	// Build vertex buffers
	posBuf.clear();
	norBuf.clear();
	texBuf.clear();
	eleBuf.clear();
	posBuf.resize(nVerts * 3);
	norBuf.resize(nVerts * 3);
	eleBuf.resize(nTriFaces * 3);
	updatePosNor();

	for (int i = 0; i <nTriFaces; i++) {
		eleBuf[3 * i] = out.trifacelist[3 * i];
		eleBuf[3 * i + 1] = out.trifacelist[3 * i + 1];
		eleBuf[3 * i + 2] = out.trifacelist[3 * i + 2];
	}

}

void FemNesme::step(double h, const Vector3d &grav) {
	M.setZero();
	v.setZero();
	F.setZero();
	K.setZero();

	for (int i = 0; i < particles.size(); i++) {
		int idx = particles[i]->i;
		double mass = particles[i]->m;

		Matrix3d A;
		A.setIdentity();
		A *= mass;
		M.block<3, 3>(3 * idx, 3 * idx) = A; // filling M
		v.segment<3>(3 * idx) = particles[i]->v; // filling v 
		F.segment<3>(3 * idx) = mass * grav; // filling f with fg
	}

	for (int itet = 0; itet < nTets; itet++) {
		int a = out.tetrahedronlist[itet * 4 + 0];
		int b = out.tetrahedronlist[itet * 4 + 1];
		int c = out.tetrahedronlist[itet * 4 + 2];
		int d = out.tetrahedronlist[itet * 4 + 3];

		MatrixXd dp(3, 3);
		Vector3d pb = particles[b]->x - particles[a]->x;
		Vector3d pc = particles[c]->x - particles[a]->x;
		Vector3d pd = particles[d]->x - particles[a]->x;
		dp.col(0) = pb;
		dp.col(1) = pc;
		dp.col(2) = pd;
		cout << "dp:" << dp << endl;

		// Compute Deformation Gradient
		Matrix3d P = (X_invs.block(3 * itet, 0, 3, 3)) *dp ;
		//cout << P << endl;

		// Compute R and Re using QR decomposition (Gram-Schmidt orthogonalization)

		Vector3d r0 = P.col(0);
		Vector3d r1 = P.col(1);
		r0 = r0.normalized();
		Vector3d r2 = r0.cross(r1);
		r2.normalized();
		r1 = r2.cross(r0);
		r1 = r1.normalized();

		Matrix3d R;
		R.col(0) = r0;
		R.col(1) = r1;
		R.col(2) = r2;
		//cout << "Rotation Matrix"<<endl<< R << endl;
		cout << "M:"<<M << endl;
		// Compute matrix Re
		MatrixXd Re(12, 12);
		Re.setZero();
		for (int i = 0; i < 4; i++) {
			Re.block<3, 3>(3 * i, 3 * i) = R;
		}

		cout << "R" << R << endl<<endl;
		Matrix3d I;
		I.setIdentity();

		// Deformation matrix
		Matrix3d Ee;
		Ee = R.transpose() * P;
		cout << "E:" << endl << Ee << endl;

		VectorXd XX(12);
		XX.setZero();
		XX.segment<3>(0) = particles[a]->x0;
		XX.segment<3>(3) = particles[b]->x0;
		XX.segment<3>(6) = particles[c]->x0;
		XX.segment<3>(9) = particles[d]->x0 ;

		// Fill the initial position x0 to vector xx
		VectorXd xx(12);
		xx.setZero();
		xx.segment<3>(0) = particles[a]->x;
		xx.segment<3>(3) = particles[b]->x;
		xx.segment<3>(6) = particles[c]->x;
		xx.segment<3>(9) = particles[d]->x;

		VectorXd deformed = Re.transpose() * xx;
		//cout << "deformed pos:" << endl << deformed << endl << endl;
		VectorXd deformed_local = deformed;
		//deformed_local.segment<3>(3) -= deformed_local.segment<3>(0);
		//deformed_local.segment<3>(6) -= deformed_local.segment<3>(0);
		//deformed_local.segment<3>(9) -= deformed_local.segment<3>(0);
		//deformed_local.segment<3>(0) -= deformed_local.segment<3>(0);

		VectorXd rest_local = XX;
		//rest_local.segment<3>(3) -= rest_local.segment<3>(0);
		//rest_local.segment<3>(6) -= rest_local.segment<3>(0);
		//rest_local.segment<3>(9) -= rest_local.segment<3>(0);
		//rest_local.segment<3>(0) -= rest_local.segment<3>(0);

		VectorXd displaced = deformed_local - rest_local;
		//cout << "disp: " << endl << displaced << endl;

		//cout << "Rotation Matrix"<<endl<< R << endl;
		// Compute matrix Re
		//MatrixXd Re(12, 12);
		//Re.setZero();
		// Filling Re
		//for (int i = 0; i < 4; i++) {
		//	Re.block<3, 3>(3 * i, 3 * i) = R;
		//}

		// Compute the coefficients of the shape functions N = alpha + beta *x + gamma * y + delta * z (for B, we don't need alpha)
		// In the local frame!!!
		double yc = pc(1);
		double xb = pb(0);
		double xc = pc(0);
		double zd = pd(2);
		double xd = pd(0);
		double yd = pd(1);

		//cout << yc << endl;
		// For a
		double beta_a = -yc * zd;
		double gamma_a = (xc * zd) - (xb * zd);
		double delta_a = yc * xd - xc * yd + xb * yd - xb * yc;

		// For b
		double beta_b = yc * zd;
		double gamma_b = -xc * zd;
		double delta_b = -yc * xd + xc * yd;

		// For c
		double beta_c = 0;
		double gamma_c = zd * xb;
		double delta_c = -yd * xb;

		// For d
		double beta_d = 0;
		double gamma_d = 0;
		double delta_d = xb * yc;

		// Fill the strain-displacement matrix Be

		MatrixXd Be(6, 12);
		Be.setZero();
		Be << beta_a, 0, 0, beta_b, 0, 0, beta_c, 0, 0, beta_d, 0, 0,
			0, gamma_a, 0, 0, gamma_b, 0, 0, gamma_c, 0, 0, gamma_d, 0,
			0, 0, delta_a, 0, 0, delta_b, 0, 0, delta_c, 0, 0, delta_d,
			gamma_a, beta_a, 0, gamma_b, beta_b, 0, gamma_c, beta_c, 0, gamma_d, beta_d, 0,
			0, delta_a, gamma_a, 0, delta_b, gamma_b, 0, delta_c, gamma_c, 0, delta_d, gamma_d,
			delta_a, 0, beta_a, delta_b, 0, beta_b, delta_c, 0, beta_c, delta_d, 0, beta_d;

		Be = 1.0 / 6.0 / volume(itet) * Be;
		cout << "StrainDisplacement Be:" << endl << Be << endl << endl;
		// cout << Be << endl;
		MatrixXd Ke(12, 6);

		Ke = Be.transpose() * D * Be;
		cout << "Stiffness Ke:" << endl << Ke << endl << endl;


		K.block(3 * a, 3 * a, 3, 3) += Ke.block(0, 0, 3, 3);
		K.block(3 * b, 3 * b, 3, 3) += Ke.block(3, 3, 3, 3);
		K.block(3 * c, 3 * c, 3, 3) += Ke.block(6, 6, 3, 3);
		K.block(3 * d, 3 * d, 3, 3) += Ke.block(9, 9, 3, 3);
		K.block(3 * a, 3 * b, 3, 3) += Ke.block(0, 3, 3, 3);
		K.block(3 * a, 3 * c, 3, 3) += Ke.block(0, 6, 3, 3);
		K.block(3 * a, 3 * d, 3, 3) += Ke.block(0, 9, 3, 3);
		K.block(3 * b, 3 * a, 3, 3) += Ke.block(3, 0, 3, 3);
		K.block(3 * b, 3 * c, 3, 3) += Ke.block(3, 6, 3, 3);
		K.block(3 * b, 3 * d, 3, 3) += Ke.block(3, 9, 3, 3);
		K.block(3 * c, 3 * a, 3, 3) += Ke.block(6, 0, 3, 3);
		K.block(3 * c, 3 * b, 3, 3) += Ke.block(6, 3, 3, 3);
		K.block(3 * c, 3 * d, 3, 3) += Ke.block(6, 9, 3, 3);
		K.block(3 * d, 3 * a, 3, 3) += Ke.block(9, 0, 3, 3);
		K.block(3 * d, 3 * b, 3, 3) += Ke.block(9, 3, 3, 3);
		K.block(3 * d, 3 * c, 3, 3) += Ke.block(9, 6, 3, 3);

		// Compute force Fe
		
		VectorXd XXX(12);
		XXX.setZero();
		XXX.segment<3>(0) = particles[a]->x0;
		XXX.segment<3>(3) = particles[b]->x0;
		XXX.segment<3>(6) = particles[c]->x0;
		XXX.segment<3>(9) = particles[d]->x0;

		VectorXd xxx(12);
		xxx.setZero();
		xxx.segment<3>(0) = particles[a]->x;
		xxx.segment<3>(3) = particles[b]->x;
		xxx.segment<3>(6) = particles[c]->x;
		xxx.segment<3>(9) = particles[d]->x;

		VectorXd Fe = Re * Ke * displaced;

		cout << "Forces:" << endl << Fe << endl << endl;
		// Fill the forces fe to vector F
		//F.segment<3>((3 * a)) -= Fe.segment<3>(3) + Fe.segment<3>(6) + Fe.segment<3>(9);
		F.segment<3>((3 * a)) += Fe.segment<3>(0);//a = 0
		
		F.segment<3>((3 * b)) += Fe.segment<3>(3);//b=3
		F.segment<3>((3 * c)) += Fe.segment<3>(6);//c=2
		F.segment<3>((3 * d)) += Fe.segment<3>(9);//d=1
	}
	//cout << "Forces "<< F << endl << endl;

	MatrixXd LHS(n, n);
	LHS.setZero();

	double damping = 0.9;

	LHS = M + h*damping*M - h*h*damping*K;

	VectorXd RHS(n);
	RHS.setZero();
	RHS = M * v + h * F ;

	VectorXd result = (LHS).ldlt().solve(RHS);

	VectorXd v_new = result;
	//cout << "velocity "<< endl<< v_new << endl<<endl;

	// Update velocity
	for (int i = 0; i < particles.size(); i++) {
		if (particles[i]->i != -1) {
			particles[i]->v = v_new.segment<3>(3 * particles[i]->i);
			cout << "v_" << particles[i]->i << endl << particles[i]->v << endl;
		}
	}

	// Update position
	for (int i = 0; i < particles.size(); i++) {
		if (particles[i]->i != -1) {
			particles[i]->x += particles[i]->v * h;
			cout << "x_" << i << endl << particles[i]->x << endl;
		}
	}

	//Collision Detection with the floor
	for (int i = 0; i < particles.size(); i++) {
		if (particles[i]->x(1) <= -2 && particles[i]->v(1) < -0.001) {
			particles[i]->x(1) = -1.99;
			particles[i]->v(1) = 0.0;
		}
	}
	updatePosNor();
}

void FemNesme::init() {
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

void FemNesme::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> p) const
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

void FemNesme::tare() {
	for (int k = 0; k < (int)particles.size(); k++) {
		particles[k]->tare();
	}
}

void FemNesme::reset() {
	for (int k = 0; k < (int)particles.size(); k++) {
		particles[k]->reset();
	}
	updatePosNor();
}

void FemNesme::updatePosNor()
{
	// Position
	for (int i = 0; i < (int)particles.size(); i++) {
		Vector3d x = particles[i]->x;
		posBuf[3 * i + 0] = x(0);
		posBuf[3 * i + 1] = x(1);
		posBuf[3 * i + 2] = x(2);
	}

	// Normal
	for (int iface = 0; iface < nTriFaces; iface++) {
		Vector3d p1 = particles[out.trifacelist[3 * iface]]->x;
		Vector3d p2 = particles[out.trifacelist[3 * iface + 1]]->x;
		Vector3d p3 = particles[out.trifacelist[3 * iface + 2]]->x;

		Vector3d e1 = p2 - p1;
		Vector3d e2 = p3 - p1;
		Vector3d normal = e1.cross(e2);

		particles[out.trifacelist[3 * iface]]->normal += normal;
		particles[out.trifacelist[3 * iface + 1]]->normal += normal;
		particles[out.trifacelist[3 * iface + 2]]->normal += normal;
	}

	for (int ipt = 0; ipt < particles.size(); ipt++) {
		Vector3d nor = particles[ipt]->normal;
		nor.normalize();
		norBuf[3 * ipt + 0] = nor(0);
		norBuf[3 * ipt + 1] = nor(1);
		norBuf[3 * ipt + 2] = nor(2);
	}

}

MatrixXd FemNesme::computeD(double lambda, double mu) {
	double v = lambda;
	double v1 = 1 + 2.0 * mu;
	double v2 = mu;

	MatrixXd D(6, 6);
	D <<
		v1, v, v, 0, 0, 0,
		v, v1, v, 0, 0, 0,
		v, v, v1, 0, 0, 0,
		0, 0, 0, v2, 0, 0,
		0, 0, 0, 0, v2, 0,
		0, 0, 0, 0, 0, v2;
	return D;
}

FemNesme::~FemNesme()
{
}
