#define TETLIBRARY

#include <iostream>
#include <tetgen.h>
#include <cmath>        // std::abs
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <omp.h>

#include "FemNesme.h"
#include "Particle.h"
#include "Program.h"
#include "GLSL.h"
#include "MatrixStack.h"
#include "Spring.h"
#include "stdlib.h"

using namespace std;
using namespace Eigen;

FemNesme::FemNesme(
	double density,
	const Eigen::Vector2d &damping)
{
	assert(density > 0.0);

	this->damping = damping;
	this->young = 1e2;
	this->poisson = 0.30;

	// Compute Lame coefficients: [mu, lambda]
	this->mu = young / (2.0 * (1.0 + poisson));
	this->lambda = young * poisson / ((1.0 + poisson)*(1.0 - 2.0 * poisson));

	double r = 0.02; // Used for collisions

	//in_2.load_ply("icosahedron"); // YOUNG: 1E1 POSSION 0.35
	//in_2.load_ply("dodecahedron"); // YOUNG: 1E0 POSSION 0.3
	//in_2.load_off("octtorus");
	//in_2.load_off("N");
	//in_2.load_ply("bunny");
	in_2.load_ply("tetrahedron"); // YOUNG: 1E2 POSSION 0.30

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

	D.resize(6, 6);

	mass.resize(nVerts);
	mass.setZero();

	volume.resize(nTets);
	volume.setZero();

	F.resize(n);
	X0.resize(n);
	x.resize(n);
	
	rotatedRestPosMatrix.resize(4 * nTets, 3);

	for (int itet = 0; itet < nTets; itet++) {
		Vector3d xa, xb, xc, xd;
		int a = out.tetrahedronlist[itet * 4];
		int b = out.tetrahedronlist[itet * 4 + 1];
		int c = out.tetrahedronlist[itet * 4 + 2];
		int d = out.tetrahedronlist[itet * 4 + 3];

		int ia = 3 * a;
		int ib = 3 * b;
		int ic = 3 * c;
		int id = 3 * d;

		xa << out.pointlist[ia], out.pointlist[ia + 1], out.pointlist[ia + 2];
		xb << out.pointlist[ib], out.pointlist[ib + 1], out.pointlist[ib + 2];
		xc << out.pointlist[ic], out.pointlist[ic + 1], out.pointlist[ic + 2];
		xd << out.pointlist[id], out.pointlist[id + 1], out.pointlist[id + 2];

		double vol = abs((xa - xd).transpose()*(xb - xd).cross(xc - xd)) / 6.0;
		volume(itet) = vol;

		double m = vol * density;

		mass(a) += m / 4;
		mass(b) += m / 4;
		mass(c) += m / 4;
		mass(d) += m / 4;
	}

	// Create particles
	for (int i = 0; i < nVerts; i++) {
		auto p = make_shared<Particle>();
		particles.push_back(p);
		p->r = r;
		p->x0 << out.pointlist[3 * i], out.pointlist[3 * i + 1], out.pointlist[3 * i + 2];
		X0.segment<3>(3 * i) = p->x0;
		p->x = p->x0;
		p->v0 << 0.0, 0.0, 0.0;
		p->v = p->v0;
		p->m = mass(i);
		p->i = i;
		p->fixed = false;
	}

	computeD(D, lambda, mu);

	for (int itet = 0; itet < nTets; itet++) {

		int a = out.tetrahedronlist[itet * 4 + 0];
		int b = out.tetrahedronlist[itet * 4 + 1];
		int c = out.tetrahedronlist[itet * 4 + 2];
		int d = out.tetrahedronlist[itet * 4 + 3];

		Matrix3d R_rest;
		computeQR(R_rest, X0, a, b, c);
		cout << R_rest << endl;

		MatrixXd rotatedRestPos(4, 3);
		rotatedRestPos.row(0).setZero(); //a
		rotatedRestPos.row(1) = R_rest * X0.segment<3>(3*b) - R_rest * X0.segment<3>(3*a);
		rotatedRestPos.row(2) = R_rest * X0.segment<3>(3*c) - R_rest * X0.segment<3>(3*a);
		rotatedRestPos.row(3) = R_rest * X0.segment<3>(3*d) - R_rest * X0.segment<3>(3*a);

		rotatedRestPosMatrix.block<4, 3>(4 * itet, 0) = rotatedRestPos;

		cout << rotatedRestPos << endl;

		//MatrixXd B_rest(6, 12);
		//computeB(B_rest, volume, itet, rotatedRestPos.row(0), rotatedRestPos.row(1), rotatedRestPos.row(2), rotatedRestPos.row(3));

		Matrix3d Rot;
		Rot.setIdentity();

		//MatrixXd BtDB(12,12);
		//computeK(BtDB, B_rest, Rot);
	
	}
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

void FemNesme::computeQR(Matrix3d &rotation, const VectorXd &x, const int &a, const int &b, const int &c) {
	Vector3d xa = x.segment<3>(3 * a);
	Vector3d xb = x.segment<3>(3 * b);
	Vector3d xc = x.segment<3>(3 * c);
	Vector3d axis_x = xb - xa;
	axis_x.normalize();

	Vector3d axis_y = xc - xa;
	axis_y.normalize();

	Vector3d axis_z = axis_x.cross(axis_y);
	axis_z.normalize();

	axis_y = axis_z.cross(axis_x);
	axis_y.normalize();

	rotation.row(0) = axis_x;
	rotation.row(1) = axis_y;
	rotation.row(2) = axis_z;
}

void FemNesme::computeB(MatrixXd &Bmatrix, VectorXd volume, int itet, Vector3d rotated_a, Vector3d rotated_b, Vector3d rotated_c, Vector3d rotated_d) {
	// Compute the coefficients of the shape functions N = alpha + beta *x + gamma * y + delta * z (for B, we don't need alpha)
	// In the local frame!!!
	double yc = rotated_c(1);
	double xb = rotated_b(0);
	double xc = rotated_c(0);
	double zd = rotated_d(2);
	double xd = rotated_d(0);
	double yd = rotated_d(1);

	// For a
	double beta_a = - yc * zd;
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

	Bmatrix.resize(6, 12);
	Bmatrix.setZero();
	Bmatrix << beta_a,       0,       0,  beta_b,       0,       0,  beta_c,       0,       0,  beta_d,       0,       0,
			        0, gamma_a,       0,       0, gamma_b,       0,       0, gamma_c,       0,       0, gamma_d,       0,
		            0,       0, delta_a,       0,       0, delta_b,       0,       0, delta_c,       0,       0, delta_d,
		      gamma_a,  beta_a,       0, gamma_b,  beta_b,       0, gamma_c,  beta_c,       0, gamma_d,  beta_d,       0,
		            0, delta_a, gamma_a,       0, delta_b, gamma_b,       0, delta_c, gamma_c,       0, delta_d, gamma_d,
		      delta_a,       0,  beta_a, delta_b,       0,  beta_b, delta_c,       0,  beta_c, delta_d,       0,  beta_d;

	Bmatrix = 1.0 / 6.0 / volume(itet) * Bmatrix;
}

void FemNesme::computeK(MatrixXd &Kmatrix, MatrixXd &KRtmatrix, MatrixXd Bmatrix, Matrix3d Rot) {
	MatrixXd Bt = Bmatrix.transpose();
	MatrixXd BtDB = Bt * D * Bmatrix;
	MatrixXd RR(12, 12);
	MatrixXd RRt(12, 12);

	RR.setZero();
	for (int i = 0; i < 4; i++) {
		RR.block<3, 3>(3 * i, 3 * i) = Rot;
	}
	RRt = RR.transpose();
	Kmatrix = RR*BtDB ;
	KRtmatrix = Kmatrix *RRt;
}

void FemNesme::step(double h, const Vector3d &grav) {
	M.setZero();
	v.setZero();
	F.setZero();
	K.setZero();
	x.setZero();

	for (int i = 0; i < particles.size(); i++) {
		int idx = particles[i]->i;
		double mass = particles[i]->m;

		Matrix3d A;
		A.setIdentity();
		A *= mass;
		M.block<3, 3>(3 * idx, 3 * idx) = A; 
		v.segment<3>(3 * idx) = particles[i]->v;
		x.segment<3>(3 * idx) = particles[i]->x;
		F.segment<3>(3 * idx) = mass * grav; 
	}

	for (int itet = 0; itet < nTets; itet++) {
		int a = out.tetrahedronlist[itet * 4 + 0];
		int b = out.tetrahedronlist[itet * 4 + 1];
		int c = out.tetrahedronlist[itet * 4 + 2];
		int d = out.tetrahedronlist[itet * 4 + 3];

		Matrix3d R_def;
		computeQR(R_def, x, a, b, c);
		cout << "R_def: " << R_def << endl << endl;

		MatrixXd rotatedCurrentPos(4, 3);
		rotatedCurrentPos.row(0).setZero(); //a
		rotatedCurrentPos.row(1) = R_def * x.segment<3>(3 * b) - R_def * x.segment<3>(3 * a);
		rotatedCurrentPos.row(2) = R_def * x.segment<3>(3 * c) - R_def * x.segment<3>(3 * a);
		rotatedCurrentPos.row(3) = R_def * x.segment<3>(3 * d) - R_def * x.segment<3>(3 * a);

		cout << "rotatedCurrentPos: " << rotatedCurrentPos << endl << endl;

		MatrixXd B_def(6, 12);
		computeB(B_def, volume, itet, rotatedCurrentPos.row(0), rotatedCurrentPos.row(1), rotatedCurrentPos.row(2), rotatedCurrentPos.row(3));

		MatrixXd RtBtDB, RtBtDBR;
		computeK(RtBtDB, RtBtDBR, B_def, R_def.transpose());
		

		VectorXd disp(12);
		disp.setZero();
	
		disp.segment<3>(3) = rotatedCurrentPos.row(1) - rotatedRestPosMatrix.row(4 * itet + 1);
		disp.segment<3>(6) = rotatedCurrentPos.row(2) - rotatedRestPosMatrix.row(4 * itet + 2);
		disp.segment<3>(9) = rotatedCurrentPos.row(3) - rotatedRestPosMatrix.row(4 * itet + 3);

		MatrixXd RR(12, 12);
		MatrixXd RRt(12, 12);

		RR.setZero();
		for (int i = 0; i < 4; i++) {
			RR.block<3, 3>(3 * i, 3 * i) = R_def;
		}
		RRt = RR.transpose();
		VectorXd Fe =  B_def.transpose() * D * B_def * disp;
		//VectorXd Fe = RtBtDB * disp;

		//cout << "Displacement disp: " << endl << disp << endl << endl;
		//cout << "StrainDisplacement Be: " << endl << B_def << endl << endl;	
		//cout << "Stiffness Ke: " << endl << RtBtDB << endl << endl;
		//cout << "Forces: " << endl << Fe << endl << endl;
		
		//F.segment<3>((3 * a)) -= Fe.segment<3>(3) + Fe.segment<3>(6) + Fe.segment<3>(9);
		F.segment<3>((3 * a)) += -R_def.transpose()*Fe.segment<3>(0);//a = 0
		F.segment<3>((3 * b)) += -R_def.transpose()*Fe.segment<3>(3);//b=3
		F.segment<3>((3 * c)) += -R_def.transpose()*Fe.segment<3>(6);//c=2
		F.segment<3>((3 * d)) += -R_def.transpose()*Fe.segment<3>(9);//d=1
	}
	
	//cout << "Forces "<< F << endl << endl;

	MatrixXd LHS(n, n);
	LHS.setZero();

	double damping = 0.9;
	LHS = M + h*damping*M;
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
			//cout << "v_" << particles[i]->i << endl << particles[i]->v << endl;
		}
	}

	// Update position
	for (int i = 0; i < particles.size(); i++) {
		if (particles[i]->i != -1) {
			particles[i]->x += particles[i]->v * h;
		}
	}

	//Collision Detection with the floor
	for (int i = 0; i < particles.size(); i++) {
		if (particles[i]->x(1) <= -2.0 && particles[i]->v(1) < -0.00001) {
			particles[i]->x(1) = -2.0;
			particles[i]->v(1) = 0.0;
		}
		//cout << "x_" << i << endl << particles[i]->x << endl;
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
	for (int iface = 0; iface < nTriFaces; iface++) {
		Vector3d p1 = particles[out.trifacelist[3 * iface + 0]]->x;
		Vector3d p2 = particles[out.trifacelist[3 * iface + 1]]->x;
		Vector3d p3 = particles[out.trifacelist[3 * iface + 2]]->x;

		//Position
		Vector3d e1 = p2 - p1;
		Vector3d e2 = p3 - p1;
		Vector3d normal = e1.cross(e2);
		normal.normalize();

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

void FemNesme::computeD(MatrixXd &Dmatrix, double lambda, double mu) {
	double v = lambda;
	double v1 = 1 + 2.0 * mu;
	double v2 = mu;

	Dmatrix <<
		v1, v, v, 0, 0, 0,
		v, v1, v, 0, 0, 0,
		v, v, v1, 0, 0, 0,
		0, 0, 0, v2, 0, 0,
		0, 0, 0, 0, v2, 0,
		0, 0, 0, 0, 0, v2;
}

FemNesme::~FemNesme()
{
}
