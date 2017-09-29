#include <iostream>

#include "Tetrahedron.h"
#include "Particle.h"
#include "Program.h"
#include "GLSL.h"
#include "MatrixStack.h"
#include "Spring.h"
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

using namespace std;
using namespace Eigen;

Tetrahedron::Tetrahedron(
	const Eigen::Vector3d &x0,
	const Eigen::Vector3d &x1,
	const Eigen::Vector3d &x2,
	const Eigen::Vector3d &x3,
	double mass,
	double stiffness,
	const Eigen::Vector2d &damping)
{
	assert(mass > 0.0);
	assert(stiffness > 0.0);
	this->steps = 0;
	this->damping = damping;
	this->young = 1e5;
	this->poisson = 0.4;

	double r = 0.02; // Used for collisions
	int nVerts = 4;

	MatrixXd X(4, 3);
	X.row(0) = x0;
	X.row(1) = x1;
	X.row(2) = x2;
	X.row(3) = x3;
		 
	double mass_n = mass / (nVerts);
	n = 0;
	// Create particles
	for (int i = 0; i < nVerts; i++) {
		auto p = make_shared<Particle>();
		particles.push_back(p);
		p->r = r;
		p->x = X.row(i); // Initialize pos of nodes
		p->v << 0.0, 0.0, 0.0;
		p->m = mass_n;
		p->i = i;
		p->fixed = false;
		n+=3;
	}
	
	// Build system matrices and vectors
	M.resize(n, n);
	K.resize(n, n);
	v.resize(n);
	f.resize(n);
	p.resize(n);
	// Precompute X_inv = [x1,x2,x3]-1, constant throughout the simulation
	MatrixXd _X_inv(3, 3);
	Vector3d xb = x1 - x0;
	Vector3d xc = x2 - x0;
	Vector3d xd = x3 - x0;
	_X_inv.col(0) = xb;
	_X_inv.col(1) = xc;
	_X_inv.col(2) = xd;
	X_inv = _X_inv.inverse(); 

	// 
	// Build vertex buffers
	posBuf.clear();
	norBuf.clear();
	texBuf.clear();
	eleBuf.clear();
	posBuf.resize(nVerts * 3);
	norBuf.resize(nVerts * 3);
	updatePosNor();

	//// Texture coordinated
	
	//texBuf = { 0, 0, 0, 1, 1, 1, 1, 0 };
	// Elements
	eleBuf =  { 
		0, 2, 1,
		1, 2, 3,
		2, 0, 3,
		3, 0, 2};

}

MatrixXd Tetrahedron::hooke(double young, double poisson) {
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


void Tetrahedron::step(double h, const Vector3d &grav) {
	M.setZero();
	K.setZero();
	v.setZero();
	f.setZero();
	p.setZero();
	for (int i = 0; i < (int)particles.size(); i++) {
		int idx = particles[i]->i;
		double mass = particles[i]->m;
		
		Matrix3d A;
		A << mass, 0, 0,
			0, mass, 0,
			0, 0, mass;
		M.block<3, 3>(3*idx, 3*idx) = A; // filling M

		Vector3d B;
		B << particles[i]->v;
		v.segment<3>(3*idx) = B; // filling v 
	    f.segment<3>(3*idx) = mass * grav; // filling f with fg
		p.segment<3>(3 * idx) = particles[i]->x;
	}

	MatrixXd dp(3, 3);
	Vector3d pb = particles[1]->x - particles[0]->x;
	Vector3d pc = particles[2]->x - particles[0]->x;
	Vector3d pd = particles[3]->x - particles[0]->x;
	dp.col(0) = pb;
	dp.col(1) = pc;
	dp.col(2) = pd;
	// Compute Deformation Gradient
	Matrix3d P = dp*(X_inv);
	Matrix3d I;
	I.setIdentity();
	MatrixXd du = P - I;

	// Compute R and Re
	Vector3d a0 = P.col(0);
	Vector3d a1 = P.col(1);
	Vector3d r0 = a0;
	r0 = r0 / r0.norm();
	Vector3d r1 = a1 - r0.dot(a1)*r0;
	r1 = r1 / r1.norm();
	Vector3d r2 = r0.cross(r1);
	Matrix3d R;
	R.col(0) = r0;
	R.col(1) = r1;
	R.col(2) = r2;

	MatrixXd Re(12, 12);
	Re.setZero();
	// Filling Re
	for (int i = 0; i < 4; i++) {
		Re.block<3, 3>(3 * i, 3 * i) = R;
	}
	//cout << Re << endl;

	// Compute the strain
	Matrix3d strain_m = 0.5 * (du + du.transpose() + du.transpose() * du);

	VectorXd strain_v(6);
	strain_v << strain_m(0, 0), strain_m(1, 1), strain_m(2, 2), 
	strain_m(0, 1), strain_m(1, 2), strain_m(0, 2);
	
	MatrixXd E = hooke(young, poisson);
	Matrix3d E1, E2; 
	E1 = E.block<3, 3>(0, 0);
	E2 = E.block<3, 3>(3, 3);

	// Compute the stress
	VectorXd stress_v = E * strain_v;

	Matrix3d stress_m;
	stress_m <<
		stress_v(0), stress_v(3), stress_v(4),
		stress_v(3), stress_v(1), stress_v(5),
		stress_v(4), stress_v(5), stress_v(2);

	MatrixXi triIndices(4, 3);
	triIndices<< 
		0, 1, 2,
		0, 2, 3,
		1, 3, 2,
		0, 3, 1;

	for (int i = 0; i < 4; i++) {
		
		int ia = triIndices(i, 0);
		int ib = triIndices(i, 1);
		int ic = triIndices(i, 2);

		Vector3d pa = particles[ia]->x;
		Vector3d pb = particles[ib]->x;
		Vector3d pc = particles[ic]->x;

		Vector3d fk = stress_m * ((pb - pa).cross(pc - pa)); // f_face
		//cout << fk << endl << endl;
		// filling f with fk
		/*for (int i = 0; i < 3; i++) {
			f(3 * ia + i) += fk(i) / 3;
			f(3 * ib + i) += fk(i) / 3;
			f(3 * ic + i) += fk(i) / 3;
		}*/
	}

	// Compute stiffness matrix Ke for each element/tet

	MatrixXd Ke(12, 12);
	Ke.setZero();

	MatrixXd ys(4, 3); // Four auxiliary vectors y0, y1, y2, y3 

	Vector3d y0, y1, y2, y3;
	y1 = X_inv.row(0).transpose();
	y2 = X_inv.row(1).transpose();
	y3 = X_inv.row(2).transpose();
	y0 = -y1 - y2 - y3;

	ys.row(0) = y0;

	ys.row(1) = y1;
	ys.row(2) = y2;
	ys.row(3) = y3;

	for (int i = 0; i < 4; i++) {
		int ii = 3 * i; // Local starting index into the 12x12 Ke matrix
		// Compute vertex i 
		Vector3d yi = ys.row(i);
		Matrix3d Ni; // Diagonal Matrix
		Ni.setZero();
		for (int idx = 0; idx < 3; idx++) {
			Ni(idx, idx) = yi(idx);
		}
		Matrix3d NiEn = Ni * E1;

		Matrix3d Si;
		Si << yi(1), 0, yi(2),
			  yi(0), yi(2), 0,
			  0, yi(1), yi(0);

		Matrix3d SiEs = Si * E2;

		for (int j = 0; j < 4; j++) {

			//Compute vertex j
			int jj = 3 * j;// Local starting index into the 12x12 Ke matrix
			Vector3d yj = ys.row(j);

			Matrix3d Nj; // Diagonal Matrix
			Nj.setZero();
			for (int idx = 0; idx < 3; idx++) {
				Nj(idx, idx) = yj(idx);
			}
			
			Matrix3d Sj;
			Sj << yj(1), 0, yj(2),
				 yj(0), yj(2), 0,
				 0, yj(1), yj(0);

			Matrix3d Kij = NiEn * Nj + SiEs * (Sj.transpose());
			Ke.block<3, 3>(ii, jj) -= Kij;


			if (i != j) {
				Ke.block<3,3>(jj,ii) -= Kij.transpose();
			}

		}
	}

	//Ke Element assembly
	//cout << Ke << endl;
	MatrixXd K(12, 12);
	K = Re * Ke* (Re.transpose());
	//cout << K << endl;
	
	f += K*p - Re*Ke*p;
	//cout << f << endl;

	MatrixXd G(9, 12);
	G.setZero();
	VectorXi fixed(3);
	fixed << 0, 1, 3;
	double damping = 0.9;
	// Fixed point constraints:

	for (int i = 0; i < fixed.size(); i++) {
		G.block<3, 3>(3 * i, 3* fixed(i)) = I; // filling G
	}

	//cout << G << endl;

	MatrixXd LHS(21, 21);
	LHS.setZero();
	LHS.block<12, 12>(0, 0) = M + h*damping*M - h*h*damping*K;
	//LHS.block<12, 12>(0, 0) = M + h*damping*M;
	MatrixXd temp(12, 9);
	temp = G.transpose();
	LHS.block<12, 9>(0, 12) = temp;
	LHS.block<9, 12>(12, 0) = G;
	
	//cout << LHS << endl;

	VectorXd RHS(21);
	RHS.setZero();
	RHS.segment<12>(0) = M*v + h*f;

	//cout << f << endl << endl;

	//VectorXd v_new = (M).ldlt().solve(M*v + h*f);
	VectorXd result = (LHS).ldlt().solve(RHS);
	cout << result(7) << endl<<endl;

	VectorXd v_new = result.segment<12>(0);
	// Update velocity
	for (int i = 0; i < (int)particles.size(); i++) {
		if (particles[i]->i != -1) {
			particles[i]->v = v_new.segment<3>(3*particles[i]->i);
		}
	}

	// Update position
	for (int i = 0; i < (int)particles.size(); i++) {
		if (particles[i]->i != -1) {
			particles[i]->x += particles[i]->v *h;
		}
		//cout << particles[i]->x << endl;
	}

	//cout << "v_new" << v_new << endl << endl;
	updatePosNor();
}


void Tetrahedron::init() {
	glGenBuffers(1, &posBufID);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size() * sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);

	//glGenBuffers(1, &posBufID);
	//glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	//glBufferData(GL_ARRAY_BUFFER, posBuf.size() * sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);

	/*glGenBuffers(1, &texBufID);
	glBindBuffer(GL_ARRAY_BUFFER, texBufID);
	glBufferData(GL_ARRAY_BUFFER, texBuf.size() * sizeof(float), &texBuf[0], GL_STATIC_DRAW);
*/
	glGenBuffers(1, &eleBufID);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, eleBuf.size() * sizeof(unsigned int), &eleBuf[0], GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	assert(glGetError() == GL_NO_ERROR);

}

void Tetrahedron::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> p) const
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

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);

	glDrawElements(GL_TRIANGLES, 12, GL_UNSIGNED_INT, (const void *)(0  * sizeof(unsigned int)));
	
	glDisableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	MV->popMatrix();
}

void Tetrahedron::tare() {
	for (int k = 0; k < (int)particles.size(); k++) {
		particles[k]->tare();
	}
}

void Tetrahedron::reset() {
	for (int k = 0; k < (int)particles.size(); k++) {
		particles[k]->reset();
	}
	updatePosNor();
}

void Tetrahedron::updatePosNor()
{
	// Position
	for (int i = 0; i < (int)particles.size(); i++) {
		Vector3d x = particles[i]->x;
		posBuf[3 * i + 0] = x(0);
		posBuf[3 * i + 1] = x(1);
		posBuf[3 * i + 2] = x(2);
	}

	// Normal
	// first triangles
	//Vector3d dx0 = particles[2]->x - particles[0]->x;
	//Vector3d dx1 = particles[1]->x - particles[0]->x;
	//Vector3d c = dx0.cross(dx1);

}

Tetrahedron::~Tetrahedron()
{
}
