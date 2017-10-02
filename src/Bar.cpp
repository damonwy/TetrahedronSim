#define TETLIBRARY

#include <iostream>
#include <tetgen.h>
#include <cmath>        // std::abs

#include "Bar.h"
#include "Particle.h"
#include "Program.h"
#include "GLSL.h"
#include "MatrixStack.h"
#include "Spring.h"
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

using namespace std;
using namespace Eigen;

Bar::Bar(
	const Eigen::Vector3d &x0,
	const Eigen::Vector3d &x1,
	const Eigen::Vector3d &x2,
	const Eigen::Vector3d &x3,
	double density,
	double height,
	const Eigen::Vector2d &damping)
{
	assert(density > 0.0);
	assert(height > 0.0);
	
	this->damping = damping;
	this->young = 1e5;
	this->poisson = 0.4;
	
	double r = 0.02; // Used for collisions
	
	//tetgenio in, out;
	tetgenio::facet *f;
	tetgenio::polygon *p;
	int i;

	// All indices start from 1.
	in.firstnumber = 1;
	in.numberofpoints = 8;
	in.pointlist = new REAL[in.numberofpoints * 3];
	in.pointlist[0] = 0;  // node 1.
	in.pointlist[1] = 0;
	in.pointlist[2] = 0;
	in.pointlist[3] = 2;  // node 2.
	in.pointlist[4] = 0;
	in.pointlist[5] = 0;
	in.pointlist[6] = 2;  // node 3.
	in.pointlist[7] = 2;
	in.pointlist[8] = 0;
	in.pointlist[9] = 0;  // node 4.
	in.pointlist[10] = 2;
	in.pointlist[11] = 0;
	// Set node 5, 6, 7, 8.
	for (i = 4; i < 8; i++) {
		in.pointlist[i * 3] = in.pointlist[(i - 4) * 3];
		in.pointlist[i * 3 + 1] = in.pointlist[(i - 4) * 3 + 1];
		in.pointlist[i * 3 + 2] = height;
	}

	in.numberoffacets = 6;
	in.facetlist = new tetgenio::facet[in.numberoffacets];
	in.facetmarkerlist = new int[in.numberoffacets];

	// Facet 1. The leftmost facet.
	f = &in.facetlist[0];
	f->numberofpolygons = 1;
	f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
	f->numberofholes = 0;
	f->holelist = NULL;
	p = &f->polygonlist[0];
	p->numberofvertices = 4;
	p->vertexlist = new int[p->numberofvertices];
	p->vertexlist[0] = 1;
	p->vertexlist[1] = 2;
	p->vertexlist[2] = 3;
	p->vertexlist[3] = 4;

	// Facet 2. The rightmost facet.
	f = &in.facetlist[1];
	f->numberofpolygons = 1;
	f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
	f->numberofholes = 0;
	f->holelist = NULL;
	p = &f->polygonlist[0];
	p->numberofvertices = 4;
	p->vertexlist = new int[p->numberofvertices];
	p->vertexlist[0] = 5;
	p->vertexlist[1] = 6;
	p->vertexlist[2] = 7;
	p->vertexlist[3] = 8;

	// Facet 3. The bottom facet.
	f = &in.facetlist[2];
	f->numberofpolygons = 1;
	f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
	f->numberofholes = 0;
	f->holelist = NULL;
	p = &f->polygonlist[0];
	p->numberofvertices = 4;
	p->vertexlist = new int[p->numberofvertices];
	p->vertexlist[0] = 1;
	p->vertexlist[1] = 5;
	p->vertexlist[2] = 6;
	p->vertexlist[3] = 2;

	// Facet 4. The back facet.
	f = &in.facetlist[3];
	f->numberofpolygons = 1;
	f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
	f->numberofholes = 0;
	f->holelist = NULL;
	p = &f->polygonlist[0];
	p->numberofvertices = 4;
	p->vertexlist = new int[p->numberofvertices];
	p->vertexlist[0] = 2;
	p->vertexlist[1] = 6;
	p->vertexlist[2] = 7;
	p->vertexlist[3] = 3;

	// Facet 5. The top facet.
	f = &in.facetlist[4];
	f->numberofpolygons = 1;
	f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
	f->numberofholes = 0;
	f->holelist = NULL;
	p = &f->polygonlist[0];
	p->numberofvertices = 4;
	p->vertexlist = new int[p->numberofvertices];
	p->vertexlist[0] = 3;
	p->vertexlist[1] = 7;
	p->vertexlist[2] = 8;
	p->vertexlist[3] = 4;

	// Facet 6. The front facet.
	f = &in.facetlist[5];
	f->numberofpolygons = 1;
	f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
	f->numberofholes = 0;
	f->holelist = NULL;
	p = &f->polygonlist[0];
	p->numberofvertices = 4;
	p->vertexlist = new int[p->numberofvertices];
	p->vertexlist[0] = 4;
	p->vertexlist[1] = 8;
	p->vertexlist[2] = 5;
	p->vertexlist[3] = 1;

	// Set 'in.facetmarkerlist'

	in.facetmarkerlist[0] = -1;
	in.facetmarkerlist[1] = -2;
	in.facetmarkerlist[2] = 0;
	in.facetmarkerlist[3] = 0;
	in.facetmarkerlist[4] = 0;
	in.facetmarkerlist[5] = 0;

	// Output the PLC to files 'barin.node' and 'barin.poly'.
	in.save_nodes("barin");
	in.save_poly("barin");

	tetrahedralize("pqz", &in, &out);

	// Output mesh to files 'barout.node', 'barout.ele' and 'barout.face'.
	out.save_nodes("barout");
	out.save_elements("barout");
	out.save_faces("barout");


	nVerts = out.numberofpoints;
	nTets = out.numberoftetrahedra;
	nFacets = out.numberoffacets;
	nTriFaces = out.numberoftrifaces;

	// Build system matrices and vectors
	n = (int)3 * nVerts;
	
	M.resize(n, n);
	K.resize(n, n);
	v.resize(n);
	mass.resize(nVerts);
	mass.setZero();
	X_invs.resize((int)3 * nTets, 3);
	F.resize(n);
	X.resize(n);

	// Compute tet mass and distribute to vertices
	std::cout << "Start computing mass and X_inv..." << endl;

	for (int itet = 0; itet < nTets; itet++) {
		Vector3d xa, xb, xc, xd;
		int a = out.tetrahedronlist[itet * 4];
		int b = out.tetrahedronlist[itet * 4 + 1];
		int c = out.tetrahedronlist[itet * 4 + 2];
		int d = out.tetrahedronlist[itet * 4 + 3];

		int ia = (int)3 * a;
		xa << out.pointlist[ia], out.pointlist[ia + 1], out.pointlist[ia + 2];

		int ib = (int)3 * b;
		xb << out.pointlist[ib], out.pointlist[ib + 1], out.pointlist[ib + 2];

		int ic = (int)3 * c;
		xc << out.pointlist[ic], out.pointlist[ic + 1], out.pointlist[ic + 2];

		int id = (int)3 * d;
		xd << out.pointlist[id], out.pointlist[id + 1], out.pointlist[id + 2];

		// Compute volume and mass of each tet
		double vol = abs((xa - xd).transpose()*(xb - xd).cross(xc - xd)) / 6;
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
		X_invs.block((int)itet * 3, 0, 3, 3) = _X_inv.inverse();

	}
	std::cout << "Finish computing mass and X_inv!" << endl;
	

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
	std::cout << "Finish creating particles!" << endl;

	// Precompute Ke for each tetrahedron
	std::cout << "Start computing Kes ... " << endl;
	Kes.resize((int)12 * nTets, 12);
	Kes.setZero();

	MatrixXd E = hooke(young, poisson);
	Matrix3d E1, E2;
	E1 = E.block<3, 3>(0, 0);
	E2 = E.block<3, 3>(3, 3);

	for (int itet = 0; itet < nTets; itet++) {

		MatrixXd ys(4, 3); // Four auxiliary vectors y0, y1, y2, y3 
		Vector3d y0, y1, y2, y3;

		y1 = X_invs.row((int)itet * 3).transpose();
		y2 = X_invs.row((int)itet * 3 + 1).transpose();
		y3 = X_invs.row((int)itet * 3 + 2).transpose();
		y0 = -y1 - y2 - y3;

		ys.row(0) = y0;
		ys.row(1) = y1;
		ys.row(2) = y2;
		ys.row(3) = y3;

		MatrixXd ke(12, 12);
		ke.setZero();

		for (int i = 0; i < 4; i++) {

			// Compute vertex i 
			int ii = 3 * i; // Local starting index into the 12x12 Ke matrix
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
			
			for (int j = i; j < 4; j++) {

				//Compute vertex j
				int jj = 3 * j; // Local starting index into the 12x12 Ke matrix
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

				ke.block<3, 3>(ii, jj) -= Kij;
				if (i != j) {
					ke.block<3, 3>(jj, ii) -= Kij.transpose();
				}
			}
		}
		Kes.block((int)12 * itet, 0, 12, 12) = ke;
		//cout << ke << endl;
	}

	std::cout << "Finish computing Kes!" << endl;


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

MatrixXd Bar::hooke(double young, double poisson) {
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


void Bar::step(double h, const Vector3d &grav) {
	M.setZero();
	K.setZero();
	v.setZero();
	F.setZero();
	X.setZero();
	
	for (int i = 0; i < (int)particles.size(); i++) {
		int idx = particles[i]->i;
		double mass = particles[i]->m;

		Matrix3d A;
		A << mass, 0, 0,
			0, mass, 0,
			0, 0, mass;
		M.block<3, 3>(3 * idx, 3 * idx) = A; // filling M

		v.segment<3>(3 * idx) = particles[i]->v; // filling v 
		F.segment<3>(3 * idx) = mass * grav; // filling f with fg
		X.segment<3>(3 * idx) = particles[i]->x;// filling X with x
	}

	for (int itet = 0; itet < nTets; itet++) {
		int a = out.tetrahedronlist[itet * 4];
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
		// Compute Deformation Gradient
		Matrix3d P = dp * (X_invs.block((int)3*itet, 0, 3, 3));
		
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
		
		//Ke Element assembly
		MatrixXd RKR(12, 12);
		MatrixXd Ke(12, 12);
		Ke = Kes.block((int)itet * 12, 0, 12, 12);
		//cout << Ke << endl;
		RKR = Re * Ke* (Re.transpose());
		//cout << RKR << endl;

		K.block(3 * a, 3 * a, 3, 3) += RKR.block(0, 0, 3, 3);
		K.block(3 * b, 3 * b, 3, 3) += RKR.block(3, 3, 3, 3);
		K.block(3 * c, 3 * c, 3, 3) += RKR.block(6, 6, 3, 3);
		K.block(3 * d, 3 * d, 3, 3) += RKR.block(9, 9, 3, 3);
		K.block(3 * a, 3 * b, 3, 3) += RKR.block(0, 3, 3, 3);
		K.block(3 * a, 3 * c, 3, 3) += RKR.block(0, 6, 3, 3);
		K.block(3 * a, 3 * d, 3, 3) += RKR.block(0, 9, 3, 3);
		K.block(3 * b, 3 * a, 3, 3) += RKR.block(3, 0, 3, 3);
		K.block(3 * b, 3 * c, 3, 3) += RKR.block(3, 6, 3, 3);
		K.block(3 * b, 3 * d, 3, 3) += RKR.block(3, 9, 3, 3);
		K.block(3 * c, 3 * a, 3, 3) += RKR.block(6, 0, 3, 3);
		K.block(3 * c, 3 * b, 3, 3) += RKR.block(6, 3, 3, 3);
		K.block(3 * c, 3 * d, 3, 3) += RKR.block(6, 9, 3, 3);
		K.block(3 * d, 3 * a, 3, 3) += RKR.block(9, 0, 3, 3);
		K.block(3 * d, 3 * b, 3, 3) += RKR.block(9, 3, 3, 3);
		K.block(3 * d, 3 * c, 3, 3) += RKR.block(9, 6, 3, 3);

		
		// Fill F with f0
		VectorXd xx(12);
		xx.segment<3>(0) = particles[a]->x0;
		xx.segment<3>(3) = particles[b]->x0;
		xx.segment<3>(6) = particles[c]->x0;
		xx.segment<3>(9) = particles[d]->x0;

		
		VectorXd RKX(12);
		RKX = Re * Ke * xx;
		
		//cout << RKX << endl<<endl;
		
		F.segment<3>((int)(3 * a)) -= RKX.segment<3>(0);
		F.segment<3>((int)(3 * b)) -= RKX.segment<3>(3);
		F.segment<3>((int)(3 * c)) -= RKX.segment<3>(6);
		F.segment<3>((int)(3 * d)) -= RKX.segment<3>(9);

	}
	cout << K*X << endl;
	F += K * X;
	cout << F << endl;
	//cout << F << endl;
	MatrixXd LHS(n, n);
	LHS.setZero();

	double damping = 0.9;
	LHS = M + h*damping*M - h*h*damping*K;
	
	VectorXd RHS(n);
	RHS.setZero();
	RHS = M*v + h*F;

	VectorXd result = (LHS).ldlt().solve(RHS);

	VectorXd v_new = result;
	//cout << v_new << endl << endl;
	//std::cout << v_new << endl << endl;
	// Compute the strain
	//Matrix3d strain_m = 0.5 * (du + du.transpose() + du.transpose() * du);

	//VectorXd strain_v(6);
	//strain_v << strain_m(0, 0), strain_m(1, 1), strain_m(2, 2),
	//	strain_m(0, 1), strain_m(1, 2), strain_m(0, 2);

	// Compute the stress
	//VectorXd stress_v = E * strain_v;

	/*Matrix3d stress_m;
	stress_m <<
		stress_v(0), stress_v(3), stress_v(4),
		stress_v(3), stress_v(1), stress_v(5),
		stress_v(4), stress_v(5), stress_v(2);

	MatrixXi triIndices(4, 3);
	triIndices <<
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
	//VectorXi fixed(0);
	//fixed << 0;
	// Fixed point constraints:
	//MatrixXd G(fixed.size()*3, 12);
	//G.setZero();
	//for (int i = 0; i < fixed.size(); i++) {
	//	G.block<3, 3>(3 * i, 3* fixed(i)) = I; // filling G
	//}
	//int numfixed = fixed.size();
	//MatrixXd LHS(12+(int)(numfixed*3), 12+(int)(numfixed*3));
	//LHS.block<12, 12>(0, 0) = M + h*damping*M;
	//MatrixXd temp(12, fixed.size());
	//temp = G.transpose();
	//LHS.block(0, 12, 12, numfixed*3) = temp;
	//LHS.block(12, 0, numfixed*3, 12) = G;
	//VectorXd RHS(12+(int)(numfixed*3));
	//VectorXd v_new = (M).ldlt().solve(M*v + h*f);

	// Update velocity
	for (int i = 0; i < (int)particles.size(); i++) {
		if (particles[i]->i != -1) {
			particles[i]->v = v_new.segment<3>(3 * particles[i]->i);
		}
	}

	// Update position
	for (int i = 0; i < (int)particles.size(); i++) {
		if (particles[i]->i != -1) {
			particles[i]->x += particles[i]->v * h;
			//cout << particles[i]->x << endl;
		}
	}
	//cout << X << endl << endl;
	 //Collision Detection with the floor
	for (int i = 0; i < (int)particles.size(); i++) {
		if (particles[i]->x(1) <= -7 && particles[i]->v(1)<0) {
			particles[i]->x(1) = -7;
			particles[i]->v(1) = 0;
		}
	}

	updatePosNor();
}


void Bar::init() {
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

void Bar::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> p) const
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


	/*int h_nor = p->getAttribute("aNor");
	glEnableVertexAttribArray(h_nor);
	glBindBuffer(GL_ARRAY_BUFFER, norBufID);
	glBufferData(GL_ARRAY_BUFFER, norBuf.size() * sizeof(float), &norBuf[0], GL_DYNAMIC_DRAW);
	glVertexAttribPointer(h_nor, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);*/


	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);

	glDrawElements(GL_TRIANGLES, 3*nTriFaces, GL_UNSIGNED_INT, (const void *)(0 * sizeof(unsigned int)));

	//glDisableVertexAttribArray(h_nor);
	glDisableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	MV->popMatrix();
}

void Bar::tare() {
	for (int k = 0; k < (int)particles.size(); k++) {
		particles[k]->tare();
	}
}

void Bar::reset() {
	for (int k = 0; k < (int)particles.size(); k++) {
		particles[k]->reset();
	}
	updatePosNor();
}

void Bar::updatePosNor()
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

Bar::~Bar()
{
}
