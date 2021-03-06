#define TETLIBRARY

#include <iostream>
#include <tetgen.h>
#include <cmath>        // std::abs
#include <math.h>
#include "ChronoTimer.h"

#include "FemMuller.h"
#include "Particle.h"
#include "Muscle.h"
#include "Program.h"
#include "GLSL.h"
#include "MatrixStack.h"
#include "Spring.h"
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

using namespace std;
using namespace Eigen;
typedef Eigen::Triplet<double> T;

FemMuller::FemMuller(
	double density,
	const Eigen::Vector2d &damping)
{
	assert(density > 0.0);
	this->damping = damping;
	this->poisson = 0.4;

	double r = 0.02;
	model = 1;
	// 0: cube 1: rabbit 

	//in_2.load_ply("icosahedron");
	//in_2.load_ply("dodecahedron");
	//in_2.load_off("octtorus");

	if (model == 0) {
		in_2.load_ply("cube");
		this->young = 1e2;
	}
	if (model == 1) {
		in_2.load_ply("bunny600");
		this->young = 1e0;
	}
	

	tetrahedralize("pq1.2Vz", &in_2, &out);

	out.save_nodes("out");
	out.save_elements("out");
	out.save_faces("out");

	nVerts = out.numberofpoints;
	nTets = out.numberoftetrahedra;
	nFacets = out.numberoffacets;
	nTriFaces = out.numberoftrifaces;

	n = 3 * nVerts;

	M.resize(n, n);
	M.setZero();

	K.resize(n, n);
	v.resize(n);

	mass.resize(nVerts);
	mass.setZero();

	volume.resize(nTets);
	volume.setZero();

	X_invs.resize(3 * nTets, 3);

	F0.resize(n);
	F0.setZero();
	F.resize(n);

	X.resize(n);
	D.resize(6, 6);

	RKR.resize(12, 12);
	Ke.resize(12, 12);
	Re.resize(12, 12);

	xx.resize(12);
	RKX.resize(12);

	// Create Muscles
	
	if (model == 0) {
		createMuscle(Vector3d(-10.0, 0.5, 0.0), Vector3d(10.0, 0.5, 0.0));
		createMuscle(Vector3d(-10.0, 0.5, 0.5), Vector3d(10.0, 0.5, 0.5));
		createMuscle(Vector3d(-10.0, 0.5, -0.5), Vector3d(10.0, 0.5, -0.5));

		createMuscle(Vector3d(-10.0, 0.0, 0.5), Vector3d(10.0, 0.0, 0.5));
		createMuscle(Vector3d(-10.0, 0.0, -0.5), Vector3d(10.0, 0.0, -0.5));

		createMuscle(Vector3d(-10.0, -0.5, -0.5), Vector3d(10.0, -0.5, -0.5));
		createMuscle(Vector3d(-10.0, -0.5, 0.5), Vector3d(10.0, -0.5, 0.5));
		createMuscle(Vector3d(-10.0, -0.5, 0.0), Vector3d(10.0, -0.5, 0.0));
		
		isTop = false;
		isBottom = false;
		isLeft = false;
		isRight = false;
	}

	if (model == 1) {
		double y = 0.4;
		double x = 10.0;

		for (int t = 0; t < 4; t++) {
			double z = -0.3;
			for (int r = 0; r < 4; r++) {
				createMuscle(Vector3d(-x, y, z), Vector3d(x, y, z));
				z += 0.2;
			}
			y += 0.15;
		}
		createMuscle(Vector3d(-10.0, 1.3, 0.0), Vector3d(10.0, 1.3, 0.0));
		createMuscle(Vector3d(-10.0, 1.3, 0.1), Vector3d(10.0, 1.3, 0.1));
		createMuscle(Vector3d(-10.0, 1.4, 0.0), Vector3d(10.0, 1.4, 0.0));
		createMuscle(Vector3d(-10.0, 1.4, 0.2), Vector3d(10.0, 1.4, 0.2));
		createMuscle(Vector3d(-10.0, 1.4, 0.4), Vector3d(10.0, 1.4, 0.4));

		createMuscle(Vector3d(-10.0, 1.1, 0.0), Vector3d(10.0, 1.1, 0.0));
		createMuscle(Vector3d(-10.0, 1.1, 0.2), Vector3d(10.0, 1.1, 0.2));
		createMuscle(Vector3d(-10.0, 1.1, 0.3), Vector3d(10.0, 1.1, 0.3));

		isTop = false;
		isBottom = false;
		isLeft = false;
		isRight = false;	
	}

	// Compute tet mass and distribute to vertices
	std::cout << "Start computing mass and X_inv..." << endl;

	for (int itet = 0; itet < nTets; itet++) {
		Vector3d xa, xb, xc, xd;
		int a = out.tetrahedronlist[itet * 4 + 0];
		int b = out.tetrahedronlist[itet * 4 + 1];
		int c = out.tetrahedronlist[itet * 4 + 2];
		int d = out.tetrahedronlist[itet * 4 + 3];
		int ia = 3 * a;
		int ib = 3 * b;
		int ic = 3 * c;
		int id = 3 * d;
		xb << out.pointlist[ib], out.pointlist[ib + 1], out.pointlist[ib + 2];
		xa << out.pointlist[ia], out.pointlist[ia + 1], out.pointlist[ia + 2];
		xc << out.pointlist[ic], out.pointlist[ic + 1], out.pointlist[ic + 2];
		xd << out.pointlist[id], out.pointlist[id + 1], out.pointlist[id + 2];

		// Compute volume and mass of each tet
		double vol = abs((xa - xd).transpose()*(xb - xd).cross(xc - xd)) / 6;
		volume(itet) = vol;

		double m = vol * density;
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

		// Compute Muscle Segments
		for (int imus = 0; imus < muscles.size(); imus++) {
			Vector3d orig = muscles[imus]->p0->x;

			double t, u, v;
			int numIntersects = 0;
			vector < shared_ptr<Particle> > nodes;
			direction = muscles[imus]->direction;

			if (rayTriangleIntersects(xa, xc, xb, direction, orig, t, u, v)) {
				numIntersects += 1;
				auto p = make_shared<Particle>();
				nodes.push_back(p);
				p->U = u;
				p->V = v;
				p->x = orig + direction * t;
				p->triIndex = Vector3i(a, c, b);
			}

			if (rayTriangleIntersects(xa, xb, xd, direction, orig, t, u, v)) {
				numIntersects += 1;
				auto p = make_shared<Particle>();
				nodes.push_back(p);
				p->U = u;
				p->V = v;
				p->x = orig + direction * t;
				p->triIndex = Vector3i(a, b, d);
			}

			if (rayTriangleIntersects(xb, xd, xc, direction, orig, t, u, v)) {
				numIntersects += 1;
				auto p = make_shared<Particle>();
				nodes.push_back(p);
				p->U = u;
				p->V = v;
				p->x = orig + direction * t;
				p->triIndex = Vector3i(b, d, c);
			}

			if (rayTriangleIntersects(xa, xc, xd, direction, orig, t, u, v)) {
				numIntersects += 1;
				auto p = make_shared<Particle>();
				nodes.push_back(p);
				p->U = u;
				p->V = v;
				p->x = orig + direction * t;
				p->triIndex = Vector3i(a, c, d);
			}

			if (numIntersects == 2) {
				auto seg = make_shared<Segment>(nodes[0], nodes[1]);
				seg->eleID = itet;
				muscles[imus]->insertSegment(seg);
			}
		}
	}

	// Init muscles
	for (int imus = 0; imus < muscles.size(); imus++) {
		muscles[imus]->init();
	}

	std::cout << "Finish computing mass and X_inv!" << endl;
	I.setIdentity();
	Vector3d grav;
	grav << 0.0, -10, 0.0;
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
		M.block<3, 3>(3 * i, 3 * i) = I * p->m;
		F0.segment<3>(3 * i) = p->m * grav;
	}
	std::cout << "Finish creating particles!" << endl;


	// Precompute Ke for each tetrahedron
	std::cout << "Start computing Kes ... " << endl;
	Kes.resize(12 * nTets, 12);
	Kes.setZero();

	MatrixXd E = hooke(young, poisson);
	D = E;
	Matrix3d E1, E2;
	E1 = E.block<3, 3>(0, 0);
	E2 = E.block<3, 3>(3, 3);

	for (int itet = 0; itet < nTets; itet++) {

		MatrixXd ys(4, 3); // Four auxiliary vectors y0, y1, y2, y3 
		Vector3d y0, y1, y2, y3;

		y1 = X_invs.row(itet * 3).transpose();
		y2 = X_invs.row(itet * 3 + 1).transpose();
		y3 = X_invs.row(itet * 3 + 2).transpose();
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
		Kes.block(12 * itet, 0, 12, 12) = ke;
	}

	std::cout << "Finish computing Kes!" << endl;

	posBuf.clear();
	norBuf.clear();
	texBuf.clear();
	eleBuf.clear();
	posBuf.resize(nTriFaces * 9);
	norBuf.resize(nTriFaces * 9);
	eleBuf.resize(nTriFaces * 3);

	updatePosNor();

	for (int i = 0; i < nTriFaces; i++) {
		for (int t = 0; t < 3; t++) {
			eleBuf[3 * i + t] = 3 * i + t;
		}
	}
}

void FemMuller::createMuscle(Eigen::Vector3d p0, Eigen::Vector3d p1) {
	auto n0 = make_shared<Particle>();
	n0->x = p0;
	auto n1 = make_shared<Particle>();
	n1->x = p1;
	auto muscle = make_shared<Muscle>(n0, n1);
	muscles.push_back(muscle);
}

void FemMuller::step(double h, const Vector3d &grav, const bool *keyToggles) {
	//ChronoTimer timer("test", 1);
	//timer.tic();

	vector<T> A_;
	K.setZero();
	v.setZero();
	X.setZero();
	F = F0;
	double damping = 0.9;
	double coef = 1 + h * damping;

	if (keyToggles[(unsigned)'t']&& model == 1) {
		if (isTop) {
			for (int tt = 0; tt < 8; tt++) {
				muscles[16+tt]->isActive = false;
			}
			
			isTop = false;
		}
		else {
			for (int tt = 0; tt < 8; tt++) {
				muscles[16 + tt]->isActive = true;
			}
			isTop = true;
		}
	}
	if (keyToggles[(unsigned)'d'] && model == 1) {
		if (isBottom) {
			for (int tt = 0; tt < 8; tt++) {
				muscles[0 + tt]->isActive = false;
			}

			isBottom = false;
		}
		else {
			for (int tt = 0; tt < 8; tt++) {
				muscles[0 + tt]->isActive = true;
			}
			isBottom = true;
		}
	}
	if (keyToggles[(unsigned)'m'] && model == 1) {
		if (isLeft) {
			for (int tt = 0; tt < 8; tt++) {
				muscles[8 + tt]->isActive = false;
			}

			isLeft = false;
		}
		else {
			for (int tt = 0; tt < 8; tt++) {
				muscles[8 + tt]->isActive = true;
			}
			isLeft = true;
		}
	}
	if (keyToggles[(unsigned)'t'] && model == 0) {
		if (isTop) {
			for (int tt = 0; tt < 3; tt++) {
				muscles[0 + tt]->isActive = false;
			}
			isTop = false;
		}
		else {
			for (int tt = 0; tt < 3; tt++) {
				muscles[0 + tt]->isActive = true;
			}
			isTop = true;
		}
	}
	if (keyToggles[(unsigned)'d'] && model == 0) {
		if (isBottom) {
			for (int tt = 0; tt < 3; tt++) {
				muscles[5 + tt]->isActive = false;
			}

			isBottom = false;
		}
		else {
			for (int tt = 0; tt < 3; tt++) {
				muscles[5 + tt]->isActive = true;
			}
			isBottom = true;
		}
	}
	
	for (int i = 0; i < particles.size(); i++) {

		v.segment<3>(3 * i) = particles[i]->v;
		X.segment<3>(3 * i) = particles[i]->x;
		for (int j = 0; j < 3; j++) {
			A_.push_back(T(3 * i + j, 3 * i + j, coef * particles[i]->m));
		}
	}

	auto mus = make_shared<Muscle>();
	auto seg = make_shared<Segment>();
	Vector3d nor_acb, f_acb, nor_abd, f_abd, nor_bdc, f_bdc, nor_acd, f_acd;
	double A_acb, A_abd, A_bdc, A_acd;

	// coef for K
	coef = -h * h * damping;

	for (int itet = 0; itet < nTets; itet++) {
		a = out.tetrahedronlist[itet * 4 + 0];
		b = out.tetrahedronlist[itet * 4 + 1];
		c = out.tetrahedronlist[itet * 4 + 2];
		d = out.tetrahedronlist[itet * 4 + 3];

		xa = particles[a]->x;
		pb = particles[b]->x - xa;
		pc = particles[c]->x - xa;
		pd = particles[d]->x - xa;

		dp.col(0) = pb;
		dp.col(1) = pc;
		dp.col(2) = pd;

		// Compute Deformation Gradient
		P = dp * (X_invs.block(3 * itet, 0, 3, 3));

		// Compute R and Re
		a0 = P.col(0);
		a1 = P.col(1);
		r0 = a0;
		r0 = r0 / r0.norm();
		r1 = a1 - r0.dot(a1) * r0;

		r1 = r1 / r1.norm();
		r2 = r0.cross(r1);

		R.col(0) = r0;
		R.col(1) = r1;
		R.col(2) = r2;
		
		Re.setZero();
		for (int i = 0; i < 4; i++) {
			Re.block<3, 3>(3 * i, 3 * i) = R;
		}

		// Compute muscle forces
		for (int imus = 0; imus < muscles.size(); imus++) {
			
			mus = muscles[imus];
			direction = mus->direction;
			U << direction(0), 0, 0, direction(1), 0, 0, direction(2), 0, 0;

			for (int iseg = 0; iseg < mus->numElements; iseg++) {

				if (itet == mus->elementIDs[iseg] && mus->isActive) {

					seg = mus->segments[iseg];
					double elongation = seg->Ld - seg->L;
					double Fseg = seg->k * elongation;

					FF << Fseg, 0, 0, 0, 0, 0, 0, 0, 0;

					stress = R * U * FF * U.transpose() * R.transpose();

					// Create Muscle Segments
					//  f0: 0 2 1  f1: 0 1 3 f2: 1 3 2 f3: 0 2 3  CCW

					// tri a c b
					nor_acb = pc.cross(particles[b]->x - particles[c]->x).normalized();
					A_acb = pc.cross(particles[b]->x - particles[c]->x).norm() * 1.0 / 2.0;
					f_acb = A_acb * stress * nor_acb / 3.0;

					F.segment<3>((3 * a)) += f_acb;
					F.segment<3>((3 * c)) += f_acb;
					F.segment<3>((3 * b)) += f_acb;

					// tri  a b d
					nor_abd = pb.cross(particles[d]->x - particles[b]->x).normalized();
					A_abd = pb.cross(particles[d]->x - particles[b]->x).norm() * 1.0 / 2.0;
					f_abd = A_abd * stress * nor_abd / 3.0;

					F.segment<3>((3 * a)) += f_abd;
					F.segment<3>((3 * b)) += f_abd;
					F.segment<3>((3 * d)) += f_abd;

					// tri b d c
					nor_bdc = (particles[d]->x - particles[b]->x).cross(particles[c]->x - particles[d]->x).normalized();
					A_bdc = (particles[d]->x - particles[b]->x).cross(particles[c]->x - particles[d]->x).norm() * 1.0 / 2.0;
					f_bdc = A_bdc * stress * nor_bdc / 3.0;

					F.segment<3>((3 * b)) -= f_bdc;
					F.segment<3>((3 * d)) -= f_bdc;
					F.segment<3>((3 * c)) -= f_bdc;

					// tri a c d
					nor_acd = pc.cross(particles[d]->x - particles[c]->x).normalized();
					A_acd = pc.cross(particles[d]->x - particles[c]->x).norm() * 1.0 / 2.0;
					f_acd = A_acd * stress * nor_acd / 3.0;

					F.segment<3>((3 * a)) -= f_acd;
					F.segment<3>((3 * c)) -= f_acd;
					F.segment<3>((3 * d)) -= f_acd;
				}
			}
		}

		//Ke Element assembly
		
		Ke = Kes.block(itet * 12, 0, 12, 12);
		RKR = Re * Ke* (Re.transpose());

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

		for (int idx = 0; idx < 3; idx++) {
			for (int jdx = 0; jdx < 3; jdx++) {
				A_.push_back(T(3 * a + idx, 3 * a + jdx, coef*RKR(0 + idx, 0 + jdx)));
				A_.push_back(T(3 * a + idx, 3 * b + jdx, coef*RKR(0 + idx, 3 + jdx)));
				A_.push_back(T(3 * a + idx, 3 * c + jdx, coef*RKR(0 + idx, 6 + jdx)));
				A_.push_back(T(3 * a + idx, 3 * d + jdx, coef*RKR(0 + idx, 9 + jdx)));

				A_.push_back(T(3 * b + idx, 3 * a + jdx, coef*RKR(3 + idx, 0 + jdx)));
				A_.push_back(T(3 * b + idx, 3 * b + jdx, coef*RKR(3 + idx, 3 + jdx)));
				A_.push_back(T(3 * b + idx, 3 * c + jdx, coef*RKR(3 + idx, 6 + jdx)));
				A_.push_back(T(3 * b + idx, 3 * d + jdx, coef*RKR(3 + idx, 9 + jdx)));

				A_.push_back(T(3 * c + idx, 3 * a + jdx, coef*RKR(6 + idx, 0 + jdx)));
				A_.push_back(T(3 * c + idx, 3 * b + jdx, coef*RKR(6 + idx, 3 + jdx)));
				A_.push_back(T(3 * c + idx, 3 * c + jdx, coef*RKR(6 + idx, 6 + jdx)));
				A_.push_back(T(3 * c + idx, 3 * d + jdx, coef*RKR(6 + idx, 9 + jdx)));

				A_.push_back(T(3 * d + idx, 3 * a + jdx, coef*RKR(9 + idx, 0 + jdx)));
				A_.push_back(T(3 * d + idx, 3 * b + jdx, coef*RKR(9 + idx, 3 + jdx)));
				A_.push_back(T(3 * d + idx, 3 * c + jdx, coef*RKR(9 + idx, 6 + jdx)));
				A_.push_back(T(3 * d + idx, 3 * d + jdx, coef*RKR(9 + idx, 9 + jdx)));
			}
		}


		//f Fill F with f0
		
		xx.segment<3>(0) = particles[a]->x0;
		xx.segment<3>(3) = particles[b]->x0;
		xx.segment<3>(6) = particles[c]->x0;
		xx.segment<3>(9) = particles[d]->x0;

		RKX = Re * Ke * xx;

		F.segment<3>((3 * a)) -= RKX.segment<3>(0);
		F.segment<3>((3 * b)) -= RKX.segment<3>(3);
		F.segment<3>((3 * c)) -= RKX.segment<3>(6);
		F.segment<3>((3 * d)) -= RKX.segment<3>(9);
	}
	
	//timer.toc();
	//timer.print();

	F += K * X;

	//MatrixXd LHS(n, n);
	//LHS.setZero();
	
	//LHS = M + h * damping * M - h * h * damping * K;

	VectorXd RHS(n);
	RHS.setZero();
	RHS = M * v + h * F;

	Eigen::SparseMatrix<double> A(n, n);
	A.setFromTriplets(A_.begin(), A_.end());
	ConjugateGradient< SparseMatrix<double> > cg;
	cg.setMaxIterations(25);
	cg.setTolerance(1e-3);
	cg.compute(A);

	//timer.toc();
	//timer.print();

	//VectorXd result = (LHS).ldlt().solve(RHS);
	VectorXd result = cg.solveWithGuess(RHS, v);
	//timer.toc();
	//timer.print();

	VectorXd v_new = result;

	// Update velocity
	for (int i = 0; i < particles.size(); i++) {
		if (particles[i]->i != -1) {
			particles[i]->v = v_new.segment<3>(3 * particles[i]->i);
		}
	}
	
	// Update position
	for (int i = 0; i <particles.size(); i++) {
		if (particles[i]->i != -1) {
			particles[i]->x += particles[i]->v * h;
		}
	}

	//timer.toc();
	//timer.print();

	//Collision Detection with the floor
	for (int i = 0; i <particles.size(); i++) {
		if (particles[i]->x(1) <= -2.0 && particles[i]->v(1)<0) {
			particles[i]->x(1) = -2.0;
			particles[i]->v(1) = 0;
		}
	}

	for (int imus = 0; imus < muscles.size(); imus++) {
		muscles[imus]->step(particles, model);
	}
	updatePosNor();
}

void FemMuller::init() {
	glGenBuffers(1, &posBufID);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size() * sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);

	glGenBuffers(1, &eleBufID);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, eleBuf.size() * sizeof(unsigned int), &eleBuf[0], GL_STATIC_DRAW);

	glGenBuffers(1, &norBufID);
	glBindBuffer(GL_ARRAY_BUFFER, norBufID);
	glBufferData(GL_ARRAY_BUFFER, norBuf.size() * sizeof(float), &norBuf[0], GL_DYNAMIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	assert(glGetError() == GL_NO_ERROR);
}

void FemMuller::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> p) const
{
	for (int imus = 0; imus < muscles.size(); imus++) {
		muscles[imus]->draw(MV, p);
	}

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

	glDisableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	MV->popMatrix();
}

void FemMuller::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> p, const shared_ptr<Program> p1, const shared_ptr<MatrixStack> P) const
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

	glDrawElements(GL_TRIANGLES, 3 * nTriFaces, GL_UNSIGNED_INT, 0);

	//(const void *)(0 * sizeof(unsigned int))
	glDisableVertexAttribArray(h_pos);
	glDisableVertexAttribArray(h_nor);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	MV->popMatrix();

	p->unbind();

	p1->bind();
	glUniformMatrix4fv(p1->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniformMatrix4fv(p1->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	for (int imus = 0; imus < muscles.size(); imus++) {
		muscles[imus]->draw(MV, p1);
	}
	p1->unbind();
}

void FemMuller::tare() {
	for (int k = 0; k < particles.size(); k++) {
		particles[k]->tare();
	}
}

void FemMuller::reset() {
	for (int k = 0; k < particles.size(); k++) {
		particles[k]->reset();
	}
	updatePosNor();
}

void FemMuller::updatePosNor()
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

MatrixXd FemMuller::hooke(double young, double poisson) {
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

bool FemMuller::rayTriangleIntersects(Vector3d v1, Vector3d v2, Vector3d v3, Vector3d dir, Vector3d pos, double &t, double &u, double &v) {

	Vector3d e1 = v2 - v1;
	Vector3d e2 = v3 - v1;

	// Calculate planes normal vector
	//cross product
	Vector3d pvec = dir.cross(e2);

	//dot product
	double det = e1.dot(pvec);

	// Ray is parallel to plane
	if (det <1e-8 && det > -1e-8) {
		return false;
	}

	double inv_det = 1 / det;

	// Distance from v1 to ray pos
	Vector3d tvec = pos - v1;
	u = (tvec.dot(pvec))*inv_det;
	if (u < 0 || u > 1) {
		return false;
	}

	Vector3d qvec = tvec.cross(e1);
	v = dir.dot(qvec) * inv_det;
	if (v<0 || u + v>1) {
		return false;
	}

	t = e2.dot(qvec) * inv_det;
	if (t > 1e-8) { return true; }
	return false;
}

FemMuller::~FemMuller()
{
}