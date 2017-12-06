#pragma once

#include <vector>
#include <memory>
#include <tetgen.h>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

class Particle;
class Spring;
class MatrixStack;
class Program;
class Muscle;
class Segment;

class FemMuller
{
public:
	FemMuller(
		double density,
		const Eigen::Vector2d &damping);

	virtual ~FemMuller();

	void step(double h, const Eigen::Vector3d &grav, const bool *keyToggles);
	void init();
	void tare();
	void reset();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
	Eigen::MatrixXd hooke(double young, double poisson);
	void updatePosNor();
	bool rayTriangleIntersects(Eigen::Vector3d v1, Eigen::Vector3d v2, Eigen::Vector3d v3, Eigen::Vector3d dir, Eigen::Vector3d pos, double &t, double &u, double &v);
	void createMuscle(Eigen::Vector3d p0, Eigen::Vector3d p1);
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p, const std::shared_ptr<Program> p1, const std::shared_ptr<MatrixStack> P) const;
private:
	double young;
	double poisson;
	int model;
	int steps;
	int n;
	int nVerts;
	int nFacets;
	int nTriFaces;
	int nTets;
	tetgenio in, out, in_2;
	bool isTop;
	bool isBottom;
	bool isLeft;
	bool isRight;

	int a, b, c, d;
	Eigen::Matrix3d dp, U, R, FF, stress, P;
	Eigen::MatrixXd RKR;
	Eigen::MatrixXd Ke;
	Eigen::MatrixXd Re;

	Eigen::VectorXd xx;
	Eigen::VectorXd RKX;
	Eigen::Vector3d xa, pb, pc, pd, a0, a1, r0, r1, r2;

	Eigen::Matrix3d I;
	Eigen::Matrix3d X_inv;
	Eigen::MatrixXd X_invs;
	std::vector< std::shared_ptr<Particle> > particles;
	std::vector < std::shared_ptr<Muscle> > muscles;

	Eigen::Vector3d direction;
	Eigen::VectorXd mass;
	Eigen::VectorXd volume;
	Eigen::VectorXd v;
	Eigen::VectorXd F0;
	Eigen::VectorXd F;
	Eigen::VectorXd X;
	Eigen::Vector2d damping;

	Eigen::MatrixXd D; // symmetric 6x6 material stiffness matrix
	Eigen::MatrixXd M;
	Eigen::MatrixXd K;
	Eigen::MatrixXd Kes;

	std::vector<unsigned int> eleBuf;
	std::vector<float> posBuf;
	std::vector<float> norBuf;
	std::vector<float> texBuf;

	unsigned eleBufID;
	unsigned posBufID;
	unsigned norBufID;
	unsigned texBufID;
};
