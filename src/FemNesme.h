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

class FemNesme
{
public:
	FemNesme(
		double density,
		const Eigen::Vector2d &damping);

	virtual ~FemNesme();

	void step(double h, const Eigen::Vector3d &grav);
	void init();
	void tare();
	void reset();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
	Eigen::MatrixXd computeD(double lambda, double mu);
	void updatePosNor();

private:
	double young;
	double poisson;
	double mu;
	double lambda;
	int n;
	int nVerts;
	int nFacets;
	int nTriFaces;
	int nTets;
	tetgenio in, out, in_2;

	Eigen::Matrix3d X_inv;
	Eigen::MatrixXd X_invs;
	std::vector< std::shared_ptr<Particle> > particles;
	Eigen::VectorXd mass;
	Eigen::VectorXd volume;
	Eigen::VectorXd v;
	Eigen::VectorXd F;
	Eigen::MatrixXd K;
	Eigen::Vector2d damping;

	Eigen::MatrixXd D; // symmetric 6x6 material stiffness matrix
	Eigen::MatrixXd M;

	std::vector<unsigned int> eleBuf;
	std::vector<float> posBuf;
	std::vector<float> norBuf;
	std::vector<float> texBuf;

	unsigned eleBufID;
	unsigned posBufID;
	unsigned norBufID;
	unsigned texBufID;
};

