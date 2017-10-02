#pragma once
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

class Bar
{
public:
	Bar(
		const Eigen::Vector3d &x0,
		const Eigen::Vector3d &x1,
		const Eigen::Vector3d &x2,
		const Eigen::Vector3d &x3,
		double density,
		double height,
		const Eigen::Vector2d &damping);

	virtual ~Bar();

	void step(double h, const Eigen::Vector3d &grav);
	void init();
	void tare();
	void reset();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
	Eigen::MatrixXd hooke(double young, double poisson);
	void updatePosNor();

private:
	double young;
	double poisson;
	int steps;
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
	Eigen::VectorXd v;
	Eigen::VectorXd F;
	Eigen::VectorXd X;
	Eigen::Vector2d damping;

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

