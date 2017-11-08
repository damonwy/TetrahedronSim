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

class FemSimit
{
public:
	FemSimit(double density, const Eigen::Vector2d &damping);
	virtual ~FemSimit();

	void step(double h, const Eigen::Vector3d &grav);
	void init();
	void tare();
	void reset();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
	Eigen::Matrix3d computeStress(Eigen::Matrix3d dp, int material_type, double mu, double lambda);
	Eigen::Matrix3d computeStressDerivative(Eigen::Matrix3d F, Eigen::Matrix3d dF, int material_type, double mu, double lambda);
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
	std::vector<std::shared_ptr<Particle>> particles;
	Eigen::Vector2d damping;

	Eigen::VectorXd mass;

	Eigen::Matrix3d Dm_1; // Dm inverse
	Eigen::MatrixXd Bm;   // Bm[e]
	Eigen::VectorXd W;    // W is the undeformed volume of Te
	Eigen::MatrixXd P; // the 1st Piola-Kirchhoff stress tensor, 3*nTets x 3

	Eigen::MatrixXd M;
	Eigen::VectorXd v;
	Eigen::VectorXd f;
	Eigen::VectorXd X;
	Eigen::MatrixXd K;
	Eigen::Matrix3d I;

	std::vector<unsigned int> eleBuf;
	std::vector<float> posBuf;
	std::vector<float> norBuf;
	std::vector<float> texBuf;

	unsigned eleBufID;
	unsigned posBufID;
	unsigned norBufID;
	unsigned texBufID;
};
