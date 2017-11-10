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

	void computeQR(Eigen::Matrix3d &rotation, const Eigen::VectorXd &p, const int &a, const int &b, const int &c);
	void computeB(Eigen::MatrixXd &Bmatrix, Eigen::VectorXd volume, int itet, Eigen::Vector3d rotated_a, Eigen::Vector3d rotated_b, Eigen::Vector3d rotated_c, Eigen::Vector3d rotated_d);
	void computeK(Eigen::MatrixXd &Kmatrix, Eigen::MatrixXd &KRtmatrix, Eigen::MatrixXd Bmatrix, Eigen::Matrix3d Rot);
	void step(double h, const Eigen::Vector3d &grav);
	void init();
	void tare();
	void reset();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
	void computeD(Eigen::MatrixXd &Dmatrix, double lambda, double mu);
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
	
	Eigen::VectorXd mass;
	Eigen::VectorXd volume;
	Eigen::MatrixXd rotatedRestPosMatrix;

	Eigen::VectorXd X0; // rest position
	Eigen::VectorXd x; // current position
	Eigen::VectorXd v;
	Eigen::VectorXd F;
	Eigen::MatrixXd K;
	
	Eigen::Vector2d damping;
	Eigen::MatrixXd D; // symmetric 6x6 material stiffness matrix
	Eigen::MatrixXd M;

	std::vector< std::shared_ptr<Particle> > particles;

	std::vector<unsigned int> eleBuf;
	std::vector<float> posBuf;
	std::vector<float> norBuf;
	std::vector<float> texBuf;

	unsigned eleBufID;
	unsigned posBufID;
	unsigned norBufID;
	unsigned texBufID;
};

