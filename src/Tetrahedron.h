#pragma once


#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

class Particle;
class Spring;
class MatrixStack;
class Program;

class Tetrahedron
{
public:
	Tetrahedron(
		const Eigen::Vector3d &x0,
		const Eigen::Vector3d &x1,
		const Eigen::Vector3d &x2,
		const Eigen::Vector3d &x3,
		double mass,
		double stiffness,
		const Eigen::Vector2d &damping);

	virtual ~Tetrahedron();

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
	Eigen::Matrix3d X_inv;
	std::vector< std::shared_ptr<Particle> > particles;
	Eigen::VectorXd v;
	Eigen::VectorXd f;
	Eigen::VectorXd p;
	Eigen::Vector2d damping;

	Eigen::MatrixXd M;
	Eigen::MatrixXd K;

	std::vector<unsigned int> eleBuf;
	std::vector<float> posBuf;
	std::vector<float> norBuf;
	std::vector<float> texBuf;

	unsigned eleBufID;
	unsigned posBufID;
	unsigned norBufID;
	unsigned texBufID;
};

