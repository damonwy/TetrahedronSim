#pragma once
#ifndef __Tet__
#define __Tet__

#include <memory>
#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

class Tet;

class Tet
{
public:
	Tet(std::shared_ptr<Particle> p0, std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2, std::shared_ptr<Particle> p3);
	virtual ~Tet();

	std::shared_ptr<Particle> p0;
	std::shared_ptr<Particle> p1;
	std::shared_ptr<Particle> p2;
	std::shared_ptr<Particle> p3;

	Eigen::Matrix3d R_rest;
	Eigen::MatrixXd rotatedRestPos;

};

#endif