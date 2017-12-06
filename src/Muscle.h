#pragma once
#ifndef __Muscle__
#define __Muscle__
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <memory>
using namespace Eigen;
class Particle;
class MatrixStack;
class Program;

class Segment
{
public:
	double L0; // initial length
	double L;  // current length
	double Ld; // designed length
	double k;  // stiffness
	int eleID; // index of element

	std::shared_ptr<Particle> p0;
	std::shared_ptr<Particle> p1;
	Segment();
	Segment(std::shared_ptr<Particle> p0, std::shared_ptr<Particle> p1);
};

class Muscle
{
public:
	int numElements; // number of muscle segments
	int type; // LONGITUDINAL
	int muscleID;

	std::shared_ptr<Particle> p0;
	std::shared_ptr<Particle> p1;
	Eigen::Vector3d direction;

	double E;
	double L;
	double L0;
	double Ld;
	int time;
	int theta;
	bool isActive;
	Muscle();
	Muscle(std::shared_ptr<Particle> p0, std::shared_ptr<Particle> p1);
	virtual ~Muscle();

	std::vector<int> elementIDs;
	std::vector<std::shared_ptr<Segment>> segments;

	void insertSegment(std::shared_ptr<Segment> segment);
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
	void init();
	void updatePos();
	void step(std::vector< std::shared_ptr<Particle> > particles, int model);

private:
	unsigned eleBufID;
	unsigned posBufID;
	std::vector<unsigned int> eleBuf;
	std::vector<float> posBuf;
};
#endif