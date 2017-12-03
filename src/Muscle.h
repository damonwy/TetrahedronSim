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

	int faceinID; // the index of face in one tet 
	int faceoutID; 

	std::shared_ptr<Particle> p0;
	std::shared_ptr<Particle> p1;
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

	double E;
	double L;
	double L0;

	double time;

	Muscle(std::shared_ptr<Particle> p0, std::shared_ptr<Particle> p1);
	virtual ~Muscle();

	std::vector<int> elementIDs;
	std::vector<std::shared_ptr<Segment>> segments;

	void insertSegment(std::shared_ptr<Segment> segment);
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;
	void init();
	void updatePos();
	void step(std::vector< std::shared_ptr<Particle> > particles);

private:
	unsigned eleBufID;
	unsigned posBufID;
	std::vector<unsigned int> eleBuf;
	std::vector<float> posBuf;
};
#endif