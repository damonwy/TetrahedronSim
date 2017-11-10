#include "Tet.h"
#include "Particle.h"

using namespace std;
using namespace Eigen;

Tet::Tet(std::shared_ptr<Particle> p0, std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2, std::shared_ptr<Particle> p3){
	this->p0 = p0;
	this->p1 = p1;
	this->p2 = p2;
	this->p3 = p3;
	R_rest.setZero();
	rotatedRestPos.resize(4,3);
	rotatedRestPos.setZero();

}

Tet::~Tet()
{

}
