#include "Muscle.h"
#include "Particle.h"
#include "Program.h"
#include "GLSL.h"
#include "MatrixStack.h"
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>

using namespace std;
using namespace Eigen;

Segment::Segment() {

}

Segment::Segment(shared_ptr<Particle> p0, shared_ptr<Particle> p1) {
	this->p0 = p0;
	this->p1 = p1;
	this->L0 = (p0->x - p1->x).norm();
	this->Ld = this->L0;
	this->L = this->L0;
	this->k = 1e3;
}

Muscle::Muscle() {

}

Muscle::Muscle(shared_ptr<Particle> p0, shared_ptr<Particle> p1) :
	E(1.0)
{
	assert(p0);
	assert(p1);
	assert(p0 != p1);
	this->p0 = p0;
	this->p1 = p1;
	this->time = 0;
	this->numElements = 0;
	this->theta = 0;
	this->isActive = false;
	this->direction = (p1->x - p0->x).normalized();
}

void Muscle::updatePos() {
	for (int i = 0; i < numElements; i++) {
		Vector3d p0 = segments[i]->p0->x;
		Vector3d p1 = segments[i]->p1->x;
		for (int t = 0; t < 3; t++) {
			posBuf[6 * i + 0 + t] = p0(t);
			posBuf[6 * i + 3 + t] = p1(t);
		}
	}
}

void Muscle::init() {
	posBuf.clear();
	eleBuf.clear();
	posBuf.resize(6 * numElements);
	eleBuf.resize(2 * numElements);
	updatePos();

	for (int i = 0; i < numElements; i++) {
		eleBuf[2 * i + 0] = 2 * i + 0;
		eleBuf[2 * i + 1] = 2 * i + 1;
	}

	glGenBuffers(1, &posBufID);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size() * sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);

	glGenBuffers(1, &eleBufID);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, eleBuf.size() * sizeof(unsigned int), &eleBuf[0], GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	assert(glGetError() == GL_NO_ERROR);
}

void Muscle::insertSegment(std::shared_ptr<Segment> segment) {
	segments.push_back(segment);
	elementIDs.push_back(segment->eleID);
	numElements += 1;
}

void Muscle::step(std::vector< std::shared_ptr<Particle> > particles, int model) {
	time += 1;
	if (time % 1 == 0 && model == 0) {
		for (int iseg = 0; iseg < numElements; iseg++) {
			this->segments[iseg]->Ld = (1.5 + 1 * sin(theta * 3.141592 / 180.0))* this->segments[iseg]->L0;
		}
		theta += 1;
	}

	if (time % 10 == 0 && model == 1) {
		for (int iseg = 0; iseg < numElements; iseg++) {
			this->segments[iseg]->Ld = (2.5 + 2 * sin(theta * 3.141592 / 180.0))* this->segments[iseg]->L0;
		}
		theta += 1;
	}

	// update muscle position
	for (int i = 0; i < numElements; i++) {
		double u0 = segments[i]->p0->U;
		double v0 = segments[i]->p0->V;
		Vector3i triIndex0 = segments[i]->p0->triIndex;
		Vector3d p0 = (1 - u0 - v0) * particles[triIndex0(0)]->x + u0 * particles[triIndex0(1)]->x + v0 * particles[triIndex0(2)]->x;
		segments[i]->p0->x = p0;

		double u1 = segments[i]->p1->U;
		double v1 = segments[i]->p1->V;
		Vector3i triIndex1 = segments[i]->p1->triIndex;
		Vector3d p1 = (1 - u1 - v1) * particles[triIndex1(0)]->x + u1 * particles[triIndex1(1)]->x + v1 * particles[triIndex1(2)]->x;
		segments[i]->p1->x = p1;

		// update muscle segment length
		segments[i]->L = (segments[i]->p0->x - segments[i]->p1->x).norm();
	}
	updatePos();
}

void Muscle::draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const {
	MV->pushMatrix();
	glUniformMatrix4fv(p->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	if (this->isActive) {
		glColor3f(0.0, 1.0, 0.0);   // green 
		glLineWidth(3.5);
	}
	else {
		glColor3f(0.0, 0.0, 1.0);   // BLUE
		glLineWidth(2.5);
	}
	glBegin(GL_LINES);
	for (int iseg = 0; iseg < numElements; iseg++) {
		glVertex3f(posBuf[6 * iseg], posBuf[6 * iseg + 1], posBuf[6 * iseg + 2]);
		glVertex3f(posBuf[6 * iseg + 3], posBuf[6 * iseg + 4], posBuf[6 * iseg + 5]);
	}
	glEnd();
	MV->popMatrix();
}

Muscle::~Muscle()
{
}