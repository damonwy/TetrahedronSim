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

Segment::Segment(shared_ptr<Particle> p0, shared_ptr<Particle> p1) {
	this->p0 = p0;
	this->p1 = p1;
	this->L0 = (p0->x - p1->x).norm();
	this->Ld = 2 * this->L0;
	this->L = this->L0;
	this->k = 1e2;
}

Muscle::Muscle(shared_ptr<Particle> p0, shared_ptr<Particle> p1) :
	E(1.0)
{
	assert(p0);
	assert(p1);
	assert(p0 != p1);
	this->p0 = p0;
	this->p1 = p1;
	this->time = 0.0;
	this->numElements = 0;	
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
 
void Muscle::step(std::vector< std::shared_ptr<Particle> > particles) {
	time += 1.0;
	if (time) {
		for (int iseg = 0; iseg < numElements; iseg++) {

		}
		this->Ld = 2 * this->L0;
	}

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
	glUniform3fv(p->getUniform("kdFront"), 1, Vector3f(1.0, 0.0, 0.0).data());
	glUniform3fv(p->getUniform("kdBack"), 1, Vector3f(1.0, 1.0, 0.0).data());
	MV->pushMatrix();

	glUniformMatrix4fv(p->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));

	int h_pos = p->getAttribute("aPos");
	glEnableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size() * sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);
	glVertexAttribPointer(h_pos, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);
	glLineWidth(3.0);
	glDrawElements(GL_LINES, 2 * numElements, GL_UNSIGNED_INT, (const void *)(0 * sizeof(unsigned int)));

	glDisableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	MV->popMatrix();
}

Muscle::~Muscle()
{
}
