#include <iostream>

#include "Scene.h"
#include "Particle.h"
//#include "Cloth.h"
//#include "Tetrahedron.h"
#include "Bar.h"
#include "Shape.h"
#include "Program.h"

using namespace std;
using namespace Eigen;

Scene::Scene() :
	t(0.0),
	h(1e-3),
	grav(0.0, 0.0, 0.0)
{
}

Scene::~Scene()
{
}

void Scene::load(const string &RESOURCE_DIR)
{
	// Units: meters, kilograms, seconds
	h = 2e-3;
	
	grav << 0.0, -10, 0.0;
	
	//int rows = 2;
	//int cols = 2;
	double mass = 166.667;
	double density = 10;
	double height = 12.0;
	//double stiffness = 1e2;
	Vector2d damping(0.0, 1.0);

	// Cloth 

	//Vector3d x00(-0.25, 0.5, 0.0);
	//Vector3d x01(0.25, 0.5, 0.0);
	//Vector3d x10(-0.25, 0.5, -0.5);
	//Vector3d x11(0.25, 0.5, -0.5);
	//Vector3d x0(0.0, 0.95, -0.05);
	//Vector3d x1(-0.1, 0.85, -0.05);
	//Vector3d x2(0.0, 0.9, 0.0);
	//Vector3d x3(0.1, 0.85, -0.05);


	// Tetrahedron

	Vector3d x0(0.0, 5.0, 1.0);
	Vector3d x1(0.0, 5.0, 0.0);
	Vector3d x2(1.0, 5.0, 0.0);
	Vector3d x3(0.0, 6.0, 0.0);


	//cloth = make_shared<Cloth>(rows, cols, x00, x01, x10, x11, mass, stiffness, damping);
	
	//tet = make_shared<Tetrahedron>(x0, x1, x2, x3, mass, stiffness, damping);
	bar = make_shared<Bar>(x0, x1, x2, x3, density, height, damping);
	//sphereShape = make_shared<Shape>();
	//sphereShape->loadMesh(RESOURCE_DIR + "sphere2.obj");
	
	//auto sphere = make_shared<Particle>(sphereShape);
	//spheres.push_back(sphere);
	//sphere->r = 0.1;
	//sphere->x = Vector3d(0.0, 0.2, 0.0);
}

void Scene::init()
{
	//sphereShape->init();
	//cloth->init();
	//tet->init();
	bar->init();
}

void Scene::tare()
{
	/*for(int i = 0; i < (int)spheres.size(); ++i) {
		spheres[i]->tare();
	}*/
	//cloth->tare();
	//tet->tare();
	bar->tare();
}

void Scene::reset()
{
	t = 0.0;
	/*for(int i = 0; i < (int)spheres.size(); ++i) {
		spheres[i]->reset();
	}*/
	//cloth->reset();
	//tet->reset();
	bar->reset();
}

void Scene::step()
{
	t += h;
	
	// Move the sphere
	/*if(!spheres.empty()) {
		auto s = spheres.front();
		Vector3d x0 = s->x;
		double radius = 0.5;
		double a = 2.0*t;
		s->x(2) = radius * sin(a);
		Vector3d dx = s->x - x0;
		s->v = dx/h;
	}*/
	
	// Simulate the cloth
	//cloth->step(h, grav, spheres);
	//tet->step(h, grav);
	bar->step(h, grav);
}

void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog) const
{
	glUniform3fv(prog->getUniform("kdFront"), 1, Vector3f(1.0, 1.0, 1.0).data());
	//for(int i = 0; i < (int)spheres.size(); ++i) {
	//	//spheres[i]->draw(MV, prog);
	//}
	//cloth->draw(MV, prog);
	//tet->draw(MV, prog);
	bar->draw(MV, prog);
}
