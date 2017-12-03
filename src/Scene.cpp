#include <iostream>

#include "Scene.h"
#include "Particle.h"
#include "Cloth.h"
#include "Shape.h"
#include "Program.h"
#include "FemTet.h"
#include "FemNesme.h"
#include "FemSimit.h"
#include "FemMuller.h"

using namespace std;
using namespace Eigen;

// 1: LINEAR 2: STVK 3: NEOHOOKEAN 4: COROTATED 6: Simit 7: Nesme

// 10: FemTet 11: FemNesme 12: FemSimit 13: FemMuller
// 20: RigidBody
#define MODEL 13

Scene::Scene() :
	t(0.0),
	h(1e-2),
	grav(0.0, 0.0, 0.0)
{
}

Scene::~Scene()
{
}

void Scene::load(const string &RESOURCE_DIR)
{

	// Units: meters, kilograms, seconds
	grav << 0.0, -10, 0.0;
	double mass = 166.667;
	double density = 0.1;
	double height = 12.0;

	Vector2d damping(0.8, 0.8);

	if (MODEL == 10) {
		h = 1e-2;
		femtet = make_shared<FemTet>(density, damping);
	}
	if (MODEL == 11) {
		h = 1e-3;
		femNesme = make_shared<FemNesme>(density, damping);
	}
	if (MODEL == 12) {
		h = 1e-3;
		femSimit = make_shared<FemSimit>(density, damping);
	}
	if (MODEL == 13) {
		h = 1e-2;
		femMuller = make_shared<FemMuller>(density, damping);
	}
}

void Scene::init()
{
	if (MODEL == 10) {
		femtet->init();
	}
	if (MODEL == 11) {
		femNesme->init();
	}

	if (MODEL == 12) {
		femSimit->init();
	}
	if (MODEL == 13) {
		femMuller->init();
	}
}

void Scene::tare()
{
}

void Scene::reset()
{
	t = 0.0;
	//femtet->reset();
	//femNesme->reset();

}

void Scene::step()
{
	t += h;

	if (MODEL == 10) {
		femtet->step(h, grav);
	}
	if (MODEL == 11) {
		femNesme->step(h, grav);
	}
	if (MODEL == 12) {
		femSimit->step(h, grav);
	}
	if (MODEL == 13) {
		femMuller->step(h, grav);
	}
}

void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog) const
{
	glUniform3fv(prog->getUniform("kdFront"), 1, Vector3f(1.0, 1.0, 1.0).data());

	if (MODEL == 10) {
		femtet->draw(MV, prog);
	}
	if (MODEL == 11) {
		femNesme->draw(MV, prog);
	}
	if (MODEL == 12) {
		femSimit->draw(MV, prog);
	}
	if (MODEL == 13) {
		femMuller->draw(MV, prog);
	}
}
