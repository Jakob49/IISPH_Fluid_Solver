
#include <vector>
#include <chrono>
#include "Particle.h"
#include "ParticleContainer.h"
#include "UniformGrid.h"
#include "Globals.h"
#include "Helper.h"


class ISPH {

public:

	void Update(ParticleContainer& pc, UniformGrid& uniform_grid);

};