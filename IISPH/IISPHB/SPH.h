#include <vector>
#include "Particle.h"
#include "ParticleContainer.h"
#include "UniformGrid.h"


class SPH {

public:

	void Update(ParticleContainer& pc, UniformGrid& uniform_grid);

};