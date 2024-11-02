#pragma once
#include <vector>
#include "Particle.h"
#include "ParticleContainer.h"
#include "UniformGrid.h"
#include "Globals.h"

class Helper
{
public:
	float CalculateDistance(sf::Vector2f a, sf::Vector2f b);
	static float ScalarProduct(sf::Vector2f a, sf::Vector2f b);
	static float MagniVector(sf::Vector2f vec);
	static sf::Vector2f gradientKernel(sf::Vector2f positionX, sf::Vector2f positionY);
	float CalculateDensity(size_t i, ParticleContainer& pc);


private:

};
