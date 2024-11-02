#pragma once
#include <vector>
#include <iostream>
#include <SFML/System/Vector2.hpp>
#include "ParticleContainer.h"

class UniformGrid
{


public:
	UniformGrid();
	void GridClear();
	void AddParticle(sf::Vector2f position, size_t index);
	std::vector<size_t> getNeighbors(size_t index, sf::Vector2f position, ParticleContainer& pc);
	std::vector<std::vector<std::vector<size_t>>> cells_;
	float CalulateDistance(sf::Vector2f a, sf::Vector2f b);
	int GiveClosestParticle(sf::Vector2f click, ParticleContainer pc);

};

