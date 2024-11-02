#pragma once

#include <vector>
#include "Particle.h"
#include "ParticleContainer.h"
#include "UniformGrid.h"
#include "Globals.h"
#include "Helper.h"

#include <SFML/Graphics.hpp>
#include <iostream>


class RigidBody {


public:
	float bodymass;                      
	sf::Vector2f centerOfMass;
	float inverseInertiaTensor;
	float torque;

	sf::Vector2f velocityCenterOfMass;
	
	std::vector<sf::Vector2f> relativePosition;
	std::vector<sf::Vector2f> forces;
	std::vector<sf::Vector2f> initialDistance;
	std::vector<int> RigidboundaryParticles;
	std::vector<int> outer;



	float angularVelocity;
	float angularMomentum;
	float angle; 
	float rotationMatrix[2][2];


	RigidBody(std::vector<int> pinidices, std::vector<int> outerLayerParticles);
	void Initialize(ParticleContainer& pc);
	void Update(ParticleContainer& pc);
	std::vector<int> particleIndices;
};