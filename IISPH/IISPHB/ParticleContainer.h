#pragma once
#include <vector>
#include "Particle.h"

class ParticleContainer
{
public:
	float test;
	int gridSize_;
	int fluidNumber_;
	bool rigidbody;
	std::vector<Particle> fluidParticles;
	std::vector<int> movingParticles;
	std::vector<int> rigidParticles;
	std::vector<int> outerLayerParticles;



	ParticleContainer();


	int getNumberOfParticles() const { return gridSize_; }
	void setParticles();
	void setParticles1();
	void setParticles2();
	void setParticles3();
	void setParticles4();
	void setParticles5();
	void setParticles6();
	void setParticles7();
	void setParticles7_2();
	void setParticles8();
	void setParticles9();
	void setParticles10();
	void setParticles11();
	void setParticles13();
	void setParticles14();
	void setParticles15();
	void setParticles16();
	void setParticles17();
	void setParticles18();
	void setParticles19();
	void setParticlesBowl();
	



	

	void moveBoundaryParticles();

};
