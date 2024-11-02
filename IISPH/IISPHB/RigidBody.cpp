#include "RigidBody.h"


RigidBody::RigidBody(std::vector<int> pindices, std::vector<int> outerLayerParticles) {
	particleIndices = pindices;
	centerOfMass = sf::Vector2f(0.0f, 0.0f);
	bodymass = 0.0f;
	inverseInertiaTensor = 0.0f;
	velocityCenterOfMass = sf::Vector2f(0.0f, 0.0f);
	angularMomentum = 0.0f;
	angularVelocity = 0.0f;
	rotationMatrix[0][0] = 1.0f;
	rotationMatrix[0][1] = 0.0f;
	rotationMatrix[1][0] = 0.0f;
	rotationMatrix[1][1] = 1.0f;
	torque = 0.0f;
	outer = outerLayerParticles;
}


void RigidBody::Initialize(ParticleContainer& pc) {
	// Mass computation
	for (auto& i : particleIndices) {

		pc.fluidParticles[i].mass_ = pc.fluidParticles[i].mass_ * 8.0f;
		bodymass = bodymass + pc.fluidParticles[i].mass_;
	}
	
	for (auto& i : particleIndices) {
		centerOfMass = centerOfMass + pc.fluidParticles[i].getPosition() * pc.fluidParticles[i].mass_;
	}
	centerOfMass = centerOfMass / bodymass;

	relativePosition.clear();
	for (auto& i : particleIndices) {
		relativePosition.push_back(pc.fluidParticles[i].getPosition() - centerOfMass);
		forces.push_back(sf::Vector2f(0.0f, 0.0f));
	}



	float inertiaTensor = 0.0f;
	for (int i = 0; i < particleIndices.size(); i++) {
		inertiaTensor +=  pc.fluidParticles[particleIndices[i]].mass_ * (relativePosition[i].x * relativePosition[i].x + relativePosition[i].y * relativePosition[i].y);
	}
	

	inverseInertiaTensor = (inertiaTensor != 0.0f) ? 1.0f / inertiaTensor : 0.0f;


	angularVelocity = inverseInertiaTensor * angularMomentum;


	sf::Vector2f totalRelativePosition(0.0f, 0.0f);

	// Sum up all the relative positions
	for (auto& relPos : relativePosition) {
		totalRelativePosition += relPos;
		std::cout << "Relative Position x  " << relPos.x << " Relative Position y  "<<  relPos.y << std::endl;
	}

	// Check if the total sum of relative positions is near zero
	float epsilon = 1e-3; // Tolerance for floating point errors
	if (std::abs(totalRelativePosition.x) < epsilon && std::abs(totalRelativePosition.y) < epsilon) {
		std::cout << "Relative positions are symmetric (sum near zero)." << std::endl;
	}
	else {
		std::cout << "Relative positions are not symmetric." << std::endl;
		std::cout << "Total relative position sum: (" << totalRelativePosition.x << ", " << totalRelativePosition.y << ")" << std::endl;
	}

	//std::cout << "number of particles" << particleIndices.size() << "first particle  " << particleIndices[0] << " center of mass x" << centerOfMass.x << " center of mass y" << centerOfMass.y << std::endl;
	//for (auto& i : relativePosition) {
	//	std::cout << " relative x" << i.x << " relative y" << i.y << std::endl;
	//}
}

void RigidBody::Update(ParticleContainer& pc) {
	sf::Vector2f accviscosity;
	sf::Vector2f accpressure;
	sf::Vector2f total;
	sf::Vector2f force;
	double acctorque = 0.0f;
	torque = 0;


	// only for the border particle and the neighbors should only be fluid particles
	int counting = 0;
	// only outer particles
	for (auto& k : outer) {
		int i = particleIndices[k];
		pc.fluidParticles[i].setColor(sf::Color(75, 75, 75));
		std::vector<size_t> neighbors = pc.fluidParticles[i].getFluidNeighbors();
		for (auto& j : neighbors) {
			if (pc.fluidParticles[j].boundary_ == false) {
				float pi = ((Helper::ScalarProduct(pc.fluidParticles[j].getVelocity() - pc.fluidParticles[i].getVelocity(), pc.fluidParticles[j].getPosition() - pc.fluidParticles[i].getPosition()))
					/ (Helper::ScalarProduct(pc.fluidParticles[j].getPosition() - pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition() - pc.fluidParticles[i].getPosition()) + 0.01 * h * h));
				accviscosity = accviscosity + (float)(pc.fluidParticles[j].mass_ / pc.fluidParticles[j].getDensity()) * pi * Helper::gradientKernel(pc.fluidParticles[j].getPosition(), pc.fluidParticles[i].getPosition());

				accpressure = accpressure + (-1.0f * pc.fluidParticles[i].mass_ * pc.fluidParticles[j].mass_ *
					(pc.fluidParticles[j].getPressure() / (p0 * p0)) * Helper::gradientKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition()));

				//accviscosity = accviscosity - pc.fluidParticles[j].getNonPAcceleration();
				//accpressure = accpressure - pc.fluidParticles[j].getPAcceleration();
				//force = force + accpressure + accviscosity;
				force = force +  accpressure + accviscosity;
			}
			accviscosity = sf::Vector2f(0.0f, 0.0f);
			accpressure =  sf::Vector2f(0.0f, 0.0f);
		}
		force = force + sf::Vector2f(0.0f, gravity * pc.fluidParticles[i].mass_);

		acctorque = relativePosition[k].x * force.y - relativePosition[k].y * force.x;

		/*std::cout << "Particle " << i << " Torque Contribution: " << acctorque << std::endl;*/
		torque = torque + acctorque;
		total = total + force;
		forces[k] = force;

		force = sf::Vector2f(0.0f, 0.0f);

		counting++;
	}
	//std::cout << "total force x " << total.x << "   total force y " << total.y << std::endl;
	//std::cout << "angular velocity" << angularVelocity << std::endl;
	//std::cout << "torque " << torque << std::endl;
	// std::cout << "center of Mass x " << centerOfMass.x <<" center of mass.y  " << centerOfMass.y << std::endl;



	// update center of mass
	centerOfMass += velocityCenterOfMass * t;


	//  Update velocity of the center of mass
	velocityCenterOfMass += (total / bodymass) * t;

	// Update the rotation matrix A

	float a11 = rotationMatrix[0][0];
	float a12 = rotationMatrix[0][1];
	float a21 = rotationMatrix[1][0];
	float a22 = rotationMatrix[1][1];

	rotationMatrix[0][0] = a11 - t * angularVelocity * a21;
	rotationMatrix[0][1] = a12 - t * angularVelocity * a22;
	rotationMatrix[1][0] = a21 + t * angularVelocity * a11;
	rotationMatrix[1][1] = a22 + t * angularVelocity * a12;


	// update the angular momentum L 
	angularMomentum = angularMomentum + t * torque;

	// update the inverse of the inertia tensor apparently no need to update
	// inverseInertiaTensor = 

	angularVelocity = inverseInertiaTensor * angularMomentum;
	int count = 0;
	for (auto& i : particleIndices) {
		// std::cout << "relative Position " << relativePosition[count].x << "   total force y " << relativePosition[count].y << std::endl;
		// update relative positions
		float old_x = relativePosition[count].x;
		float old_y = relativePosition[count].y;
		relativePosition[count].x = old_x * rotationMatrix[0][0] + old_y * rotationMatrix[0][1];
		relativePosition[count].y = old_x * rotationMatrix[1][0] + old_y * rotationMatrix[1][1];
		// update boundary positions
		pc.fluidParticles[i].setBoundaryPosition(centerOfMass + relativePosition[count]);

		// update velocities

		count++;
	}

	// Gram Schmidt
	//a11 = rotationMatrix[0][0];
	//a12 = rotationMatrix[0][1];
	//a21 = rotationMatrix[1][0];
	//a22 = rotationMatrix[1][1];


	//float length1 = sqrt(a11 * a11 + a21 * a21);
	//float b11 = a11 / length1;
	//float b21 = a21 / length1;

	//
	//float dotProduct = b11 * a12 + b21 * a22;
	//float a12_prime = a12 - dotProduct * b11;
	//float a22_prime = a22 - dotProduct * b21;

	//
	//float length2 = sqrt(a12_prime * a12_prime + a22_prime * a22_prime);
	//float b12 = a12_prime / length2;
	//float b22 = a22_prime / length2;

	//
	//rotationMatrix[0][0] = b11;
	//rotationMatrix[0][1] = b12;
	//rotationMatrix[1][0] = b21;
	//rotationMatrix[1][1] = b22;


}