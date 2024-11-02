#include "Particle.h"
#include <iostream>

#define _USE_MATH_DEFINES
#include<math.h>
#include <cmath>

#include "Globals.h"


Particle::Particle(sf::Vector2f position, bool bound)
{
	mass_ = mass;
	position_ = position;
	velocity_ = sf::Vector2f(0, 0);
	acceleration_ = sf::Vector2f(0, 0);
	npacceleration_ = sf::Vector2f(0, 0);
	pacceleration_ = sf::Vector2f(0, 0); 
	predictedVelocity_ = sf::Vector2f(0, 0);
	particlePressureAcceleration = sf::Vector2f(0, 0);
	boundaryPressureAcceleration = sf::Vector2f(0, 0);
	transformedPosition = sf::Vector2f(0, 0);


	kernel_ = 0;
	color_ = sf::Color::Blue;
	radius_ = 0.505f; // 0.505f
	pressure_ = 0.0f;
	density_ = p0;
	boundary_ = bound;
	fluidNeighbors_;
	densityError_ = 0;

	aff_ = 0;
	ap_ = 0;
	sf_ = 0;


	interpolateColors = { sf::Color::Green, sf::Color::Red, sf::Color::Green };

	if (bound)
	{
		color_ = sf::Color(75, 75, 75);
	}
	shape_.setRadius(radius_);
	shape_.setPosition(position_);
	shape_.setFillColor(color_);


}

void Particle::Update()
{
	shape_.setPosition(position_);
	shape_.setFillColor(color_);
}

void Particle::setColor(sf::Color color)
{
	color_ = color;
	shape_.setFillColor(color_);
}


sf::CircleShape Particle::getShape() const {
	return shape_;
}

std::vector<size_t> Particle::getFluidNeighbors()
{
	return fluidNeighbors_;
}


void Particle::setDensity(float temp)
{
	density_ = temp;

	densityError_ = (density_ - p0) / p0;

}

void Particle::setPressure(float pressure)
{
	pressure_ = pressure;
}
void Particle::setPredictedVelocity() {
	predictedVelocity_ = velocity_ + t * npacceleration_;

}


void Particle::setVelocity()
{
	velocity_ = velocity_ + t * acceleration_;


	if (std::isinf(velocity_.x) || std::isinf(velocity_.y)) {
		std::cout << "here aff  " << std::endl;
		std::cout << "at particle " << "   aff   " << velocity_.x << "    sf:   " << velocity_.y << "    pressure   " << std::endl;
	}
}
void Particle::setBoundaryVelocity(sf::Vector2f vel) {
	velocity_ = vel;
}
void Particle::setBoundaryPosition(sf::Vector2f pos) {
	//std::cout << "pos.x" << pos.x << "pos.y" << pos.y << std::endl;
	position_ = pos;
	UpdateShape();
}


void Particle::setAcceleration()
{
	acceleration_ = npacceleration_ + pacceleration_;

}

void Particle::setNonPAcceleration(sf::Vector2f acc)
{
	npacceleration_ = acc;
}

void Particle::setPAcceleration(sf::Vector2f acc)
{
	pacceleration_ = acc;
}


void Particle::setPosition()
{
	position_ = position_ + t * velocity_;
	UpdateShape();
}





void Particle::UpdateShape()
{
	shape_.setPosition(position_);
}

void Particle::setNeighborIndices(const std::vector<size_t>& neighborIndices)
{
	fluidNeighbors_.clear();
	fluidNeighbors_ = neighborIndices;
}
