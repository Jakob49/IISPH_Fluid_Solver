#include "SPH.h"
#include <SFML/Graphics.hpp>
#include <vector>
#include <math.h>
#include <iostream>
#define _USE_MATH_DEFINES
#include<math.h>
#include <cmath>
#include <thread>
#include <fstream>
#include <chrono>
#include <omp.h>

#include "Particle.h"
#include "Globals.h"
#include "ParticleContainer.h"
#include "UniformGrid.h"



float CalculateDistance2(sf::Vector2f a, sf::Vector2f b)
{
    return std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2));
}

float ScalarProduct2(sf::Vector2f a, sf::Vector2f b)
{
    return a.x * b.x + a.y * b.y;
}

float MagnitudeVector2(sf::Vector2f vec) {
    return std::sqrt(vec.x * vec.x + vec.y * vec.y);
}


std::vector<size_t> FindNeighborIndices2(Particle& particle, std::vector<Particle>& particles)
{
    std::vector<size_t> neighborIndices;

    const float epsilon = 0.01f;

    for (size_t i = 0; i < particles.size(); ++i) {
        float distance = CalculateDistance2(particle.getPosition(), particles[i].getPosition());
        if (distance <= 2 * h + epsilon) {
            neighborIndices.push_back(i);
        }
    }
    return neighborIndices;
}


sf::Vector2f gradientKernel2(sf::Vector2f positionX, sf::Vector2f positionY) {

    float r = std::sqrt(std::pow(positionX.x - positionY.x, 2)
        + std::pow(positionX.y - positionY.y, 2)); // distance of the particles
    float d = r / h;

    if (d == 0)
    {
        return sf::Vector2f(0.0f, 0.0f);
    }

    float t1 = std::max(1 - d, 0.0f);
    float t2 = std::max(2 - d, 0.0f);


    sf::Vector2f w;
    w.x = alpha * ((positionX.x - positionY.x) / (d * h)) * (-3 * t2 * t2 + 12 * t1 * t1);
    w.y = alpha * ((positionX.y - positionY.y) / (d * h)) * (-3 * t2 * t2 + 12 * t1 * t1);
    return w;
}




float CalculateDensity2(size_t i, ParticleContainer& pc)
{

    float kernel = 0.0f;
    float temp = 0;
    std::vector<size_t> neighbors = pc.fluidParticles[i].getFluidNeighbors();
    for (auto& j : neighbors)
    {
        float r = std::sqrt(std::pow(pc.fluidParticles[i].getPosition().x - pc.fluidParticles[j].getPosition().x, 2)
            + std::pow(pc.fluidParticles[i].getPosition().y - pc.fluidParticles[j].getPosition().y, 2)); // distance of the particles
        float d = r / h;
        float t1 = std::max(1 - d, 0.0f);
        float t2 = std::max(2 - d, 0.0f);
        float w = alpha * (t2 * t2 * t2 - 4 * t1 * t1 * t1);

        kernel = kernel + w;
        temp = temp + pc.fluidParticles[j].mass_ * w;

    }

    return temp;
}

sf::Vector2f CalculateAcceleration2(size_t i, ParticleContainer& pc)
{
    if (pc.fluidParticles[i].boundary_)
    {
        return sf::Vector2f(0, 0);
    }
    std::vector<size_t> neighbors = pc.fluidParticles[i].getFluidNeighbors();

    sf::Vector2f pressureAcceleration;
    sf::Vector2f viscosityAcceleration;
    sf::Vector2f kernel;

    sf::Vector2f acceleration;

    for (auto& j : neighbors)
    {

        float r = std::sqrt(std::pow(pc.fluidParticles[i].getPosition().x - pc.fluidParticles[j].getPosition().x, 2)
            + std::pow(pc.fluidParticles[i].getPosition().y - pc.fluidParticles[j].getPosition().y, 2)); // distance of the particles
        float d = r / h;

        sf::Vector2f w;
        if (d == 0)
        {
            w = sf::Vector2f(0.0f, 0.0f);
        }
        else
        {
            float t1 = std::max(1 - d, 0.0f);
            float t2 = std::max(2 - d, 0.0f);


            w.x = alpha * ((pc.fluidParticles[i].getPosition().x - pc.fluidParticles[j].getPosition().x) / (d * h)) * (-3 * t2 * t2 + 12 * t1 * t1);
            w.y = alpha * ((pc.fluidParticles[i].getPosition().y - pc.fluidParticles[j].getPosition().y) / (d * h)) * (-3 * t2 * t2 + 12 * t1 * t1);
        }


        kernel = kernel + w;


        if (pc.fluidParticles[j].boundary_)
        {
            pressureAcceleration = pressureAcceleration + pc.fluidParticles[j].mass_ *
                (pc.fluidParticles[i].getPressure() / pc.fluidParticles[i].getDensity() * pc.fluidParticles[i].getDensity() +
                    pc.fluidParticles[i].getPressure() / pc.fluidParticles[i].getDensity() * pc.fluidParticles[i].getDensity()) * gradientKernel2(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition());


        }
        else
        {
            pressureAcceleration = pressureAcceleration + pc.fluidParticles[j].mass_ *
                (pc.fluidParticles[i].getPressure() / pc.fluidParticles[i].getDensity() * pc.fluidParticles[i].getDensity() +
                    pc.fluidParticles[j].getPressure() / pc.fluidParticles[j].getDensity() * pc.fluidParticles[j].getDensity()) * gradientKernel2(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition());

        }


        viscosityAcceleration = viscosityAcceleration + (float)(pc.fluidParticles[j].mass_ / pc.fluidParticles[j].getDensity() *
            ((ScalarProduct2(pc.fluidParticles[i].getVelocity() - pc.fluidParticles[j].getVelocity(), pc.fluidParticles[i].getPosition() - pc.fluidParticles[j].getPosition()))
            / (ScalarProduct2(pc.fluidParticles[i].getPosition() - pc.fluidParticles[j].getPosition(), pc.fluidParticles[i].getPosition() - pc.fluidParticles[j].getPosition()) + 0.01 * h * h))) * gradientKernel2(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition());

    }


    viscosityAcceleration = 2 * v * viscosityAcceleration;


    acceleration = viscosityAcceleration + sf::Vector2f(0, gravity * pc.fluidParticles[i].mass_);



    return acceleration;
}

sf::Vector2f CalculatePressureAcceleration2(size_t i, ParticleContainer& pc) {

    if (pc.fluidParticles[i].boundary_)
    {
        return sf::Vector2f(0, 0);
    }
    std::vector<size_t> neighbors = pc.fluidParticles[i].getFluidNeighbors();


    sf::Vector2f pressureAcceleration;
    for (auto& j : neighbors)
    {

        if (pc.fluidParticles[j].boundary_)
        {
            pressureAcceleration = pressureAcceleration + pc.fluidParticles[j].mass_ *
                (pc.fluidParticles[i].getPressure() / pc.fluidParticles[i].getDensity() * pc.fluidParticles[i].getDensity() +
                    pc.fluidParticles[i].getPressure() / pc.fluidParticles[i].getDensity() * pc.fluidParticles[i].getDensity()) * gradientKernel2(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition());
        }
        else
        {
            pressureAcceleration = pressureAcceleration + pc.fluidParticles[j].mass_ *
                (pc.fluidParticles[i].getPressure() / pc.fluidParticles[i].getDensity() * pc.fluidParticles[i].getDensity() +
                    pc.fluidParticles[j].getPressure() / pc.fluidParticles[j].getDensity() * pc.fluidParticles[j].getDensity()) * gradientKernel2(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition());

        }
    }
    return  -1.0f * pressureAcceleration;
}


void SPH::Update(ParticleContainer& pc, UniformGrid& uniform_grid) {

    auto start_neighbor = std::chrono::high_resolution_clock::now();
//#pragma omp parallel for 
    for (int i = 0; i < pc.fluidParticles.size(); i++)
    {
        if (pc.fluidParticles[i].boundary_ == false)
        {
            // pc.fluidParticles[i].setNeighborIndices(FindNeighborIndices(pc.fluidParticles[i], pc.fluidParticles));
            pc.fluidParticles[i].setNeighborIndices(uniform_grid.getNeighbors(i, pc.fluidParticles[i].getPosition(), pc));

        }
    }
    auto end_neighbor = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_neighbor = end_neighbor - start_neighbor;
    
    // set pressure and density loop

    auto start_pressuredensity = std::chrono::high_resolution_clock::now();
//#pragma omp parallel for 
    for (int i = 0; i < pc.fluidParticles.size(); i++)
    {
        if (pc.fluidParticles[i].boundary_ == false)
        {
            pc.fluidParticles[i].setDensity(CalculateDensity2(i, pc));
            pc.fluidParticles[i].setPressure(std::max(k * (pc.fluidParticles[i].getDensity() / p0 - 1), 0.0f));
        }
    }
    auto end_pressuredensity = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_pressuredensity = end_pressuredensity - start_pressuredensity;
    // file9 << "pressuredensity: " << elapsed_pressuredensity.count() << " seconds" << std::endl;


    auto start_grid = std::chrono::high_resolution_clock::now();
    uniform_grid.GridClear();
    auto end_grid = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_grid = end_grid - start_grid + elapsed_neighbor;
    // file9 << "Neighbor: " << elapsed_grid.count() << " seconds" << std::endl;

    auto start_acceleration = std::chrono::high_resolution_clock::now();
//#pragma omp parallel for 
    for (int i = 0; i < pc.fluidParticles.size(); i++)
    {
        if (pc.fluidParticles[i].boundary_ == false)
        {


            pc.fluidParticles[i].setNonPAcceleration(CalculateAcceleration2(i, pc));
            pc.fluidParticles[i].setPredictedVelocity();
            pc.fluidParticles[i].setPAcceleration(CalculatePressureAcceleration2(i, pc));

        }
    }
    auto end_acceleration = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_acceleration = end_acceleration - start_acceleration;
    // file9 << "acceleration: " << elapsed_acceleration.count() << " seconds" << std::endl;


    auto start_integration = std::chrono::high_resolution_clock::now();
//#pragma omp parallel for 
    for (int i = 0; i < pc.fluidParticles.size(); i++) {
        if (pc.fluidParticles[i].boundary_ == false) {
            pc.fluidParticles[i].setAcceleration();
            pc.fluidParticles[i].setVelocity();
            pc.fluidParticles[i].setPosition();
        }
        uniform_grid.AddParticle(pc.fluidParticles[i].getPosition(), i);

        if (!pc.fluidParticles[i].boundary_) {
            if (pc.fluidParticles[i].getDensity() < 1.2f) {
                pc.fluidParticles[i].setColor(sf::Color(0, 30, 210));
            }
            else if (pc.fluidParticles[i].getDensity() < 1.205f) {
                pc.fluidParticles[i].setColor(sf::Color(0, 100, 210));
            }
            else if (pc.fluidParticles[i].getDensity() < 1.21f) {
                pc.fluidParticles[i].setColor(sf::Color(0, 150, 210));
            }
            else if (pc.fluidParticles[i].getDensity() < 1.215f) {
                pc.fluidParticles[i].setColor(sf::Color(0, 200, 210));
            }
            else if (pc.fluidParticles[i].getDensity() < 1.22f) {
                pc.fluidParticles[i].setColor(sf::Color(0, 200, 150));
            }
            else if (pc.fluidParticles[i].getDensity() < 1.225f) {
                pc.fluidParticles[i].setColor(sf::Color(0, 200, 50));
            }
            else if (pc.fluidParticles[i].getDensity() < 1.23f) {
                pc.fluidParticles[i].setColor(sf::Color(50, 250, 0));
            }
            else if (pc.fluidParticles[i].getDensity() < 1.235f) {
                pc.fluidParticles[i].setColor(sf::Color(150, 200, 0));
            }
            else if (pc.fluidParticles[i].getDensity() < 1.24f) {
                pc.fluidParticles[i].setColor(sf::Color(250, 100, 0));
            }
            else {
                pc.fluidParticles[i].setColor(sf::Color(200, 0, 0));
            }
        }

    }
    auto end_integration = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_integration = end_integration - start_integration;
    // file9 << "integration: " << elapsed_integration.count() << " seconds" << std::endl;
}



