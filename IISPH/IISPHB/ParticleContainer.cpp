#include "ParticleContainer.h"

#include "Particle.h"
#include <iostream>
#define _USE_MATH_DEFINES
#include<math.h>
#include <cmath>

#include "Globals.h"
ParticleContainer::ParticleContainer()
{
    
     //setParticles15();
    //setParticles6();
    
    // setParticles11();
    //setParticles2();


    // setParticles8();

    // setParticles10();

    // test for comparing physical states
     // setParticles16();


    //for SESPH
    // setParticles17();


    // setParticles();

    // setParticles7_2();

    // for IIPSH
    // setParticles11();

    // rigid body
    //setParticles15();


    // 80 x 40
    setParticles8();


    // setParticles11();

    // 100 000
    // setParticles13();

    //setParticles1();
    //setParticlesBowl();

    // adjustable height
    // setParticles18();

    //setParticles8();

}
void ParticleContainer::setParticles() {

    // sf::View view(sf::FloatRect(15, 15, 60, 30));

    int x = 20;
    int y = 10;

    fluidNumber_ = x * y;
    rigidbody = false;
    for (int i = 0; i < x; i++)
    {
        for (int j = 0; j < y; j++)
        {
            sf::Vector2f position(29 * h - i * h, 40 * h - j * h);
            Particle particle(position, false);
            fluidParticles.push_back(particle);
        }
    }


    // boundary
    // left
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < 14; j++)
        {


            sf::Vector2f position(9 * h - i * h, 40 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }



    // right
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < 14; j++)
        {

            sf::Vector2f position(30 * h + i * h, 40 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // bottom
    for (int i = 0; i < 22; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            sf::Vector2f position(i * h + 9 * h, 41 * h + j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }
    // bottom
    for (int i = 0; i < 22; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            sf::Vector2f position(i * h + 9 * h, 26 * h + j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

}

void ParticleContainer::setParticles1() {

    int x = 10;
    int y = 20;

    fluidNumber_ = x * y;
    rigidbody = false;

    for (int i = 0; i < x; i++)
    {
        for (int j = 0; j < y; j++)
        {
            sf::Vector2f position(21 * h + i * h, 40 * h - j * h);
            Particle particle(position, false);
            fluidParticles.push_back(particle);
        }
    }


    // boundary
    // left
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < 20; j++)
        {


            sf::Vector2f position(20 * h - i * h, 40 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }



    // right
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < 20; j++)
        {

            sf::Vector2f position(51 * h + i * h, 40 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // bottom
    for (int i = 0; i < 32; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            sf::Vector2f position(i * h + 20 * h, 41 * h + j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // top
    for (int i = 0; i < 32; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            sf::Vector2f position(i * h + 20 * h, 20 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }
}


// "Water Tower"
void ParticleContainer::setParticles2() {
    int x = 40;
    int y = 80;

    fluidNumber_ = x * y;
    rigidbody = false;

    
    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            sf::Vector2f position(21 * h + i * h, 101 * h - j * h);  
            Particle particle(position, false);
            fluidParticles.push_back(particle);
        }
    }

    // boundary particles
    // left boundary
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 91; j++) {
            sf::Vector2f position(20 * h - i * h, 101 * h - j * h);  
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // right boundary
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 91; j++) {
            sf::Vector2f position(61 * h + i * h, 101 * h - j * h); 
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // bottom boundary
    for (int i = 0; i < 46; i++) {
        for (int j = 0; j < 3; j++) {
            sf::Vector2f position(i * h + 18 * h, 102 * h + j * h);  
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // top boundary 
    for (int i = 0; i < 46; i++) {
        for (int j = 0; j < 3; j++) {
            sf::Vector2f position(i * h + 18 * h, 10 * h - j * h);  
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

}


void ParticleContainer::setParticles3() {

    int x = 40;
    int y = 80;

    fluidNumber_ = x * y;
    rigidbody = false;

    
    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            sf::Vector2f position(21 * h + i * h, 101 * h - j * h);  
            Particle particle(position, false);
            fluidParticles.push_back(particle);
        }
    }

    // boundary particles
    // left boundary
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 91; j++) {
            sf::Vector2f position(20 * h - i * h, 101 * h - j * h);  
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // right boundary
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 91; j++) {
            sf::Vector2f position(111 * h + i * h, 101 * h - j * h);  
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // bottom boundary
    for (int i = 0; i < 94; i++) {
        for (int j = 0; j < 2; j++) {
            sf::Vector2f position(i * h + 19 * h, 102 * h + j * h);  
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // top boundary 
    for (int i = 0; i < 94; i++) {
        for (int j = 0; j < 2; j++) {
            sf::Vector2f position(i * h + 19 * h, 10 * h - j * h);  
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

}

void ParticleContainer::setParticles4() {

    int x = 10;
    int y = 20;

    fluidNumber_ = x * y;
    rigidbody = false;

    for (int i = 0; i < x; i++)
    {
        for (int j = 0; j < y; j++)
        {
            sf::Vector2f position(21 * h + i * h, 40 * h - j * h);
            Particle particle(position, false);
            fluidParticles.push_back(particle);
        }
    }


    // boundary
    // left
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 20; j++)
        {


            sf::Vector2f position(20 * h - i * h, 40 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }



    // right
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 20; j++)
        {

            sf::Vector2f position(51 * h + i * h, 40 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // bottom
    for (int i = 0; i < 34; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            sf::Vector2f position(i * h + 19 * h, 41 * h + j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // top
    for (int i = 0; i < 34; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            sf::Vector2f position(i * h + 19 * h, 20 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }




}
void ParticleContainer::moveBoundaryParticles() {
    

    int size = fluidParticles.size() - 1;
    sf::Vector2f center(120.0f * h, 279.0f * h);
    const float angularSpeed = 10.0f;
    for (int i = 0; i++; i < 50) {
        sf::Vector2f r = fluidParticles[size - i].getPosition() - center;

        
        sf::Vector2f tangentialVelocity(-r.y, r.x);

        
        float length = std::sqrt(tangentialVelocity.x * tangentialVelocity.x + tangentialVelocity.y * tangentialVelocity.y);
        if (length != 0.0f) {
            tangentialVelocity.x /= length;
            tangentialVelocity.y /= length;
        }

        
        tangentialVelocity.x *= angularSpeed * length;
        tangentialVelocity.y *= angularSpeed * length;
        fluidParticles[size - i].setBoundaryVelocity(tangentialVelocity);

        // fluidParticles[i].setBoundaryPosition(fluidParticles[i].getPosition() + t * tangentialVelocity);
        
    }

    std::cout << "particle number " << size << " position  x" << fluidParticles[size].getPosition().x << " position  x" << fluidParticles[size].getPosition().y << std::endl;

}

void ParticleContainer::setParticles5() {

    int x = 100;
    int y = 100;

    fluidNumber_ = x * y;
    rigidbody = false;

    for (int i = 0; i < x; i++)
    {
        for (int j = 0; j < y; j++)
        {
            sf::Vector2f position(21 * h + i * h, 340 * h - j * h);
            Particle particle(position, false);
            fluidParticles.push_back(particle);
        }
    }


    // boundary
    // left
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 120; j++)
        {


            sf::Vector2f position(20 * h - i * h, 340 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }



    // right
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 120; j++)
        {

            sf::Vector2f position(253 * h + i * h, 340 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // bottom
    for (int i = 0; i < 236; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            sf::Vector2f position(i * h + 19 * h, 341 * h + j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // top
    for (int i = 0; i < 236; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            sf::Vector2f position(i * h + 19 * h, 220 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }
    // moving boundary
    int count = 0;
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < 50; j++)
        {
            movingParticles.push_back(fluidParticles.size() + count);
            sf::Vector2f position(i * h + 140 * h, 341 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }
}



void ParticleContainer::setParticles6() {

    int x = 50;
    int y = 110;

    fluidNumber_ = x * y;
    rigidbody = false;

    for (int i = 0; i < x; i++)
    {
        for (int j = 0; j < y; j++)
        {
            sf::Vector2f position(21 * h + i * h, 340 * h - j * h);
            Particle particle(position, false);
            fluidParticles.push_back(particle);
        }
    }


    // boundary
    // left
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 120; j++)
        {


            sf::Vector2f position(20 * h - i * h, 340 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }



    // right
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 120; j++)
        {

            sf::Vector2f position(71 * h + i * h, 340 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // bottom
    for (int i = 0; i < 53; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            sf::Vector2f position(i * h + 19 * h, 341 * h + j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // top
    for (int i = 0; i < 53; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            sf::Vector2f position(i * h + 19 * h, 220 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }
    //rigid body

    int count = 0;
    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            movingParticles.push_back(fluidParticles.size());
            sf::Vector2f position(i * h + 40 * h, 225 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
            if (i == 0 || i == 7 || j == 0 || j == 3)
            {
                outerLayerParticles.push_back(count); // Store the index of the outer particle
                std::cout << "Adding particle at index: " << count << std::endl;

            }
            count++;
        }
    }
    rigidbody = true;


}



// for the thesis the water tower for iisph
void ParticleContainer::setParticles7() {

    int x = 100;
    int y = 200;



    fluidNumber_ = x * y;
    rigidbody = false;

    for (int i = 0; i < x; i++)
    {
        for (int j = 0; j < y; j++)
        {
            sf::Vector2f position(21 * h + i * h, 340 * h - j * h);
            Particle particle(position, false);
            fluidParticles.push_back(particle);
        }
    }


    // boundary
    // left
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < 220; j++)
        {
            sf::Vector2f position(20 * h - i * h, 340 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }



    // right
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < 220; j++)
        {

            sf::Vector2f position(121 * h + i * h, 340 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // bottom
    for (int i = 0; i < 102; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            sf::Vector2f position(i * h + 20 * h, 341 * h + j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // top
    for (int i = 0; i < 102; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            sf::Vector2f position(i * h + 20 * h, 120 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

}


// for comparion of sesph and iisph
void ParticleContainer::setParticles8() {

    int x = 40;
    int y = 80;

    fluidNumber_ = x * y;
    rigidbody = false;

    
    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            sf::Vector2f position(21 * h + i * h, 101 * h - j * h);  
            Particle particle(position, false);
            fluidParticles.push_back(particle);
        }
    }

    // boundary particles
    // left boundary
    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < 92; j++) {
            sf::Vector2f position(20 * h - i * h, 101 * h - j * h);  
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // right boundary
    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < 92; j++) {
            sf::Vector2f position(61 * h + i * h, 101 * h - j * h);  
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // bottom boundary
    for (int i = 0; i < 42; i++) {
        for (int j = 0; j < 1; j++) {
            sf::Vector2f position(i * h + 20 * h, 102 * h + j * h);  
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // top boundary 
    for (int i = 0; i < 42; i++) {
        for (int j = 0; j < 1; j++) {
            sf::Vector2f position(i * h + 20 * h, 10 * h - j * h); 
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

}

// testing the boundary
void ParticleContainer::setParticles9() {

    int x = 98;
    int y = 198;



    fluidNumber_ = x * y;
    rigidbody = false;

    for (int i = 0; i < x; i++)
    {
        for (int j = 0; j < y; j++)
        {
            sf::Vector2f position(22 * h + i * h, 339 * h - j * h);
            Particle particle(position, false);
            fluidParticles.push_back(particle);
        }
    }


    // boundary
    // left
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < 220; j++)
        {
            sf::Vector2f position(20 * h - i * h, 340 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }



    // right
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < 220; j++)
        {

            sf::Vector2f position(121 * h + i * h, 340 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // bottom
    for (int i = 0; i < 104; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            sf::Vector2f position(i * h + 19 * h, 341 * h + j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // top
    for (int i = 0; i < 104; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            sf::Vector2f position(i * h + 19 * h, 120 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

}


// for comparion of sesph and iisph
void ParticleContainer::setParticles10() {

    int x = 40;
    int y = 80;

    fluidNumber_ = x * y;
    rigidbody = false;

    
    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            sf::Vector2f position(21 * h + i * h, 101 * h - j * h);  
            Particle particle(position, false);
            fluidParticles.push_back(particle);
        }
    }

    // boundary particles
    // left boundary
    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < 92; j++) {
            sf::Vector2f position(20 * h - i * h, 101 * h - j * h); 
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // right boundary
    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < 92; j++) {
            sf::Vector2f position(101 * h + i * h, 101 * h - j * h);  
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // bottom boundary
    for (int i = 0; i < 82; i++) {
        for (int j = 0; j < 1; j++) {
            sf::Vector2f position(i * h + 20 * h, 102 * h + j * h);  
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // top boundary 
    for (int i = 0; i < 82; i++) {
        for (int j = 0; j < 1; j++) {
            sf::Vector2f position(i * h + 20 * h, 10 * h - j * h);  
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

}


// breaking damm 20000
void ParticleContainer::setParticles11() {

    int x = 100;
    int y = 200;


    fluidNumber_ = x * y;
    rigidbody = false;

    for (int i = 0; i < x; i++)
    {
        for (int j = 0; j < y; j++)
        {
            sf::Vector2f position(21 * h + i * h, 340 * h - j * h);
            Particle particle(position, false);
            fluidParticles.push_back(particle);
        }
    }


    // boundary
    // left
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < 220; j++)
        {
            sf::Vector2f position(20 * h - i * h, 340 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }



    // right
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < 220; j++)
        {

            sf::Vector2f position(221 * h + i * h, 340 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // bottom
    for (int i = 0; i < 202; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            sf::Vector2f position(i * h + 20 * h, 341 * h + j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // top
    for (int i = 0; i < 202; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            sf::Vector2f position(i * h + 20 * h, 120 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

}

// for the thesis the water tower for iisph
void ParticleContainer::setParticles7_2() {

    int x = 100;
    int y = 200;



    fluidNumber_ = x * y;
    rigidbody = false;

    for (int i = 0; i < x; i++)
    {
        for (int j = 0; j < y; j++)
        {
            sf::Vector2f position(21 * h + i * h, 340 * h - j * h);
            Particle particle(position, false);
            fluidParticles.push_back(particle);
        }
    }


    // boundary
    // left
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 220; j++)
        {
            sf::Vector2f position(20 * h - i * h, 340 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }



    // right
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 220; j++)
        {

            sf::Vector2f position(121 * h + i * h, 340 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // bottom
    for (int i = 0; i < 104; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            sf::Vector2f position(i * h + 19 * h, 341 * h + j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // top
    for (int i = 0; i < 104; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            sf::Vector2f position(i * h + 19 * h, 120 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

}


// breaking damm 100000
void ParticleContainer::setParticles13() {

    int x = 250;
    int y = 400;


    fluidNumber_ = x * y;
    rigidbody = false;

    for (int i = 0; i < x; i++)
    {
        for (int j = 0; j < y; j++)
        {
            sf::Vector2f position(141 * h + i * h, 540 * h - j * h);
            Particle particle(position, false);
            fluidParticles.push_back(particle);
        }
    }


    // boundary
    // left
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < 420; j++)
        {
            sf::Vector2f position(20 * h - i * h, 540 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }



    // right
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < 420; j++)
        {

            sf::Vector2f position(501 * h + i * h, 540 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // bottom
    for (int i = 0; i < 482; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            sf::Vector2f position(i * h + 20 * h, 541 * h + j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // top
    for (int i = 0; i < 482; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            sf::Vector2f position(i * h + 20 * h, 120 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

}



// for the thesis the water tower for iisph
void ParticleContainer::setParticles14() {

    int x = 100;
    int y = 400;



    fluidNumber_ = x * y;
    rigidbody = false;

    for (int i = 0; i < x; i++)
    {
        for (int j = 0; j < y; j++)
        {
            sf::Vector2f position(21 * h + i * h, 540 * h - j * h);
            Particle particle(position, false);
            fluidParticles.push_back(particle);
        }
    }


    // boundary
    // left
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < 420; j++)
        {
            sf::Vector2f position(20 * h - i * h, 540 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }



    // right
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < 420; j++)
        {

            sf::Vector2f position(121 * h + i * h, 540 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // bottom
    for (int i = 0; i < 102; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            sf::Vector2f position(i * h + 20 * h, 541 * h + j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // top
    for (int i = 0; i < 102; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            sf::Vector2f position(i * h + 20 * h, 120 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

}

void ParticleContainer::setParticles15() {

    int x = 100;
    int y = 200;



    fluidNumber_ = x * y;
    rigidbody = false;

    for (int i = 0; i < x; i++)
    {
        for (int j = 0; j < y; j++)
        {
            sf::Vector2f position(21 * h + i * h, 340 * h - j * h);
            Particle particle(position, false);
            fluidParticles.push_back(particle);
        }
    }


    // boundary
    // left
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < 260; j++)
        {
            sf::Vector2f position(20 * h - i * h, 340 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }



    // right
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < 260; j++)
        {

            sf::Vector2f position(121 * h + i * h, 340 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // bottom
    for (int i = 0; i < 102; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            sf::Vector2f position(i * h + 20 * h, 341 * h + j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // top
    for (int i = 0; i < 102; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            sf::Vector2f position(i * h + 20 * h, 80 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }
    //rigid body

    int count = 0;
    for (int i = 0; i < 36; i++)
    {
        for (int j = 0; j < 22; j++)
        {
            movingParticles.push_back(fluidParticles.size());
            // sf::Vector2f position(i * h + 55 * h, 125 * h - j * h);
            sf::Vector2f position(i * h + 55 * h, 130 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
            if (i == 0 || i == 35 || j == 0 || j == 21)
            {
                outerLayerParticles.push_back(count); // Store the index of the outer particle
                std::cout << "Adding particle at index: " << count << std::endl;
               
            }
            count++;

        }
    }
    std::cout << " size of fluid particles " << fluidParticles.size() << std::endl;
    rigidbody = true;

}

// "Water Tower"
void ParticleContainer::setParticles16() {
    int x = 40;
    int y = 80;

    fluidNumber_ = x * y;
    rigidbody = false;

    
    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            sf::Vector2f position(21 * h + i * h, 101 * h - j * h);  
            Particle particle(position, false);
            fluidParticles.push_back(particle);
        }
    }

    // boundary particles
    // left boundary
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 226; j++) {
            sf::Vector2f position(20 * h - i * h, 236 * h - j * h);  
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // right boundary
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 226; j++) {
            sf::Vector2f position(61 * h + i * h, 236 * h - j * h);  
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // left slope
    float a = 0.5f;
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 3; j++) {
            sf::Vector2f position(i * h + 21 * h, 102 * h + j * h + a);  
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
        a = a + 0.5f;
    }

    // right slope
    float b = 0.5f;
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 3; j++) {
            sf::Vector2f position( - i * h + 60 * h, 102 * h + j * h + b);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
        b = b + 0.5f;
    }


    // top boundary 
    for (int i = 0; i < 46; i++) {
        for (int j = 0; j < 3; j++) {
            sf::Vector2f position(i * h + 18 * h, 10 * h - j * h);  
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // left triangle
    float c = 0.0f;
    for (int i = 0; i < 18; i++) {
        for (int j = 0; j < 3; j++) {
            sf::Vector2f position(i * h + 41 * h, 122 * h + j * h + c);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
        c = c + 0.8f;
    }

    // right triangle
    float d = 0.0f;
    for (int i = 0; i < 18; i++) {
        for (int j = 0; j < 3; j++) {
            sf::Vector2f position(-i * h + 40 * h, 122 * h + j * h + d);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
        d = d + 0.8f;
    }

    // left slope
    float e = 0.5f;
    for (int i = 0; i < 18; i++) {
        for (int j = 0; j < 3; j++) {
            sf::Vector2f position(i * h + 21 * h, 142 * h + j * h + e);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
        e = e + 0.5f;
    }

    // right slope
    float f = 0.5f;
    for (int i = 0; i < 18; i++) {
        for (int j = 0; j < 3; j++) {
            sf::Vector2f position(-i * h + 60 * h, 142 * h + j * h + f);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
        f = f + 0.5f;
    }

    // bottom boundary
    for (int i = 0; i < 46; i++) {
        for (int j = 0; j < 3; j++) {
            sf::Vector2f position(i * h + 18 * h, 236 * h + j * h);  
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }



}


void ParticleContainer::setParticles17() {


    int x = 100;
    int y = 200;


    fluidNumber_ = x * y;
    rigidbody = false;

    for (int i = 0; i < x; i++)
    {
        for (int j = 0; j < y; j++)
        {
            sf::Vector2f position(21 * h + i * h, 340 * h - j * h);
            Particle particle(position, false);
            fluidParticles.push_back(particle);
        }
    }


    // boundary
    // left
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 220; j++)
        {
            sf::Vector2f position(20 * h - i * h, 340 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }



    // right
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 220; j++)
        {

            sf::Vector2f position(221 * h + i * h, 340 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // bottom
    for (int i = 0; i < 204; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            sf::Vector2f position(i * h + 19 * h, 341 * h + j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // top
    for (int i = 0; i < 204; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            sf::Vector2f position(i * h + 19 * h, 120 * h - j * h);
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

}


// breaking damm 20000
void ParticleContainer::setParticles18() {


    int x = 80;
    int y = 250;

    fluidNumber_ = x * y;
    rigidbody = false;

   
    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            sf::Vector2f position(21 * h + i * h, 361 * h - j * h); 
            Particle particle(position, false);
            fluidParticles.push_back(particle);
        }
    }

    // boundary particles
    // left boundary
    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < 352; j++) {
            sf::Vector2f position(20 * h - i * h, 361 * h - j * h);  
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // right boundary
    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < 352; j++) {
            sf::Vector2f position(101 * h + i * h, 361 * h - j * h);  
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // bottom boundary
    for (int i = 0; i < 82; i++) {
        for (int j = 0; j < 1; j++) {
            sf::Vector2f position(i * h + 20 * h, 362 * h + j * h);  
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

    // top boundary 
    for (int i = 0; i < 82; i++) {
        for (int j = 0; j < 1; j++) {
            sf::Vector2f position(i * h + 20 * h, 10 * h - j * h); 
            Particle particle(position, true);
            fluidParticles.push_back(particle);
        }
    }

}

