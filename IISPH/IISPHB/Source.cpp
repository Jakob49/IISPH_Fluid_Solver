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
#include "SPH.h"
#include "Helper.h"

#include "imconfig-SFML.h"
#include "imgui.h"
#include "imgui-SFML.h"
#include "ISPH.h"
#include "RigidBody.h"


// Function to normalize a vector
sf::Vector2f normalize(const sf::Vector2f& vector) {
    float length = std::sqrt(vector.x * vector.x + vector.y * vector.y);
    if (length != 0)
        return sf::Vector2f(vector.x / length, vector.y / length);
    else
        return vector;
}


void drawArrow(sf::RenderWindow& window, const sf::Vector2f& startPos, const sf::Vector2f& forceVector) {
    sf::Vector2f endPos = startPos + forceVector;

   
    sf::Vertex line[] = {
        sf::Vertex(startPos),
        sf::Vertex(endPos)
    };
    window.draw(line, 2, sf::Lines);

    
    float vectorLength = std::sqrt(forceVector.x * forceVector.x + forceVector.y * forceVector.y);

    
    float arrowSize = vectorLength * 0.5f;
    float arrowWidth = arrowSize * 0.3f; 

    
    sf::Vector2f direction = normalize(forceVector);
    sf::Vector2f perpendicular(-direction.y, direction.x);

    
    sf::Vector2f arrowLeft = endPos - direction * arrowSize + perpendicular * arrowWidth;
    sf::Vector2f arrowRight = endPos - direction * arrowSize - perpendicular * arrowWidth;

    sf::Vertex arrow[] = {
        sf::Vertex(endPos),
        sf::Vertex(arrowLeft),
        sf::Vertex(endPos),
        sf::Vertex(arrowRight)
    };
    window.draw(arrow, 4, sf::Lines);
}
float Kernel(sf::Vector2f positionX, sf::Vector2f positionY) {
    float r = std::sqrt(std::pow(positionX.x - positionY.x, 2)
        + std::pow(positionX.y - positionY.y, 2)); // distance of the particles
    float d = r / h;
    float t1 = std::max(1 - d, 0.0f);
    float t2 = std::max(2 - d, 0.0f);
    float w = alpha * (t2 * t2 * t2 - 4 * t1 * t1 * t1);

    return w;
}


void setMass(ParticleContainer& pc, UniformGrid uf) {
    for (size_t i = 0; i < pc.fluidParticles.size(); i++) {
        uf.AddParticle(pc.fluidParticles[i].getPosition(), i);
    }
    for (size_t i = 0; i < pc.fluidParticles.size(); i++) {
        pc.fluidParticles[i].setNeighborIndices(uf.getNeighbors(i, pc.fluidParticles[i].getPosition(), pc));
    }
    float m = 0.0f;
    for (size_t i = 0; i < pc.fluidParticles.size(); i++) {
        if (pc.fluidParticles[i].boundary_ == false) {
            std::vector<size_t> neighbors = pc.fluidParticles[i].getFluidNeighbors();
           
            for (auto& j : neighbors) {
                m = m + Kernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition());
            }
            pc.fluidParticles[i].mass_ =  2 * (p0 / m);
            std::cout << "Particle number " << i << " mass:   " << 2*  p0 /m << std::endl;
            m = 0.0f;
        }

    }
    uf.GridClear();
}


int main() {
    
    sf::ContextSettings settings;
    settings.antialiasingLevel = 8;
    sf::RenderWindow window(sf::VideoMode(1500, 750), "SPH", sf::Style::Default, settings);
    ImGui::SFML::Init(window);
    sf::Clock deltaClock;
    


    std::ofstream file("C:/Users/Jakob/Documents/Studium/BachelorArbeit/plots/data/iterations.txt");
    std::ofstream file2("C:/Users/Jakob/Documents/Studium/BachelorArbeit/plots/data/MaxVelocity.txt");
    std::ofstream file3("C:/Users/Jakob/Documents/Studium/BachelorArbeit/plots/data/DensityError.txt");
    std::ofstream file4("C:/Users/Jakob/Documents/Studium/BachelorArbeit/plots/data/Performance.txt");
    std::ofstream file5("C:/Users/Jakob/Documents/Studium/BachelorArbeit/plots/data/boundary_pressure.txt");
    std::ofstream file6("C:/Users/Jakob/Documents/Studium/BachelorArbeit/plots/data/acceleration_average.txt");
    std::ofstream file7("C:/Users/Jakob/Documents/Studium/BachelorArbeit/plots/data/addaptive_timestep.txt");
    std::ofstream file10("C:/Users/Jakob/Documents/Studium/BachelorArbeit/plots/data/pressure_average.txt");






    ParticleContainer pc;
    SPH sphg;
    ISPH isph;

    UniformGrid uniform_grid;
    std::vector<int> particleIndices;
    for (size_t i = pc.fluidParticles.size() - 16; i < pc.fluidParticles.size(); i++)
    {
        particleIndices.push_back(i);
    }
    
    RigidBody r(pc.movingParticles, pc.outerLayerParticles);

    r.Initialize(pc);
    

    float cfl;
    float maxtime = 0.07f;
    float deltat;

    for (size_t i = 0; i < pc.fluidParticles.size(); i++)
    {
        uniform_grid.AddParticle(pc.fluidParticles[i].getPosition(), i);
    }

    sf::View view(sf::FloatRect(15, 15, 60, 30));
    // sf::View view(sf::FloatRect(15, 75, 180, 90));
    const float moveAmount = 1.0f; 
    const float zoomAmount = 0.01f; 

    view.setCenter(30.f, 30.f);
    view.setSize(300.f, 150.f);


    window.setView(view);
    sf::Texture texture;
    texture.create(window.getSize().x, window.getSize().y);
    float averageDensityError = 0;
    float accDensity = 0;
    float maxVelocity = 0.0f;
    float densityDeviation = 0.0f;
    float pressure_avg = 0.0f;

    
    bool showplot = false;
    bool reset = false;

    bool steps = false;

    bool scenario_one = false;
    bool scenario_two = false;
    bool scenario_three = false;
    bool scenario_four = false;
    bool scenario_five = false;
    bool scenario_six = false;
    bool scenario_seven = false;
    bool scenario_eight = false;
    bool scenario_ten = false;
    bool scenario_eleven = false;
    bool scenario_seven_two = false;
    bool scenario_thirteen = false;
    bool scenario_fourteen = false;
    bool scenario_fifteen = false;

    bool color_selected = false;
    bool color_selected_neighbors = false;
    int sleepy = 0;
    bool pause = true;
    int selected_particle = -1;
    bool moveboundary = false;
    bool right = false;

    bool vis = false;


    int count_compressed = 0;
    int count_boundary = 0;
    float densityerror = 0.0f;
    float acc_pressure = 0.0f;
    float acc_velocity = 0.0f;
    float svel = 0.0f;

    int selectedMethod = 0;

    int step = 0;

    // setMass(pc, uniform_grid);


    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            ImGui::SFML::ProcessEvent(event);
            if (event.type == sf::Event::Closed)
                window.close();
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up)) {
            view.move(0, -moveAmount); // Move view up
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down)) {
            view.move(0, moveAmount); // Move view down
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Left)) {
            view.move(-moveAmount, 0); // Move view left
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right)) {
            view.move(moveAmount, 0); // Move view right
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Z)) {
            view.zoom(1.0f - zoomAmount); // Zoom in
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Z)) {
            view.zoom(1.0f - zoomAmount); // Zoom in
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::U)) {
            view.zoom(1.0f + zoomAmount); // Zoom out
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::N)) {
            pause = false; // Zoom out
        }

        if (event.type == sf::Event::MouseButtonPressed) {
            if (event.mouseButton.button == sf::Mouse::Left) {
                sf::Vector2i mousePos = sf::Mouse::getPosition(window);
                sf::Vector2f worldPos = window.mapPixelToCoords(mousePos, view);
                int s = uniform_grid.GiveClosestParticle(worldPos, pc);
                if ( s != -1) {
                    selected_particle = s;
                }

            }

        }
        
        ImGui::SFML::Update(window, deltaClock.restart());
        
        sf::Vector2u windowSize = window.getSize();

        
        ImGui::SetNextWindowPos(ImVec2(windowSize.x - 300, 0));
        //if (predDensityErrors.size() > 100) {
        //    predDensityErrors.erase(predDensityErrors.begin());
        //}

        //if (DensityErrors.size() > 100) {
        //    predDensityErrors.erase(predDensityErrors.begin());
        //}
        
        if (predDensityErrors.size() > 1000) {
            predDensityErrors.erase(predDensityErrors.begin());
        }

        
        if (DensityErrors.size() > 1000) {
            DensityErrors.erase(DensityErrors.begin());
        }


        ImGui::Begin("Control Panel");

        // Section for reset and scenario checkboxes
        ImGui::Text("General Controls");
        if (ImGui::CollapsingHeader("Scenarios")) {
            ImGui::Checkbox("Reset", &reset);
            ImGui::Checkbox("Breaking Dam", &scenario_one);
            ImGui::Checkbox("Water Tower 80 x 40 double", &scenario_two);
            ImGui::Checkbox("Another Dam 80 x 40 double", &scenario_three);
            ImGui::Checkbox("Dam boundary two", &scenario_four);
            ImGui::Checkbox("Dam 10000", &scenario_five);
            ImGui::Checkbox("rigid", &scenario_six);
            ImGui::Checkbox("tower 20000 ", &scenario_seven);
            ImGui::Checkbox("tower 20000 double layer ", &scenario_seven_two);
            ImGui::Checkbox("Water tower 80 x 40", &scenario_eight);
            ImGui::Checkbox("Breaking Dam 80 x 40", &scenario_ten);
            ImGui::Checkbox("breaking dam 20000 ", &scenario_eleven);
            ImGui::Checkbox("breaking dam 100000", &scenario_thirteen);
            ImGui::Checkbox("tower 30000", &scenario_fourteen);
            ImGui::Checkbox("rigid body 20000", &scenario_fifteen);
        }
        ImGui::Separator();

        ImGui::Checkbox("Pause", &pause);
        ImGui::Checkbox("steps", &steps);
        ImGui::Checkbox("vis", &vis);
        ImGui::Checkbox("Show Plot", &showplot);
        ImGui::Checkbox("Color Selected", &color_selected);
        ImGui::Checkbox("Color Selected Neighbors", &color_selected_neighbors);
        ImGui::Checkbox("move boundary", &moveboundary);
        ImGui::Checkbox("SESPH", &sph);
        

        // Section for sliders
        if (ImGui::CollapsingHeader("Adjust Parameters")) {
            ImGui::Text("Adjust Parameters:");
            ImGui::SliderFloat("Time Step (t)", &t, 0.00001f, 0.1f, "%.6f");
            ImGui::SliderFloat("Omega", &omega, 0.2f, 1.0f);
            ImGui::SliderFloat("viscosity", &v, 0.2f, 1.0f);
            ImGui::SliderFloat("surface tension",&stension, 0.01f, 2.0f);
            ImGui::SliderFloat("stiffness k", &k, 20000, 1000000);
            
            // ImGui::SliderInt("sleepy", &sleepy, 0, 300);
            ImGui::SliderFloat("Pressure Gamma", &boundary_gamma, 0.1f, 1.0f);
            ImGui::SliderFloat("Color Max Value", &maxDefinedVelocity, 2.0f, 100.0f);
            ImGui::Separator();
        }

        ImGui::Separator();
        ImGui::Text("Boundary Handling Method:");
        if (ImGui::RadioButton("Presssure Mirroring", selectedMethod == 0)) {
            selectedMethod = 0;
            mirroring = true;
            mls = false;

        }
        if (ImGui::RadioButton("SPH Pressure Extrapolation", selectedMethod == 1)) {
            mirroring = false;
            mls = false;
            selectedMethod = 1;
        }
        if (ImGui::RadioButton("MLS Pressure Extrapolation", selectedMethod == 2)) {

            mirroring = false;
            mls = true;
            selectedMethod = 2;
        }


        // Section for iteration and time step info
        ImGui::Text("Simulation Info:");
        ImGui::Text("Iteration Count: %d", iter);
        ImGui::Text("Simulation steps: %d", step);
        ImGui::Text("Time Step: %.4f", t);
        // ImGui::Text("Maximal velocity: %.2f", maxVel1);
        ImGui::Separator();

        // Section for particle properties
        ImGui::Text("Particle Properties:");
        ImGui::Text("Particle Number: %d", selected_particle);
        if (selected_particle != -1) {
            std::vector<size_t> neighbors = pc.fluidParticles[selected_particle].getFluidNeighbors();
            ImGui::Text("Particle Neighbors Size: %d", neighbors.size());
            ImGui::Text("Particle Density: %f", pc.fluidParticles[selected_particle].getDensity());
            ImGui::Text("Particle Pressure: %f", pc.fluidParticles[selected_particle].getPressure());
            ImGui::Text("Particle sf: %f", pc.fluidParticles[selected_particle].sf_);
            ImGui::Text("Particle ap: %f", pc.fluidParticles[selected_particle].ap_);
            ImGui::Text("Particle aff: %f", pc.fluidParticles[selected_particle].aff_);
            ImGui::Text("Particle mass: %f", pc.fluidParticles[selected_particle].mass_);


            ImGui::Text("Particle Velocity");
            ImGui::Text("X: %.2f, Y: %.2f", pc.fluidParticles[selected_particle].getVelocity().x, -1.0f * pc.fluidParticles[selected_particle].getVelocity().y);

            ImGui::Text("Fluid Pressure Acceleration");
            ImGui::Text("X: %.2f, Y: %.2f", pc.fluidParticles[selected_particle].particlePressureAcceleration.x, -1.0f * pc.fluidParticles[selected_particle].particlePressureAcceleration.y);

            ImGui::Text("Boundary Pressure Acceleration");
            ImGui::Text("X: %.2f, Y: %.2f", pc.fluidParticles[selected_particle].boundaryPressureAcceleration.x, -1.0f * pc.fluidParticles[selected_particle].boundaryPressureAcceleration.y);
        }

        
        ImGui::End();

        if (showplot) {
            ImGui::Begin("Plots");

            ImGui::Text("Density Errors:");
            ImGui::PlotLines("Pred Density Errors", predDensityErrors.data(), predDensityErrors.size());
            ImGui::PlotLines("Density Errors", DensityErrors.data(), DensityErrors.size());

            ImGui::End();
        }

        
        if (reset) {
            
            pc.fluidParticles.clear();
            pc.setParticles();
            uniform_grid.GridClear();
            for (size_t i = 0; i < pc.fluidParticles.size(); i++)
            {
                uniform_grid.AddParticle(pc.fluidParticles[i].getPosition(), i);
            }
            reset = false; 
        }
        if (scenario_one)
        {
           
            pc.fluidParticles.clear();
            pc.setParticles1();
            uniform_grid.GridClear();
            for (size_t i = 0; i < pc.fluidParticles.size(); i++)
            {
                uniform_grid.AddParticle(pc.fluidParticles[i].getPosition(), i);
            }
            scenario_one = false;
        }
        else if (scenario_two)
        {
            
            pc.fluidParticles.clear();
            pc.setParticles2();
            uniform_grid.GridClear();
            for (size_t i = 0; i < pc.fluidParticles.size(); i++)
            {
                uniform_grid.AddParticle(pc.fluidParticles[i].getPosition(), i);
            }
            scenario_two = false;
            view.setCenter(40.f, 30.f);
            view.setSize(100.f, 50.f);
        }
        else if (scenario_three)
        {
            
            pc.fluidParticles.clear();
            pc.setParticles3();
            uniform_grid.GridClear();
            for (size_t i = 0; i < pc.fluidParticles.size(); i++)
            {
                uniform_grid.AddParticle(pc.fluidParticles[i].getPosition(), i);
            }
            scenario_three = false;
            view.setCenter(40.f, 30.f);
            view.setSize(100.f, 50.f);
        }
        else if (scenario_four)
        {
            pc.fluidParticles.clear();
            pc.setParticles4();
            uniform_grid.GridClear();
            for (size_t i = 0; i < pc.fluidParticles.size(); i++)
            {
                uniform_grid.AddParticle(pc.fluidParticles[i].getPosition(), i);
            }
            scenario_four = false;
            view.setCenter(40.f, 30.f);
            view.setSize(100.f, 50.f);
        }

        else if (scenario_five)
        {
            pc.fluidParticles.clear();
            pc.setParticles5();
            uniform_grid.GridClear();
            for (size_t i = 0; i < pc.fluidParticles.size(); i++)
            {
                uniform_grid.AddParticle(pc.fluidParticles[i].getPosition(), i);
            }
            scenario_five = false;
            view.setCenter(220.f, 310.f);
            view.setSize(640.f, 320.f);
        }
        else if (scenario_six)
        {
            pc.fluidParticles.clear();
            pc.setParticles6();
            uniform_grid.GridClear();
            for (size_t i = 0; i < pc.fluidParticles.size(); i++)
            {
                uniform_grid.AddParticle(pc.fluidParticles[i].getPosition(), i);
            }
            scenario_six = false;
            view.setCenter(220.f, 310.f);
            view.setSize(640.f, 320.f);
        }
        else if (scenario_seven)
        {
            pc.fluidParticles.clear();
            pc.setParticles7();
            uniform_grid.GridClear();
            for (size_t i = 0; i < pc.fluidParticles.size(); i++)
            {
                uniform_grid.AddParticle(pc.fluidParticles[i].getPosition(), i);
            }
            scenario_seven = false;
            view.setCenter(30.f, 20.f);
            view.setSize(250.f, 125.f);
        }

        else if (scenario_eight)
        {
            pc.fluidParticles.clear();
            pc.setParticles8();
            uniform_grid.GridClear();
            for (size_t i = 0; i < pc.fluidParticles.size(); i++)
            {
                uniform_grid.AddParticle(pc.fluidParticles[i].getPosition(), i);
            }
            scenario_eight = false;
            view.setCenter(30.f, 20.f);
            view.setSize(250.f, 125.f);
        }
        else if (scenario_ten)
        {
            pc.fluidParticles.clear();
            pc.setParticles10();
            uniform_grid.GridClear();
            for (size_t i = 0; i < pc.fluidParticles.size(); i++)
            {
                uniform_grid.AddParticle(pc.fluidParticles[i].getPosition(), i);
            }
            scenario_ten = false;
            view.setCenter(30.f, 20.f);
            view.setSize(250.f, 125.f);
        }
        else if (scenario_eleven)
        {
            pc.fluidParticles.clear();
            pc.setParticles11();
            uniform_grid.GridClear();
            for (size_t i = 0; i < pc.fluidParticles.size(); i++)
            {
                uniform_grid.AddParticle(pc.fluidParticles[i].getPosition(), i);
            }
            scenario_eleven = false;
            view.setCenter(30.f, 20.f);
            view.setSize(250.f, 125.f);
        }
        else if (scenario_seven_two)
        {
            pc.fluidParticles.clear();
            pc.setParticles7_2();
            uniform_grid.GridClear();
            for (size_t i = 0; i < pc.fluidParticles.size(); i++)
            {
                uniform_grid.AddParticle(pc.fluidParticles[i].getPosition(), i);
            }
            scenario_seven_two = false;
            view.setCenter(160.f, 250.f);
            view.setSize(640.f, 320.f);
        }
        else if (scenario_thirteen)
        {
            pc.fluidParticles.clear();
            pc.setParticles13();
            uniform_grid.GridClear();
            for (size_t i = 0; i < pc.fluidParticles.size(); i++)
            {
                uniform_grid.AddParticle(pc.fluidParticles[i].getPosition(), i);
            }
            scenario_thirteen = false;
            view.setCenter(160.f, 250.f);
            view.setSize(640.f, 320.f);
        }
        else if (scenario_fourteen)
        {
            pc.fluidParticles.clear();
            pc.setParticles14();
            uniform_grid.GridClear();
            for (size_t i = 0; i < pc.fluidParticles.size(); i++)
            {
                uniform_grid.AddParticle(pc.fluidParticles[i].getPosition(), i);
            }
            scenario_fourteen = false;
            view.setCenter(160.f, 250.f);
            view.setSize(640.f, 320.f);
            }
        else if (scenario_fifteen)
        {
            pc.fluidParticles.clear();
            pc.setParticles15();
            uniform_grid.GridClear();
            for (size_t i = 0; i < pc.fluidParticles.size(); i++)
            {
                uniform_grid.AddParticle(pc.fluidParticles[i].getPosition(), i);
            }
            scenario_fifteen = false;
            view.setCenter(160.f, 250.f);
            view.setSize(640.f, 320.f);
            }


        
       
        if (pause == false && sph) {
            
            
            auto start_sph = std::chrono::high_resolution_clock::now();
            sphg.Update(pc, uniform_grid);
            auto end_sph = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed_solver = end_sph - start_sph;
            // file9 << "overall time, " << elapsed_solver.count() << " seconds" << t <<  std::endl;
            file9 <<  elapsed_solver.count() << ", " << t << std::endl;

        }
        if (pause == false && sph == false) {
            auto start_solver = std::chrono::high_resolution_clock::now();
            isph.Update(pc, uniform_grid);
            auto end_solver = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed_solver = end_solver - start_solver;
            // file8 << "overall time, " << elapsed_solver.count() << " seconds" << t << std::endl;
            file8 << elapsed_solver.count() << ", " << t << std::endl;
        }

        if (pc.rigidbody) {
            for (auto& i : r.particleIndices) {
                pc.fluidParticles[i].setNeighborIndices(uniform_grid.getNeighbors(i, pc.fluidParticles[i].getPosition(), pc));
                if (vis) {
                    pc.fluidParticles[i].setColor(sf::Color::Transparent);
                }
            }
            if (vis) {
                for (auto& k : r.outer) {
                    int i = r.particleIndices[k];
                    pc.fluidParticles[i].setColor(sf::Color(75, 75, 75));
                }
            }

        }
        
        if (pause == false) {

            if (pc.rigidbody) {
                r.Update(pc);
            }
        }

        if (moveboundary) {
            int size = pc.fluidParticles.size() - 1;
            const float velocity = -2.0f;
            const float leftBoundary = 50.0f;
            const float rightBoundary = 210.0f;

            // Check if the last particle in the list has reached the boundary
            if ((pc.fluidParticles[size - 1].getPosition().x > rightBoundary && right) ||
                (pc.fluidParticles[size - 1].getPosition().x < leftBoundary && !right)) {
                // Reverse direction
                c = -c;
                right = !right;
            }

            
            for (auto& j : pc.movingParticles){
          
                sf::Vector2f position = pc.fluidParticles[j].getPosition();
                position.x += c * velocity * t; 
                pc.fluidParticles[j].setBoundaryPosition(position);
            }
        }
        window.clear(sf::Color::White);

        if (color_selected && selected_particle != -1) {
            pc.fluidParticles[selected_particle].setColor(sf::Color(sf::Color::Green));
        }
        if (color_selected_neighbors && selected_particle != -1) {
            std::vector<size_t> neighbors = pc.fluidParticles[selected_particle].getFluidNeighbors();
            for (auto& i : neighbors) {
                pc.fluidParticles[i].setColor(sf::Color(sf::Color::Red));
            }
            pc.fluidParticles[selected_particle].setColor(sf::Color(sf::Color::Green));
        }
        

 



        for (size_t i = 0; i < pc.fluidParticles.size(); i++) {

            window.draw(pc.fluidParticles[i].getShape());

        }
        //if (pc.rigidbody) {
        //    // temporaly to debug rigid body
        //    sf::CircleShape centerOfMassPoint(0.1f); 
        //    centerOfMassPoint.setFillColor(sf::Color::Red);
        //    centerOfMassPoint.setPosition(r.centerOfMass);

        //    for (const auto& relPos : r.relativePosition) {
        //        sf::CircleShape point(0.1f);  
        //        point.setFillColor(sf::Color::Blue);
        //        point.setPosition(r.centerOfMass + relPos);  // Add center of mass to get actual position in window
        //        window.draw(point);
        //    }
        //    window.draw(centerOfMassPoint);
        //    // Draw an arrow for each force vector
        //    for (int i = 0; i < r.particleIndices.size(); i++) {
        //        drawArrow(window, sf::Vector2f(pc.fluidParticles[r.particleIndices[i]].getPosition().x + h / 2 - 0.05f, pc.fluidParticles[r.particleIndices[i]].getPosition().y + h / 2 - 0.05f), r.forces[i]);
        //    }
        //    for (const auto& k : r.outer) {
        //        int i = r.particleIndices[k];
        //        sf::CircleShape point(0.2f); 
        //        point.setFillColor(sf::Color::Magenta);
        //        point.setPosition(pc.fluidParticles[i].getPosition());
        //        window.draw(point);
        //        
        //    }

        //}

        window.setView(view);

        if (steps) {
            pause = true;
        }
        



        ImGui::SFML::Render(window);
 



        window.display();
        
        texture.update(window);

        //if (step % 10 == 0) {
        //    sf::Image screenshot = texture.copyToImage();
        //    screenshot.saveToFile("C:/Users/Jakob/Documents/Studium/BachelorArbeit/animation/" + std::to_string(step) + ".png");
        //}
        //sf::Image screenshot = texture.copyToImage();
        //screenshot.saveToFile("C:/Users/Jakob/Documents/Studium/BachelorArbeit/animation/" + std::to_string(step) + ".png");


        for (size_t i = 0; i < pc.fluidParticles.size(); i++) {
            if (pc.fluidParticles[i].boundary_ == false) {
                if (pc.fluidParticles[i].getDensity() > p0) {
                    densityerror += std::abs(pc.fluidParticles[i].getDensity() - p0) ;
                    count_compressed++;
                }
                //if (pc.fluidParticles[i].ap_ >= pc.fluidParticles[i].sf_) {
                //    densityerror += pc.fluidParticles[i].ap_ - pc.fluidParticles[i].sf_;
                //    count_compressed++;
                //}
                //densityerror += (pc.fluidParticles[i].getDensity() - p0);
                //count_compressed++;
                acc_velocity = acc_velocity + Helper::MagniVector(pc.fluidParticles[i].getVelocity());
                float current_velocity = Helper::MagniVector(pc.fluidParticles[i].getVelocity());
                if (current_velocity > svel) {
                    svel = current_velocity;
                }

                pressure_avg = pressure_avg + pc.fluidParticles[i].getPressure();

            }
            else if (pc.fluidParticles[i].boundary_ &&  pc.fluidParticles[i].getPressure() > 0.0f)  {

                acc_pressure = acc_pressure + pc.fluidParticles[i].getPressure();
                count_boundary++;
            }
        }


        // Check if count_compressed is greater than 0 to avoid division by zero
        if (count_compressed > 0) {
            densityerror = densityerror / count_compressed / p0 * 100;
        }
        else {
            densityerror = 0.0f;  
        }

        if (count_boundary > 0) {
            acc_pressure = acc_pressure / count_boundary;
        }
        else {
            acc_pressure = 0.0f;  
        }
        acc_velocity = acc_velocity / pc.fluidNumber_;
        pressure_avg = pressure_avg / pc.fluidNumber_;


        if (step % 400 == 0 && step != 0) {
            omega += 0.1;
            scenario_eight = true;
        }
        

        count_compressed = 0;
        count_boundary = 0;
        bool cfl_vio = false;

        if (t >= 0.4f *  h/ svel ) {
            cfl_vio = true;
        }

        //texture.update(window);
        //sf::Image screenshot = texture.copyToImage();
        //screenshot.saveToFile("C:/Users/Jakob/Documents/Studium/BachelorArbeit/IISPHB/animation/" + std::to_string(step) + ".png");
        if (pause != true) {
            step++;

            // file << iter << ", " << step << ", " << t << std::endl; 
            file << iter << ", " << step << ", " << omega << std::endl;
            file2 << svel << ", " << step << ", " << t << ", cfl " << cfl_vio <<std::endl;
            // file2 << maxVel1 << ", " << step << ", " << omega << std::endl;
            file3 << densityerror << ", " << step << ", t =" << t << ", k = " << k << std::endl;
            file5 << acc_pressure << ", " << step << ", " << t << std::endl;
            file6 << acc_velocity << ", " << step << ", " << t << std::endl;
            // file7 << t << ", " << step << ", " << omega << std::endl;
            // file10 << pressure_avg << ", " << step << t << std::endl;
            file10 << pressure_avg << ", " << step << ", " << omega << std::endl;
        }
        
        densityerror = 0.0f;
        acc_pressure = 0.0f;
        acc_velocity = 0.0f;
        pressure_avg = 0.0f;
        maxVel1 = 0.0f;
        svel = 0.0f;


    }
    ImGui::SFML::Shutdown();

    file.close();
    file2.close();
    file3.close();
    file4.close();
    file5.close();
    file7.close();
    file8.close();

    return 0;
}