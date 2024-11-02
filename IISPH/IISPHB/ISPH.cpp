#include "ISPH.h"



float CalculateDistance(sf::Vector2f a, sf::Vector2f b)
{
    return std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2));
}


float MagnitudeVector(sf::Vector2f vec) {
    return std::sqrt(vec.x * vec.x + vec.y * vec.y);
}


std::vector<size_t> FindNeighborIndices(Particle& particle, std::vector<Particle>& particles)
{
    std::vector<size_t> neighborIndices;

    const float epsilon = 0.01f;

    for (size_t i = 0; i < particles.size(); ++i) {
        float distance = CalculateDistance(particle.getPosition(), particles[i].getPosition());
        if (distance <= 2 * h + epsilon) {
            neighborIndices.push_back(i);
        }
    }
    return neighborIndices;
}

float CalculateKernel(sf::Vector2f positionX, sf::Vector2f positionY) {
    float r = std::sqrt(std::pow(positionX.x - positionY.x, 2)
        + std::pow(positionX.y - positionY.y, 2)); // distance of the particles
    float d = r / h;
    float t1 = std::max(1 - d, 0.0f);
    float t2 = std::max(2 - d, 0.0f);
    float w = alpha * (t2 * t2 * t2 - 4 * t1 * t1 * t1);
    
    return w;
}

sf::Color CalculateColor(float velocityMagnitude) {
    //sf::Color color;
    //float normalizedVelocity = velocityMagnitude / maxDefinedVelocity;

    
    //if (normalizedVelocity > 1.0f) {
    //    normalizedVelocity = 1.0f;
    //}
    //if (normalizedVelocity < 0.01f) {
    //    return sf::Color(0, 0, 255); 
    //}

    //if (normalizedVelocity < 0.5f) {
    
    //    float factor = normalizedVelocity / 0.5f;  
    //    color.r = static_cast<sf::Uint8>(0);  
    //    color.g = static_cast<sf::Uint8>(255 * factor); 
    //    color.b = 255; 
    //}
    //else if (normalizedVelocity < 0.7f) {
    
    //    float factor = (normalizedVelocity - 0.5f) / 0.2f;  
    //    color.r = static_cast<sf::Uint8>(0);  
    //    color.g = 255; 
    //    color.b = static_cast<sf::Uint8>(255 * (1.0f - factor));  
    //}
    //else if (normalizedVelocity < 0.9f) {
    
    //    float factor = (normalizedVelocity - 0.7f) / 0.2f;  
    //    color.r = static_cast<sf::Uint8>(255 * factor); 
    //    color.g = 255;  
    //    color.b = static_cast<sf::Uint8>(0);  
    //}
    //else {
    
    //    float factor = (normalizedVelocity - 0.9f) / 0.1f; 
    //    color.r = 255;  
    //    color.g = static_cast<sf::Uint8>(255 * (1.0f - factor));  
    //    color.b = static_cast<sf::Uint8>(0);  
    //}
    sf::Color color;
    float normalizedVelocity = velocityMagnitude / maxDefinedVelocity;

    
    if (normalizedVelocity > 1.0f) {
        normalizedVelocity = 1.0f;
    }
    if (normalizedVelocity < 0.1f) {
        return sf::Color::Blue;  
    }

   
    float factor = pow(normalizedVelocity, 1.5f);  

    color.r = static_cast<sf::Uint8>(0);  
    color.g = static_cast<sf::Uint8>(255 * factor);  
    color.b = 255;  

    return color;
    
}


float Calculate_Boundary_Pressure(size_t i, ParticleContainer& pc) {
    
    float left = 0.0f;
    sf::Vector2f right = sf::Vector2f(0.0f, 0.0f);
    float bottom = 0.0f;

    float beta = 0.0f;
    float gam = 0.0f;

    float a = 0.0f;
    float c = 0.0f;
    sf::Vector2f nom = sf::Vector2f(0.0f, 0.0f);
    sf::Vector2f b = sf::Vector2f(0.0f, 0.0f);
    float denom = 0.0f;
    float a_nom = 0.0f;
    float a_denom = 0.0f;
    float M[2][2] = { {1.0f, 0.0f}, {0.0f, 1.0f} };




    std::vector<size_t> neighbors = pc.fluidParticles[i].getFluidNeighbors();
    // pressure extrapolation
    if (mls == false) {
        for (auto& j : neighbors) {
            if (pc.fluidParticles[j].boundary_ == false) {
                float w = CalculateKernel(pc.fluidParticles[j].getPosition(), pc.fluidParticles[i].getPosition());

                left = left + pc.fluidParticles[j].getPressure() * w;
                right = right + pc.fluidParticles[j].getDensity() * (pc.fluidParticles[i].getPosition() - pc.fluidParticles[j].getPosition()) * w;
                bottom = bottom + w;
            }
        }
        float res = left + Helper::ScalarProduct(sf::Vector2f(0.0f, gravity), right) / bottom;
        return std::max(res, 0.0f);

    }

    // mls
    else {
        // calculate db
        for (auto& j : neighbors) {
            if (pc.fluidParticles[j].boundary_ == false) {
                float w = CalculateKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition());
                float v = mass / p0;
                nom = nom + v * pc.fluidParticles[j].getPosition() * w;
                denom = denom + v * w;
            }
        } 

        // transform positions
        sf::Vector2f db = nom / denom;
        
        pc.fluidParticles[i].transformedPosition = pc.fluidParticles[i].getPosition() - db;
        for (auto& j : neighbors) {
            if (pc.fluidParticles[j].boundary_ == false) {
                pc.fluidParticles[j].transformedPosition = pc.fluidParticles[j].getPosition() - db;
            }
        }
        // calculate alpha
        for (auto& j : neighbors) {
            if (pc.fluidParticles[j].boundary_ == false) {
                float w = CalculateKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition());
                float v = mass / p0;
                a_nom = a_nom + v * pc.fluidParticles[j].getPressure() * w;
                a_denom = a_denom + v * w;
            }
        }
        a = a_nom / a_denom;
        // calculate b
        for (auto& j : neighbors) {
            if (pc.fluidParticles[j].boundary_ == false) {
                float w = CalculateKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition());
                float v = mass / p0;
                b = b + pc.fluidParticles[j].transformedPosition * pc.fluidParticles[j].getPressure() * v * w;
            }
        }
        // calculate matrix m 
        for (auto& j : neighbors) {
            if (pc.fluidParticles[j].boundary_ == false) {
                float w = CalculateKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition());
                float v = mass / p0;
                sf::Vector2f transf_pos = pc.fluidParticles[j].transformedPosition;
                M[0][0] = M[0][0] + transf_pos.x * transf_pos.x * v * w;
                M[0][1] = M[0][1] + transf_pos.x * transf_pos.y * v * w;
                M[1][0] = M[1][0] + transf_pos.y * transf_pos.x * v * w;
                M[1][1] = M[1][1] + transf_pos.y * transf_pos.y * v * w;
            }
        }

        // calculate inverse 

        // calculate determinante
        float det = M[0][0] * M[1][1] - M[0][1] * M[1][0];

        float invDet = 1.0f / det; // Calculate the inverse of the determinant

        float M00 = M[0][0];
        float M01 = M[0][1];
        float M10 = M[1][0];
        float M11 = M[1][1];

        M[0][0] = M11 * invDet;
        M[0][1] = -M01 * invDet;
        M[1][0] = -M10 * invDet;
        M[1][1] = M00 * invDet;

        sf::Vector2f product = sf::Vector2f(0.0f, 0.0f);
        product.x = M[0][0] * b.x + M[0][1] * b.y;
        product.y = M[1][0] * b.x + M[1][1] * b.y;

        float press = a + product.x * pc.fluidParticles[i].transformedPosition.x + product.y * pc.fluidParticles[i].transformedPosition.y;

        return std::max(press, 0.0f);
    }
    return 0.0f;
}


sf::Vector2f CalculatePressureAcceleration(size_t i, ParticleContainer& pc) {

    sf::Vector2f pressureacc = sf::Vector2f(0.0f, 0.0f);
    sf::Vector2f pressureaccfluid = sf::Vector2f(0.0f, 0.0f);
    sf::Vector2f pressureaccboundary = sf::Vector2f(0.0f, 0.0f);
    std::vector<size_t> neighbors = pc.fluidParticles[i].getFluidNeighbors();

    if (mirroring) {


        for (auto& j : neighbors) {

            if (pc.fluidParticles[j].boundary_)
            {
                
                pressureaccboundary = pressureaccboundary + pc.fluidParticles[j].mass_ *
                    ((pc.fluidParticles[i].getPressure() / (p0 * p0)) +
                        (pc.fluidParticles[i].getPressure() / (p0 * p0))) * Helper::gradientKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition());
            }
            else
            {
                pressureaccfluid = pressureaccfluid + pc.fluidParticles[j].mass_ *
                    ((pc.fluidParticles[i].getPressure() / (p0 * p0)) +
                        (pc.fluidParticles[j].getPressure() / (p0 * p0))) * Helper::gradientKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition());

            }

        }


        pressureacc = -1.0f * pressureaccfluid + -1.0f * boundary_gamma * pressureaccboundary;

        // std::cout << "at particle first loop " << i << "x:   " << pressureacc.x << "y:   " << pressureacc.y << "   aff:    " << pc.fluidParticles[i].aff_ << "      sf   " << pc.fluidParticles[i].sf_ << " iter  " << count << "  pressure  " << pc.fluidParticles[i].getPressure() << "neighbor count " << neighbors.size() << std::endl;
        pc.fluidParticles[i].boundaryPressureAcceleration = -1.0f * boundary_gamma * pressureaccboundary;
        pc.fluidParticles[i].particlePressureAcceleration = -1.0f * pressureaccfluid;
        
        return pressureacc;

    }
    else
    {
        for (auto& j : neighbors) {
            if (pc.fluidParticles[j].boundary_)
            {
               
                pc.fluidParticles[j].setPressure(Calculate_Boundary_Pressure(j, pc));
                
                pressureaccboundary = pressureaccboundary + pc.fluidParticles[j].mass_ *
                    ((pc.fluidParticles[i].getPressure() / (p0 * p0)) +
                        (pc.fluidParticles[j].getPressure() / (p0 * p0))) * Helper::gradientKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition());
            }
            else
            {
                pressureaccfluid = pressureaccfluid + pc.fluidParticles[j].mass_ *
                    ((pc.fluidParticles[i].getPressure() / (p0 * p0)) +
                        (pc.fluidParticles[j].getPressure() / (p0 * p0))) * Helper::gradientKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition());

            }

        }


        pressureacc = -1.0f * pressureaccfluid + -1.0f * pressureaccboundary;

        // std::cout << "at particle first loop " << i << "x:   " << pressureacc.x << "y:   " << pressureacc.y << "   aff:    " << pc.fluidParticles[i].aff_ << "      sf   " << pc.fluidParticles[i].sf_ << " iter  " << count << "  pressure  " << pc.fluidParticles[i].getPressure() << "neighbor count " << neighbors.size() << std::endl;
        pc.fluidParticles[i].boundaryPressureAcceleration = -1.0f * boundary_gamma * pressureaccboundary;
        pc.fluidParticles[i].particlePressureAcceleration = -1.0f * pressureaccfluid;
        return pressureacc;
    }
}


float CalculateDensity(size_t i, ParticleContainer& pc)
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
        // float w = alpha * (t2 * t2 * t2 - 4 * t1 * t1 * t1);


        float w = CalculateKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition());
        temp = temp + pc.fluidParticles[j].mass_ * w;

    }
    if (neighbors.size() == 0) {
        return p0;
    }
    if (temp == 0.0f) {
        return p0;
    }

    return temp;
}

sf::Vector2f CalculateAcceleration(size_t i, ParticleContainer& pc)
{
    if (pc.fluidParticles[i].boundary_)
    {
        return sf::Vector2f(0, 0);
    }
    std::vector<size_t> neighbors = pc.fluidParticles[i].getFluidNeighbors();

    sf::Vector2f pressureAcceleration;
    sf::Vector2f viscosityAcceleration;
    sf::Vector2f kernel;

    sf::Vector2f surfacetension;
    sf::Vector2f acceleration;

    for (auto& j : neighbors)
    {

        viscosityAcceleration = viscosityAcceleration + (float)((pc.fluidParticles[j].mass_ / pc.fluidParticles[j].getDensity()) *
            ((Helper::ScalarProduct(pc.fluidParticles[i].getVelocity() - pc.fluidParticles[j].getVelocity(), pc.fluidParticles[i].getPosition() - pc.fluidParticles[j].getPosition()))
                / (Helper::ScalarProduct(pc.fluidParticles[i].getPosition() - pc.fluidParticles[j].getPosition(), pc.fluidParticles[i].getPosition() - pc.fluidParticles[j].getPosition()) + 0.01 * h * h))) * Helper::gradientKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition());

        surfacetension = surfacetension + mass * (pc.fluidParticles[i].getPosition() - pc.fluidParticles[j].getPosition()) * CalculateKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition());


    }


    surfacetension = surfacetension * -stension / mass;


    viscosityAcceleration = 2 * v * viscosityAcceleration;

    acceleration = viscosityAcceleration + sf::Vector2f(0, gravity * pc.fluidParticles[i].mass_) + surfacetension;


    return acceleration;
}


void  sourceterm(size_t i, ParticleContainer& pc) {


    float sf = 0;
    float f = 0.0f;
    float b = 0.0f;
    std::vector<size_t> neighbors = pc.fluidParticles[i].getFluidNeighbors();
    for (auto& j : neighbors)
    {

        if (pc.fluidParticles[j].boundary_ == false) {
            f = f + pc.fluidParticles[j].mass_ * Helper::ScalarProduct(pc.fluidParticles[i].getPredictedVelocity() - pc.fluidParticles[j].getPredictedVelocity(), Helper::gradientKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition()));
        }
        if (pc.fluidParticles[j].boundary_) {
            b = b + pc.fluidParticles[j].mass_ * Helper::ScalarProduct(pc.fluidParticles[i].getPredictedVelocity(), Helper::gradientKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition()));
        }

    }

    sf = p0 - pc.fluidParticles[i].getDensity() - t * f - t * b;
    //sf = p0 - pc.fluidParticles[i].getDensity();

    pc.fluidParticles[i].sf_ = sf;

}


void diagonalElement(size_t i, ParticleContainer& pc) {
    std::vector<size_t> neighbors = pc.fluidParticles[i].getFluidNeighbors();
    std::vector<size_t> neighborsneighbors;
    float fluid = 0.0f;
    float border = 0.0f;
    sf::Vector2f fluidresult = sf::Vector2f(0.0f, 0.0f);
    sf::Vector2f boundaryresult = sf::Vector2f(0.0f, 0.0f);
    float acc = 0.0f;
    float acctwo = 0.0f;
    float accthree = 0.0f;
    float facc = 0.0f;
    sf::Vector2f diagElement = sf::Vector2f(0.0f, 0.0f);
    float aff = 0.0f;
    sf::Vector2f result = sf::Vector2f(0.0f, 0.0f);

    for (auto& j : neighbors) {
        if (pc.fluidParticles[j].boundary_ == false) {
            neighborsneighbors = pc.fluidParticles[i].getFluidNeighbors();
            for (auto& jtwo : neighborsneighbors) {
                if (pc.fluidParticles[jtwo].boundary_ == false) {
                    fluidresult = fluidresult + (pc.fluidParticles[jtwo].mass_ / (p0 * p0)) * Helper::gradientKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[jtwo].getPosition());
                }
                else {
                    boundaryresult = boundaryresult + (pc.fluidParticles[jtwo].mass_ / (p0 * p0)) * Helper::gradientKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[jtwo].getPosition());

                }
            }
            acc = acc + pc.fluidParticles[j].mass_ * Helper::ScalarProduct((-1.0f * fluidresult + -2.0f * boundary_gamma * boundaryresult), Helper::gradientKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition()));
            boundaryresult = sf::Vector2f(0.0f, 0.0f);
            fluidresult = sf::Vector2f(0.0f, 0.0f);
        }
    }



    for (auto& j : neighbors) {
        if (pc.fluidParticles[j].boundary_ == false) {
            acctwo = acctwo + pc.fluidParticles[j].mass_ * Helper::ScalarProduct((pc.fluidParticles[j].mass_ / (p0 * p0)) * Helper::gradientKernel(pc.fluidParticles[j].getPosition(), pc.fluidParticles[i].getPosition())
                , Helper::gradientKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition()));
        }
    }
    boundaryresult = sf::Vector2f(0.0f, 0.0f);
    fluidresult = sf::Vector2f(0.0f, 0.0f);


    for (auto& j : neighbors) {
        if (pc.fluidParticles[j].boundary_) {
            neighborsneighbors = pc.fluidParticles[i].getFluidNeighbors();
            for (auto& jtwo : neighborsneighbors) {
                if (pc.fluidParticles[jtwo].boundary_ == false) {
                    fluidresult = fluidresult + (pc.fluidParticles[jtwo].mass_ / (p0 * p0)) * Helper::gradientKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[jtwo].getPosition());
                }
                else {
                    boundaryresult = boundaryresult + (pc.fluidParticles[jtwo].mass_ / (p0 * p0)) * Helper::gradientKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[jtwo].getPosition());

                }
            }
            accthree = accthree + pc.fluidParticles[j].mass_ * Helper::ScalarProduct((-1.0f * fluidresult + -2.0f * boundary_gamma * boundaryresult), Helper::gradientKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition()));
            boundaryresult = sf::Vector2f(0.0f, 0.0f);
            fluidresult = sf::Vector2f(0.0f, 0.0f);
        }
    }


    aff = t * t * acc + t * t * acctwo + accthree * t * t;

    pc.fluidParticles[i].aff_ = aff;

}


void solver(ParticleContainer& pc) {
    // for all pressure = 0
    size_t l = 0;
    sf::Vector2f pressureacc = sf::Vector2f(0.0f, 0.0f);
    sf::Vector2f pressureaccfluid = sf::Vector2f(0.0f, 0.0f);
    sf::Vector2f pressureaccboundary = sf::Vector2f(0.0f, 0.0f);

    float predrhoerror = 10.0f;
    float rhoerror = 0.0f;

    float Ap = 0.0f;
    int count = 0;
    int num = 0;
    float pressure = 0.0f;


    while (std::abs(predrhoerror) > 0.1f && count < 500) {
        predrhoerror = 0.0f;
        rhoerror = 0.0f;
        for (int i = 0; i < pc.fluidParticles.size(); i++) {

            if (pc.fluidParticles[i].boundary_ == false) {

                pc.fluidParticles[i].setPAcceleration(CalculatePressureAcceleration(i, pc));
 
            }
        }
        for (int i = 0; i < pc.fluidParticles.size(); i++) {
            if (pc.fluidParticles[i].boundary_ == false) {
                std::vector<size_t> neighbors = pc.fluidParticles[i].getFluidNeighbors();
                float Apfluid = 0.0f;
                float Apboundary = 0.0f;
                for (auto& j : neighbors) {
                    if (pc.fluidParticles[j].boundary_ == false)
                    {
                        Apfluid = Apfluid + pc.fluidParticles[j].mass_ * Helper::ScalarProduct((pc.fluidParticles[i].getPAcceleration() - pc.fluidParticles[j].getPAcceleration()), Helper::gradientKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition()));

                    }
                    else
                    {
                        Apboundary = Apboundary + pc.fluidParticles[j].mass_ * Helper::ScalarProduct(pc.fluidParticles[i].getPAcceleration(), Helper::gradientKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition()));
                    }
                }

                Ap = t * t * Apfluid + t * t * Apboundary;
                pc.fluidParticles[i].ap_ = Ap;

                // std::cout << "at particle " << i << "   aff   " << pc.fluidParticles[i].aff_ << "    sf:   " << pc.fluidParticles[i].sf_ << "    pressure   " << pc.fluidParticles[i].getPressure() << "   ap " << pc.fluidParticles[i].ap_ << "  iter " << count << std::endl;

                Ap = 0.0f;
                Apfluid = 0.0f;
                Apboundary = 0.0f;


                if (pc.fluidParticles[i].aff_ != 0.0f) {
                    pressure = std::max(pc.fluidParticles[i].getPressure() + omega * ((pc.fluidParticles[i].sf_ - pc.fluidParticles[i].ap_) / pc.fluidParticles[i].aff_), 0.0f);

                    pc.fluidParticles[i].setPressure(pressure);
                }
                pressure = 0.0f;

                float predDensityError = pc.fluidParticles[i].ap_ - pc.fluidParticles[i].sf_;
                if (pc.fluidParticles[i].ap_ >= pc.fluidParticles[i].sf_) {
                    predrhoerror = predrhoerror + ((predDensityError));
                    rhoerror = rhoerror + std::abs(pc.fluidParticles[i].getDensity() - p0);
                    num++;
                }


            }
        }

        predrhoerror = ((predrhoerror / (float)num) / p0) * 100.0f;
        rhoerror = ((rhoerror / float(num)) / p0) * 100.0f;

        // std::cout << predrhoerror << std::endl;
        file10 << predrhoerror << ", " << count << ", " << omega << std::endl;
        file11 << rhoerror << ", " << count << std::endl;
        // std::cout << num << "iteration" << count << std::endl;
        num = 0; 
        
        
        count++;
    }


    sf::Vector2f acc = sf::Vector2f(0.0f, 0.0f);

    for (int i = 0; i < pc.fluidParticles.size(); i++) {

        if (pc.fluidParticles[i].boundary_ == false) {

            pc.fluidParticles[i].setPAcceleration(CalculatePressureAcceleration(i, pc));
        }
    }
    // std::cout << "dichtefehler       " << predrhoerror << "  iteration  " << count << std::endl;
    iter = count;
    predDensityErrors.push_back(predrhoerror);
    DensityErrors.push_back(rhoerror);
}


void pressure_acc(ParticleContainer& pc, int i) {
    sf::Vector2f pressureacc = sf::Vector2f(0.0f, 0.0f);
    sf::Vector2f pressureaccfluid = sf::Vector2f(0.0f, 0.0f);
    sf::Vector2f pressureaccboundary = sf::Vector2f(0.0f, 0.0f);
    sf::Vector2f acc = sf::Vector2f(0.0f, 0.0f);

    std::vector<size_t> neighbors = pc.fluidParticles[i].getFluidNeighbors();
    for (auto& j : neighbors) {

        if (pc.fluidParticles[j].boundary_)
        {
            pressureaccboundary = pressureaccboundary + pc.fluidParticles[j].mass_ *
                ((pc.fluidParticles[i].getPressure() / (p0 * p0)) +
                    (pc.fluidParticles[i].getPressure() / (p0 * p0))) * Helper::gradientKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition());
        }
        else
        {
            pressureaccfluid = pressureaccfluid + pc.fluidParticles[j].mass_ *
                ((pc.fluidParticles[i].getPressure() / (p0 * p0)) +
                    (pc.fluidParticles[j].getPressure() / (p0 * p0))) * Helper::gradientKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition());

        }
        if (i == 10) {
            acc = acc + Helper::gradientKernel(pc.fluidParticles[i].getPosition(), pc.fluidParticles[j].getPosition());
        }
    }

    pressureacc = -1.0f * pressureaccfluid + -1.0f * boundary_gamma * pressureaccboundary;

    // std::cout << "at particle first loop " << i << "x:   " << pressureacc.x << "y:   " << pressureacc.y << "   aff:    " << pc.fluidParticles[i].aff_ << "      sf   " << pc.fluidParticles[i].sf_ << " iter  " << count << "  pressure  " << pc.fluidParticles[i].getPressure() << "neighbor count " << neighbors.size() << std::endl;
    pc.fluidParticles[i].boundaryPressureAcceleration = -1.0f * boundary_gamma * pressureaccboundary;
    pc.fluidParticles[i].particlePressureAcceleration = -1.0f * pressureaccfluid;


    // pc.fluidParticles[i].setPAcceleration(pressureacc);



}



void ISPH::Update(ParticleContainer& pc, UniformGrid& uniform_grid) {

    // neighborhoodsearch
    auto start_ns = std::chrono::high_resolution_clock::now();
#pragma omp parallel for 
    for (int i = 0; i < pc.fluidParticles.size(); i++)
    {
        if (mirroring == false) {
            pc.fluidParticles[i].setNeighborIndices(uniform_grid.getNeighbors(i, pc.fluidParticles[i].getPosition(), pc));
        }
        else {
            if (pc.fluidParticles[i].boundary_ == false) {
                pc.fluidParticles[i].setNeighborIndices(uniform_grid.getNeighbors(i, pc.fluidParticles[i].getPosition(), pc));
            }
        }
    }
    auto end_ns = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_ns = end_ns - start_ns;
    

    // set pressure and density loop

    auto start_pd = std::chrono::high_resolution_clock::now();
#pragma omp parallel for 
    for (int i = 0; i < pc.fluidParticles.size(); i++)
    {
        if (pc.fluidParticles[i].boundary_ == false)
        {
            pc.fluidParticles[i].setDensity(CalculateDensity(i, pc));
            
            
        }
        pc.fluidParticles[i].setPressure(0.0f);
        pc.fluidParticles[i].setPAcceleration(sf::Vector2f(0.0f, 0.0f));
    }
    auto end_pd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_pd = end_pd - start_pd;
    // file8 << "Density: " << elapsed_pd.count() << " seconds" << std::endl;

    // set pressure and density loop
    auto start_pv = std::chrono::high_resolution_clock::now();
#pragma omp parallel for 
    for (int i = 0; i < pc.fluidParticles.size(); i++)
    {
        if (pc.fluidParticles[i].boundary_ == false)
        {
            pc.fluidParticles[i].setNonPAcceleration(CalculateAcceleration(i, pc));
            pc.fluidParticles[i].setPredictedVelocity();
        }
    }
    auto end_pv = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_pv = end_pv - start_pv;
    // file8 << "Predicted Velocity: " << elapsed_pv.count() << " seconds" << std::endl;

    // sourceterm and diagonal element


    auto start_sd = std::chrono::high_resolution_clock::now();
#pragma omp parallel for 
    for (int i = 0; i < pc.fluidParticles.size(); i++)
    {
        if (pc.fluidParticles[i].boundary_ == false)
        {
            sourceterm(i, pc);
            diagonalElement(i, pc);
        }
    }
    auto end_sd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_sd = end_sd - start_sd;
    // file8 << "Source term and diagonal Element: " << elapsed_sd.count() << " seconds" << std::endl;



    auto start_solver = std::chrono::high_resolution_clock::now();
    solver(pc);
    auto end_solver = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_solver = end_solver - start_solver;
    // file8 << "Solver: " << elapsed_solver.count() << " seconds" << std::endl;
    


    auto start_ng = std::chrono::high_resolution_clock::now();
    uniform_grid.GridClear();
    auto end_ng = std::chrono::high_resolution_clock::now();
    elapsed_ns = elapsed_ns + start_ng - end_ng;
    // file8 << "Neighborhood: " << elapsed_ns.count() << " seconds" << std::endl;




    int size = pc.fluidParticles.size() - 1;

   

    maxVel1 = 0.0f;

    auto start_loop = std::chrono::high_resolution_clock::now();  

    std::chrono::duration<double> total_without_add_particle(0);  


    for (size_t i = 0; i < pc.fluidParticles.size(); i++)
    {
        auto start_particle = std::chrono::high_resolution_clock::now();
        if (pc.fluidParticles[i].boundary_ == false)
        {
            
            pc.fluidParticles[i].setColor(sf::Color(sf::Color::Blue));
            pc.fluidParticles[i].setAcceleration();
            pc.fluidParticles[i].setVelocity();
            pc.fluidParticles[i].setPosition();
        }
        else
        {
            pc.fluidParticles[i].setColor(sf::Color(75, 75, 75));
        }
        
        if (MagnitudeVector(pc.fluidParticles[i].getVelocity()) > maxVel1)
        {
            maxVel1 = MagnitudeVector(pc.fluidParticles[i].getVelocity());
        }

        if (pc.fluidParticles[i].boundary_ == false) {
            float velocityMagnitude = MagnitudeVector(pc.fluidParticles[i].getVelocity());
            sf::Color particleColor = CalculateColor(velocityMagnitude);
            pc.fluidParticles[i].setColor(particleColor);

        }
        auto end_particle_before_add = std::chrono::high_resolution_clock::now();


        uniform_grid.AddParticle(pc.fluidParticles[i].getPosition(), i);

        auto end_particle_after_add = std::chrono::high_resolution_clock::now();
        total_without_add_particle += (end_particle_before_add - start_particle);

    }
    // file8 << "Position Update: " << total_without_add_particle.count() << std::endl;

    //float cfl = 0.4f;
    //float deltat = 0.0f;
    //float maxtime = 0.02f;

    //if (maxVel1 > 0) {
    //    deltat = (h * cfl) / maxVel1;
    //}
    //else {
    //    deltat = maxtime;
    //}


    //t = std::min(maxtime, deltat);

    //t = std::max(t, 0.0001f);

}
