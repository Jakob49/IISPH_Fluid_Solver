#pragma once
#pragma once
#include <SFML/Graphics.hpp>
#include <vector>


class Particle
{
private:
    sf::CircleShape shape_;
    sf::Vector2f velocity_;
    sf::Vector2f position_;
    sf::Vector2f acceleration_;
    sf::Vector2f npacceleration_;
    sf::Vector2f pacceleration_;
    sf::Vector2f predictedVelocity_;


    std::vector<size_t> fluidNeighbors_;
    // std::vector<size_t>& boundaryNeighbors_;
    std::vector<sf::Color> interpolateColors;

    float density_;
    float densityError_;

    float pressure_;
    float kernel_;

    float radius_;







public:
    Particle(sf::Vector2f position, bool bound);
    float mass_;

    sf::Vector2f particlePressureAcceleration;
    sf::Vector2f boundaryPressureAcceleration;
    sf::Vector2f transformedPosition;

    float aff_;
    float ap_;
    float sf_;

    sf::Color color_;
    bool boundary_;

    void Update();
    void UpdateShape();
    void setBoundaryVelocity(sf::Vector2f vel);
    void setBoundaryPosition(sf::Vector2f pos);

    void setNeighborIndices(const std::vector<size_t>& neighborIndices);



    sf::CircleShape getShape() const;

    sf::Vector2f getPosition() const { return position_; }
    sf::Vector2f getVelocity() const { return velocity_; }
    sf::Vector2f getNonPAcceleration() const { return npacceleration_; }
    sf::Vector2f getPAcceleration() const { return pacceleration_; }
    sf::Vector2f getPredictedVelocity() const { return predictedVelocity_; }


    float getDensity() const { return density_; }
    float getPressure() const { return pressure_; }
    float getDensityError() const { return densityError_; }

    sf::Color getColor() const { return color_; }
    std::vector<size_t> getFluidNeighbors();


    void setColor(sf::Color color);
    void setDensity(float temp);
    void setPressure(float pressure);
    void setPredictedVelocity();
    void setVelocity();
    void setPosition();
    void setAcceleration();
    void setNonPAcceleration(sf::Vector2f acc);
    void setPAcceleration(sf::Vector2f acc);



};
