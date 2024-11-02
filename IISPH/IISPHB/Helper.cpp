#include "Helper.h"

sf::Vector2f Helper::gradientKernel(sf::Vector2f positionX, sf::Vector2f positionY) {
    float dx = positionX.x - positionY.x;
    float dy = positionX.y - positionY.y;
    float r = std::sqrt(dx * dx + dy * dy);
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

float Helper::ScalarProduct(sf::Vector2f a, sf::Vector2f b) {
    return a.x * b.x + a.y * b.y;
}


float Helper::MagniVector(sf::Vector2f vec) {
    return std::sqrt(vec.x * vec.x + vec.y * vec.y);
}