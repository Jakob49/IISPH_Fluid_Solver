#include "UniformGrid.h"
#include "ParticleContainer.h"
#include "Globals.h"

#include <iostream>

UniformGrid::UniformGrid()
{
	size_t numCellsX = std::ceil(1000.0 / 2 * h);
	size_t numCellsY = std::ceil(800.0 / 2 * h);
	cells_.resize(numCellsX);
	for (size_t i = 0; i < numCellsX; i++) {
		cells_[i].resize(numCellsY);
	}

}

void UniformGrid::AddParticle(sf::Vector2f position, size_t index)
{
	int cellx = static_cast<int>(position.x / 2 * h);
	int celly = static_cast<int>(position.y / 2 * h);

	if (cellx >= 0 && cellx < cells_.size() && celly >= 0 && celly < cells_[cellx].size()) {
		cells_[cellx][celly].push_back(index);
	}
}

void UniformGrid::GridClear()
{
	for (int i = 0; i < cells_.size(); ++i) {
		for (int j = 0; j < cells_[i].size(); ++j) {
			cells_[i][j].clear();
		}
	}
}

std::vector<size_t> UniformGrid::getNeighbors(size_t index, sf::Vector2f position, ParticleContainer& pc)
{

	std::vector<size_t> potentialNeighbors;
	std::vector<size_t> neighborIndices;

	const float epsilon = 0.0000001f;

	int cellx = static_cast<int>(position.x / 2 * h);
	int celly = static_cast<int>(position.y / 2 * h);
	for (int i = -1; i <= 1; i++)
	{
		for (int j = -1; j <= 1; j++)
		{
			int neighborCellX = cellx + i;
			int neighborCellY = celly + j;

			// Bounds checking
			if (neighborCellX >= 0 && neighborCellX < cells_.size() &&
				neighborCellY >= 0 && neighborCellY < cells_[neighborCellX].size()) {

				for (auto& particleIndex : cells_[neighborCellX][neighborCellY]) {
					potentialNeighbors.push_back(particleIndex);

				}
			}
		}
	}



	for (size_t i : potentialNeighbors) {
		float distance = CalulateDistance(pc.fluidParticles[index].getPosition(), pc.fluidParticles[i].getPosition());
		if (distance <= 2 * h + epsilon) {
			neighborIndices.push_back(i);
		}
	}
	return neighborIndices;

}
float UniformGrid::CalulateDistance(sf::Vector2f a, sf::Vector2f b)
{
	return std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2));
}

int UniformGrid::GiveClosestParticle(sf::Vector2f click, ParticleContainer pc) {
	float minDistance = std::numeric_limits<float>::max(); // Initialize to a very large value
	int closestParticleIndex = -1;

	int cellx = static_cast<int>(click.x / 2 * h);
	int celly = static_cast<int>(click.y / 2 * h);
	for (int i = -1; i <= 1; i++) {
		for (int j = -1; j <= 1; j++) {
			int neighborCellX = cellx + i;
			int neighborCellY = celly + j;

			// Bounds checking
			if (neighborCellX >= 0 && neighborCellX < cells_.size() &&
				neighborCellY >= 0 && neighborCellY < cells_[neighborCellX].size()) {

				for (auto& particleIndex : cells_[neighborCellX][neighborCellY]) {
					float distance = CalulateDistance(click, pc.fluidParticles[particleIndex].getPosition());
					if (distance < minDistance) {
						minDistance = distance;
						closestParticleIndex = particleIndex;
					}
				}
			}
		}
	}

	return closestParticleIndex;
}
