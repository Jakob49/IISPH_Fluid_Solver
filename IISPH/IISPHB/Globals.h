#pragma once
#include <vector>
#include <fstream>

// Global variables
extern int globalVariable;

extern const float h;
extern const float p0;
extern float k;
extern const float alpha;
extern float v;
extern const float gravity;
extern const float mass;
extern float t;
extern float c;
extern int iter;
extern const float gravity;
extern float boundary_gamma;
extern float omega;
extern float maxVel1;
extern float stension;
extern bool sph;



extern bool mirroring;
extern bool mls;
extern float maxDefinedVelocity;


extern std::vector<float> predDensityErrors;
extern std::vector<float> DensityErrors;
extern std::ofstream file8;
extern std::ofstream file9;
extern std::ofstream file10;
extern std::ofstream file11;