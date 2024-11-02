#include "Globals.h"

#define _USE_MATH_DEFINES
#include<math.h>
// global constants
const float h = 1.01f;
const float p0 = 1.2f; //1.2
float k = 200000;
const float alpha = 5 / (14 * M_PI * h * h);
float v = 1.8f; //0.08
const float mass = p0 * h * h;
float t = 0.008; 
int iter = 0;
const float gravity = 9.81f;
float boundary_gamma = 0.7f; // 0.7
float omega = 0.1f;
float c = 1.0f;
float maxDefinedVelocity = 70.0f;
bool sph = false;

float stension = 0.0f;


bool mirroring = true;
bool mls = false;

float maxVel1 = 0.0f;

std::vector<float> predDensityErrors;
std::vector<float> DensityErrors;

std::ofstream file8("C:/Users/Jakob/Documents/Studium/BachelorArbeit/plots/data/timing_IISPH.txt");
std::ofstream file9("C:/Users/Jakob/Documents/Studium/BachelorArbeit/plots/data/timing_SPH.txt");
std::ofstream file10("C:/Users/Jakob/Documents/Studium/BachelorArbeit/plots/data/convergence.txt");
std::ofstream file11("C:/Users/Jakob/Documents/Studium/BachelorArbeit/plots/data/rhoerror.txt");