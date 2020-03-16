#include <array>
#include <cmath>
#include <random>
#include <chrono>
#include <vector>
#include <iostream>
#include <iomanip>

#define DEBUG() std::cout<<"BREAK POINT: LINE "<<__LINE__<<" IN "<<__FILE__<<std::endl

using vecd3 = std::array<double,3>;

template<typename T>
using vec3 = std::array<T,3>;

namespace pic
{

void printVec3(std::array<float, 3> v);
void printVec2(std::array<float, 2> v);

float getDiff(float const x0, float const x1, float const x);

float linearInterpolation(float const c0, float const c1, float const z);
float bilinearIterpolation(std::array<float, 4> attr, float y, float z);
float trilinearInterpolation(std::array<float, 8> attr, float x, float y, float z);
std::array<float, 8> getWeights(float x, float y, float z);

double dot(vecd3 const& v1, vecd3 const& v2);
vecd3 scale(vecd3 const& v1, double s);
vecd3 const project(vecd3 const& v1, vecd3 const& v2);
double length(vecd3 const& v1);
vecd3 const add(vecd3 const& v1, vecd3 const& v2);
vecd3 const minus(vecd3 const& v1, vecd3 const& v2);

int intersectRaySphere(vecd3 const& pointPos, vecd3 const& pointVel, 
	vecd3 const& sphereCenter, double sphereRad, double &t, vecd3 &q);
}