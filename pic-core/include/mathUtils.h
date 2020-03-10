#include <array>
#include <cmath>
#include <random>
#include <chrono>
#include <vector>
#include <iostream>
#include <iomanip>

#define DEBUG() std::cout<<"BREAK POINT: LINE "<<__LINE__<<" IN "<<__FILE__<<std::endl


namespace pic
{

void printVec3(std::array<float, 3> v);
void printVec2(std::array<float, 2> v);


float getDiff(float const x0, float const x1, float const x);

float linearInterpolation(float const c0, float const c1, float const z);
float bilinearIterpolation(std::array<float, 4> attr, float y, float z);
float trilinearInterpolation(std::array<float, 8> attr, float x, float y, float z);
std::array<float, 8> getWeights(float x, float y, float z);
std::array<float, 3> cellSpacePos(float x, float y, float z, float cellLen_x, 
	float cellLen_y, float cellLen_z, int cell_i,int cell_j, int cell_k);

}