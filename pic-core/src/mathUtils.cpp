#include "mathUtils.h"

namespace pic
{

void printVec3(std::array<float, 3> v)
{
	std::cout<<"x: "<<std::setprecision(22)<<v[0]
		   <<", y: "<<std::setprecision(22)<<v[1]
		   <<", z: "<<std::setprecision(22)<<v[2]<<std::endl;
}

void printVec2(std::array<float, 2> v)
{
	std::cout<<"y: "<<v[0]<<", z: "<<v[1]<<std::endl;
}

float getDiff(float const x0, float const x1, float const x)
{
	printVec3({x0, x1, x});
	if ((x1 - x0) == 0){
		DEBUG();
	}

	return (x - x0)/(x1 - x0);
}

float linearInterpolation(float const c0, float const c1, float const z)
{
	return c0 * (1.f - z) + c1 * z;
}

float bilinearIterpolation(std::array<float, 4> attr, float y, float z)
{
	auto [c00, c01, c10, c11] = attr;

	float const c0 = linearInterpolation(c00, c10, y);
	float const c1 = linearInterpolation(c01, c11, y);

	return linearInterpolation(c0, c1, z);
}

float trilinearInterpolation(std::array<float, 8> attr, float x, float y, float z)
{
	auto [c000, c100, c001, c101, c010, c110, c011, c111] = attr;

	float const c00 = linearInterpolation(c000, c100, x);
	float const c01 = linearInterpolation(c001, c101, x);
	float const c10 = linearInterpolation(c010, c110, x);
	float const c11 = linearInterpolation(c011, c111, x);

	return bilinearIterpolation({c00, c01, c10, c11}, y, z);
}

std::array<float, 8> getWeights(float x, float y, float z)
{
	float w000 = (1 - z) * (1 - y) *  (1 - x);
	float w100 = z * (1 - y) *  (1 - x);
	float w001 = (1 - z) *  (1 - y) * x;
	float w101 = z *  (1 - y) * x;
	float w010 = (1 - z) * y *  (1 - x);
	float w110 = z * y *  (1 - x);
	float w011 = (1 - z) *  y * x;
	float w111 = z *  y * x;
	
	return {w000, w100, w001, w101, w010, w110, w011, w111};
}

double dotProd(auto const& v1, auto const& v2)
{
	return {v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]}; 
}

std::array<float,3> scale(auto const& v1, double s)
{ 
	return {v1[0] * s + v1[1] * s + v1[2] * s}; 
}

std::array<float,3> const project(auto const& v1, auto const& v2)
{ 
	return scale(v1, dotProd(v1, v2) / dotProd(v1, v1));
}

double length(auto const& v1)
{
	return sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
}

}