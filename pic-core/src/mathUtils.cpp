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
	float w000 = (1 - z) * (1 - y) * (1 - x);
	float w100 = z * (1 - y) * (1 - x);
	float w001 = (1 - z) * (1 - y) * x;
	float w101 = z * (1 - y) * x;
	float w010 = (1 - z) * y * (1 - x);
	float w110 = z * y * (1 - x);
	float w011 = (1 - z) * y * x;
	float w111 = z * y * x;

	return {w000, w100, w001, w101, w010, w110, w011, w111};
}

double dot(vecd3 const& v1, vecd3 const& v2)
{
	return {v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]}; 
}

vecd3 scale(vecd3 const& v1, double s)
{ 
	return {v1[0] * s + v1[1] * s + v1[2] * s}; 
}

vecd3 const project(vecd3 const& v1, vecd3 const& v2)
{ 
	return scale(v1, dot(v1, v2) / dot(v1, v1));
}

double length(vecd3 const& v1)
{
	return sqrt(dot(v1, v1));
}

vecd3 const add(vecd3 const& v1, vecd3 const& v2)
{
	return {v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]};
}

vecd3 const minus(vecd3 const& v1, vecd3 const& v2)
{
	return {v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]};
}

int intersectRaySphere(vecd3 const& pointPos, vecd3 const& pointVel,
	vecd3 const& sphereCenter, double sphereRad, double &t, vecd3 &newPos) 
{
	vecd3 const& m =  (pointPos, sphereCenter);
	double m_len = length(m); 
	std::cout<<m_len<<std::endl;
	if(m_len > sphereRad)
	{
		std::cout<<"Outside Sphere"<<std::endl;
	}
	double b = dot(m, pointVel); 
	double c = dot(m, m) - sphereRad * sphereRad; 
	// Exit if r’s origin outside s (c > 0) and r pointing away from s (b > 0) 
	if (c > 0.0f && b > 0.0f){
		return 0;
	} 
	double discr = b*b - c; 

	// A negative discriminant corresponds to ray missing sphere 
	if (discr < 0.0f) return 0;

	// Ray now found to intersect sphere, compute smallest t value of intersection
	t = -b - sqrt(discr); 

	// If t is negative, ray started inside sphere so clamp t to zero 
	//if (t < 0.0f) t = 0.0f; 
	newPos = add(pointPos, scale(pointVel, t));
	DEBUG();

	return 1;
}

void intersectRaySphere_2()
{
	// Calculate ray start's offset from the sphere center
	vecd3 p = s - c;
	
	float rSquared = r * r;
	float p_d = dot(p, d);
	
	// The sphere is behind or surrounding the start point.
	if(p_d > 0 || dot(p, p) < rSquared)
	 return NO_COLLISION;
	
	// Flatten p into the plane passing through c perpendicular to the ray.
	// This gives the closest approach of the ray to the center.
	vecd3 a = p - p_d * d;
	
	float aSquared = dot(a, a);
	
	// Closest approach is outside the sphere.
	if(aSquared > rSquared)
	  return NO_COLLISION;
	
	// Calculate distance from plane where ray enters/exits the sphere.    
	float h = sqrt(rSquared - aSquared);
	
	// Calculate intersection point relative to sphere center.
	vecd3 i = a - h * d;
	
	vecd3 intersection = c + i;
	vecd3 normal = i/r;
	// We've taken a shortcut here to avoid a second square root.
	// Note numerical errors can make the normal have length slightly different from 1.
	// If you need higher precision, you may need to perform a conventional normalization.
	
	return (intersection, normal);
}

vecd3 raySphereIntersection(vecd3 rayStart, vecd3 rayEnd, 
	vecd3 sphereCenter, double sphereRadius)
{
	vecd3 rayStartOffset = minus(rayStart, sphereCenter);
	vecd3 rayEndOffset = minus(rayEnd, sphereCenter);

	vecd3 
	if 
	float sphereRadSquared = sphereRadius * sphereRadius;

}

}