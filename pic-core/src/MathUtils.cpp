#include "MathUtils.h"

namespace pic
{

void printVec3(std::array<double, 3> v){
	std::cout<<"x: "<<std::setprecision(22)<<v[0]
		  <<", y: "<<std::setprecision(22)<<v[1]
		  <<", z: "<<std::setprecision(22)<<v[2]<<std::endl;
}

void printVec2(std::array<double, 2> v){
	std::cout<<"y: "<<v[0]<<", z: "<<v[1]<<std::endl;
}

double getDiff(double const x, double const x0, double const x1){
	return (x - x0)/(x1 - x0);
}
Vector3d getDiff(Vector3d const x, Vector3d const x0, Vector3d const x1){
	auto v0 = (x - x0);
	auto v1 = (x1 - x0);
	DEBUG_VAR(x);
	DEBUG_VAR(x0);
	DEBUG_VAR(x1);
	DEBUG_VAR(v0);
	DEBUG_VAR(v1);
	assert(x1[0] > x0[0]);
	assert(x1[1] > x0[1]);
	assert(x1[2] > x0[2]);
	assert(v1[0] > 0);
	assert(v1[1] > 0);
	assert(v1[2] > 0);
	assert(x[0] >= x0[0]);
	assert(x[1] >= x0[1]);
	assert(x[2] >= x0[2]);
	return {v0[0]/v1[0], v0[1]/v1[1], v0[2]/v1[2]};
}

double linearInterpolation(double const c0, double const c1, double const z){
	return c0 * (1.f - z) + c1 * z;
}

double bilinearInterpolation(tuple4d const & attr, Vector2d const pos){
	auto [c00, c01, c10, c11] = attr;
	double c0 = linearInterpolation(c00, c10, pos[0]);
	double c1 = linearInterpolation(c01, c11, pos[0]);
	return linearInterpolation(c0, c1, pos[1]);
}
double bilinearInterpolation(tuple4d const & attr, double const y, double const z){
	auto [c00, c01, c10, c11] = attr;
	double c0 = linearInterpolation(c00, c10, y);
	double c1 = linearInterpolation(c01, c11, y);
	return linearInterpolation(c0, c1, z);
}

double trilinearInterpolation(tuple8d const & attr, double const x, double const y, double const z){
	auto [c000, c100, c001, c101, c010, c110, c011, c111] = attr;
	double c00 = linearInterpolation(c000, c100, x);
	double c01 = linearInterpolation(c001, c101, x);
	double c10 = linearInterpolation(c010, c110, x);
	double c11 = linearInterpolation(c011, c111, x);
	return bilinearInterpolation({c00, c01, c10, c11}, y, z);
}
double trilinearInterpolation(tuple8d const & attr, Vector3d const & pos){
	auto [c000, c100, c001, c101, c010, c110, c011, c111] = attr;
	double c00 = linearInterpolation(c000, c100, pos[0]);
	double c01 = linearInterpolation(c001, c101, pos[0]);
	double c10 = linearInterpolation(c010, c110, pos[0]);
	double c11 = linearInterpolation(c011, c111, pos[0]);
	return bilinearInterpolation({c00, c01, c10, c11}, pos[1], pos[2]);
}

//Must be in 0 < pos < 1
inline  tuple8d __getWeights(double const x, double const y, double const z){
	DEBUG_VAR(x);
	DEBUG_VAR(y);
	DEBUG_VAR(z);
	assert((1 > x) && (x >= 0));
	assert((1 > y) && (y >= 0));
	assert((1 > z) && (z >= 0));
	double w000 = (1.f - z) * (1.f - y) * (1.f - x);
	double w100 = x * (1.f - y) * (1.f - z);
	double w001 = (1.f - x) * (1.f - y) * z;
	double w101 = x * (1.f - y) * z;
	double w010 = (1.f - x) * y * (1.f - z);
	double w110 = x * y * (1.f - z);
	double w011 = (1.f - x) * y * z;
	double w111 = x * y * z;
	//std::cout<<"\nw000: "<<w000<<"\nw100: "<<w100<<"\nw001: "<<w001<<"\nw101: "<<w101<<"\nw010: "<<w010<<"\nw110: "<<w110<<"\nw011: "<<w011<<"\nw111: "<<w111<<std::endl;
	return {w000, w100, w001, w101, w010, w110, w011, w111};
}

tuple8d getWeights(double const x, double const y, double const z){
	return __getWeights(x ,y ,z);
}
tuple8d getWeights(Vector3d const & pos){
	return __getWeights(pos[0], pos[1], pos[2]);
}
tuple8d getWeights(Vector3d const & worldSpacePos, Vector3d const & c000Pos, Vector3d const & c111Pos){
	return getWeights(getDiff(worldSpacePos, c000Pos, c111Pos));
}

inline tuple8d __getDistWeights(Vector3d const & pos){
	//DEBUG_VAR(x);
	//DEBUG_VAR(y);
	//DEBUG_VAR(z);
	//assert((1 > x) && (x >= 0));
	//assert((1 > y) && (y >= 0));
	//assert((1 > z) && (z >= 0));
	double w000 = (pos - Vector3d(0,0,0)).norm();
	double w100 = (pos - Vector3d(1,0,0)).norm();
	double w001 = (pos - Vector3d(0,0,1)).norm();
	double w101 = (pos - Vector3d(1,0,1)).norm();
	double w010 = (pos - Vector3d(0,1,0)).norm();
	double w110 = (pos - Vector3d(1,1,0)).norm();
	double w011 = (pos - Vector3d(0,1,1)).norm();
	double w111 = (pos - Vector3d(1,1,1)).norm();
	//std::cout<<"\nw000: "<<w000<<"\nw100: "<<w100<<"\nw001: "<<w001<<"\nw101: "<<w101<<"\nw010: "<<w010<<"\nw110: "<<w110<<"\nw011: "<<w011<<"\nw111: "<<w111<<std::endl;
	return {w000, w100, w001, w101, w010, w110, w011, w111};
}


tuple8d getDistWeights(Vector3d const & pos){
	return __getDistWeights(pos);
}

tuple8d getDistWeights(Vector3d const & worldSpacePos, Vector3d const & c000Pos, Vector3d const & c111Pos){
	return __getDistWeights(getDiff(worldSpacePos, c000Pos, c111Pos));
}



Vector3d const project(Vector3d const& v1, Vector3d const& v2){ 
	return v1 * (v1.dot(v2) / v1.norm());
}

int intersectRaySphere(Vector3d const& pointPos, Vector3d const& pointVel,
	Vector3d const& sphereCenter, double sphereRad, double &t, Vector3d &newPos) 
{
	Vector3d const& sphereSpacePos = pointPos - sphereCenter;

	double b = sphereSpacePos.dot(pointVel); 
	double c = sphereSpacePos.squaredNorm() - sphereRad * sphereRad; 

	double discr = b*b - c; 

	// A negative discriminant corresponds to ray missing sphere 
	if (discr < 0.0f) return 0;

	// Ray now found to intersect sphere, compute smallest t value of intersection
	t = -b - sqrt(discr); 

	// If t is negative, ray started inside sphere so clamp t to zero 
	//if (t < 0.0f) t = 0.0f; 
	newPos = pointPos + pointVel * t;

	return 1;
}

//void intersectRaySphere_2()
//{
//	// Calculate ray start's offset from the sphere center
//	vecd3 p = s - c;
//	
//	double rSquared = r * r;
//	double p_d = dot(p, d);
//	
//	// The sphere is behind or surrounding the start point.
//	if(p_d > 0 || dot(p, p) < rSquared)
//	 return NO_COLLISION;
//	
//	// Flatten p into the plane passing through c perpendicular to the ray.
//	// This gives the closest approach of the ray to the center.
//	vecd3 a = p - p_d * d;
//	
//	double aSquared = dot(a, a);
//	
//	// Closest approach is outside the sphere.
//	if(aSquared > rSquared)
//	  return NO_COLLISION;
//	
//	// Calculate distance from plane where ray enters/exits the sphere.    
//	double h = sqrt(rSquared - aSquared);
//	
//	// Calculate intersection point relative to sphere center.
//	vecd3 i = a - h * d;
//	
//	vecd3 intersection = c + i;
//	vecd3 normal = i/r;
//	// We've taken a shortcut here to avoid a second square root.
//	// Note numerical errors can make the normal have length slightly different from 1.
//	// If you need higher precision, you may need to perform a conventional normalization.
//	
//	return (intersection, normal);
//}
//
//vecd3 raySphereIntersection(vecd3 rayStart, vecd3 rayEnd, 
//	vecd3 sphereCenter, double sphereRadius)
//{
//	vecd3 rayStartOffset = minus(rayStart, sphereCenter);
//	vecd3 rayEndOffset = minus(rayEnd, sphereCenter);
//
//	vecd3 
//	if 
//	double sphereRadSquared = sphereRadius * sphereRadius;
//
//}

}