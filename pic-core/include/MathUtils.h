#pragma once
#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <array>
#include <cmath>
#include <random>
#include <vector>
#include <algorithm>

#include <chrono>

#include <string>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cassert>

#include <Eigen/Eigen>
#include <Eigen/Core>
#include <Eigen/CholmodSupport>

namespace pic
{

#define DEBUG() std::cout<<"BREAK POINT: LINE "<<__LINE__<<" IN "<<__FILE__<<std::endl

using namespace Eigen;
typedef Matrix<size_t, 3, 1> Vector3st;
typedef std::array<double, 2> tuple2d;
typedef std::array<double, 3> tuple3d;
typedef std::array<double, 4> tuple4d;
typedef std::array<double, 8> tuple8d;

typedef std::array<int, 2> tuple2i;
typedef std::array<int, 3> tuple3i;
typedef std::array<int, 4> tuple4i;
typedef std::array<int, 5> tuple5i;
typedef std::array<int, 6> tuple6i;
typedef std::array<int, 7> tuple7i;
typedef std::array<int, 8> tuple8i;

typedef std::array<size_t, 2> tuple2st;
typedef std::array<size_t, 3> tuple3st;
typedef std::array<size_t, 4> tuple4st;
typedef std::array<size_t, 6> tuple6st;
typedef std::array<size_t, 8> tuple8st;

typedef Map<Vector3d> RefVector3d;

void printVec3(Vector3d v);
void printVec2(Vector2d v);

double getDiff(double const x, double const x0, double const x1);
Vector3d getDiff(Vector3d const x, Vector3d const x0, Vector3d const x1);

double linearInterpolation(double const c0, double const c1, double const z);

double bilinearInterpolation(tuple4d const & attr, Vector2d const & pos);
double bilinearInterpolation(tuple4d const & attr, double const y, double const z);

double trilinearInterpolation(tuple8d const & attr, Vector3d const & pos);
double trilinearInterpolation(tuple8d const & attr, double const x, double const y, double const z);

//Weights are gotten from scaled down 
tuple8d getWeights(double const x, double const y, double const z);
tuple8d getWeights(Vector3d const pos);
tuple8d getWeights(Vector3d const pos, Vector3d c000pos, Vector3d c111pos);

Vector3d const project(Vector3d const& v1, Vector3d const& v2);

template<typename T>
decltype(auto)
getFromIdcs(std::vector<T> v, tuple8i i){
	return std::array<T, 8>{v[i[0]], v[i[1]], v[i[2]], v[i[3]], v[i[4]], v[i[5]], v[i[6]], v[i[7]]};
}

template<typename T, size_t N>
std::ostream &operator<<(std::ostream &out, const std::array<T, N>& arr){
    out<<"i{\n";
    for(auto & e : arr)
    	out<<"\n"<<e;
    out<<"\n}\n";
    return out;
}

template <typename T, size_t N1, size_t N2>
decltype(auto)
operator+(const std::array<T, N1>& ob1, const std::array<T, N2>& ob2){
    std::array<T, N1 + N2> ob3;
    for (int i = 0; i < N1; ++i)ob3[i] = ob1[i];
    for (int i = 0; i < N2; ++i)ob3[i + N1] = ob2[i];
    return ob3; 
}


int intersectRaySphere(Vector3d const& pointPos, Vector3d const& pointVel, 
	Vector3d const& sphereCenter, double sphereRad, double &t, Vector3d &q);

static double const GRAV_Y = -9.8f;
static double const waterDensity = 1.f;
}
#endif /* MATH_UTILS_H */
