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
#include <system_error>

#include <Eigen/Dense>
#include <Eigen/Eigen>
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/CholmodSupport>
#include <Eigen/Geometry>

namespace pic{

#define EZ_MODE 0

#if EZ_MODE
#   define DEBUG() std::cout<<"BREAK POINT: LINE "<<__LINE__<<" IN "<<__FILE__<<std::endl;
#   define DEBUG_VAR(var) std::cout<<"BREAK POINT: LINE "<<__LINE__<<" IN "<<__FILE__<<" FOR VAR: "<<#var<<" = "<<var<<std::endl;
#   define DEBUG_MSG(msg) std::cout<<"BREAK POINT: LINE "<<__LINE__<<" IN "<<__FILE__<<" ~ WITH MSG: "<<msg<<std::endl;
#   define DEBUG_ASS(Expr, Msg) \
    __M_Assert(#Expr, Expr, __FILE__, __LINE__, Msg)
#else
#   define DEBUG(...) do {} while (0)
#   define DEBUG_VAR(var) do {} while (0)
#   define DEBUG_MSG(msg) do {} while (0)
#   define DEBUG_ASS(Expr, Msg)
#endif

template<typename T>
void __M_Assert(const char* expr_str, bool expr, const char* file, int line, T msg){
    if (!expr){
        std::cerr << "Assert failed:\t" << msg << "\n"
            << "Expected:\t" << expr_str << "\n"
            << "Source:\t\t" << file << ", line " << line << "\n";
        abort();
    }
}

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

typedef Map<Vector3d> Vector3dRef;

void printVec3(Vector3d v);
void printVec2(Vector2d v);


template<typename T, size_t N>
std::ostream &operator<<(std::ostream &out, const std::array<T, N>& arr){
    out<<"array<"<<N<<">{";
    for(auto & e : arr)
        out<<"\n"<<e;
    out<<"\n}\n";
    return out;
}

template<typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T>& vec){
    out<<"vector.size()="<<vec.size()<<"{\n";
    for(auto & e : vec)
        out<<"\n"<<e;
    out<<"\n}\n";
    return out;
}

double getDiff(double const x, double const x0, double const x1);
Vector3d getDiff(Vector3d const x, Vector3d const x0, Vector3d const x1);

double linearInterpolation(double const c0, double const c1, double const z);

double bilinearInterpolation(tuple4d const & attr, Vector2d const & pos);
double bilinearInterpolation(tuple4d const & attr, double const y, double const z);

double trilinearInterpolation(tuple8d const & attr, Vector3d const & pos);
double trilinearInterpolation(tuple8d const & attr, double const x, double const y, double const z);

//Weights are gotten from scaled down 
tuple8d getWeights(double const x, double const y, double const z);
tuple8d getWeights(Vector3d const & pos);
tuple8d getWeights(Vector3d const & worldSpacePos, Vector3d const & c000Pos, Vector3d const & c111Pos);

//tuple8d getDistWeights(double const x, double const y, double const z);
tuple8d getDistWeights(Vector3d const pos);
tuple8d getDistWeights(Vector3d const & worldSpacePos, Vector3d const & c000Pos, Vector3d const & c111Pos);

Vector3d const project(Vector3d const& v1, Vector3d const& v2);

template<typename T>
decltype(auto)
getFromIdcs(std::vector<T> v, tuple8i i){
	return std::array<T, 8>{v[i[0]], v[i[1]], v[i[2]], v[i[3]], v[i[4]], v[i[5]], v[i[6]], v[i[7]]};
}

template <typename T, size_t N1, size_t N2>
decltype(auto)
operator+(const std::array<T, N1>& ob1, const std::array<T, N2>& ob2){
    std::array<T, N1 + N2> ob3;
    for (int i = 0; i < N1; ++i)ob3[i] = ob1[i];
    for (int i = 0; i < N2; ++i)ob3[i + N1] = ob2[i];
    return ob3; 
}


template<typename T>
void checkForNan(std::vector<T> v){
    for (auto & i : v)
        assert(!std::isnan(i));
}

template<typename T>
void checkForNan(Matrix<T, Dynamic, 1> M){
    for (int i = 0; i < M.rows(); ++i)
        assert(!std::isnan(M(i, 0)));
}


template<typename T>
void checkForNan(SparseMatrix<T> M){
    for (int i = 0; i < M.rows(); ++i)
        for (int j = 0; j < M.cols(); ++j)
            assert(!std::isnan(M.coeffRef(i, j)));
}

int intersectRaySphere(Vector3d const& pointPos, Vector3d const& pointVel, 
	Vector3d const& sphereCenter, double sphereRad, double &t, Vector3d &q);
/**
template<typename G, void(G::*F)()>
void gridPrint(G& grid, uint ni, uint nj, uint nk){

    for ( uint i = 0 ; i < ni ; i++){
        for ( uint j = 0 ; j < nj ; j++){
            for ( uint k = 0 ; k < nk ; k++){
                std::cout<<(grid.*F)(i, j, k)<<" ";
            }
            std::cout<<std::endl;
        }
        std::cout<<std::endl;
    }
}
**/
static double const GRAV_Y = -9.8f * 3;
//static double const GRAV_Y = 0;
static double const waterDensity = 1.f;
static double const invWaterDensity = 1.f/waterDensity;
static double const airPressure = 1.f;
static double const FPIC_RATIO = 0.05f;

}
#endif /* MATH_UTILS_H */
