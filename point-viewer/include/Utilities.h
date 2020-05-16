#ifndef UTILITIES_H
#define UTILITIES_H

#define _USE_MATH_DEFINES

//std lib includes
#include <map>
#include <string>
#include <memory>
#include <vector>
#include <iostream>
#include <random>
#include <chrono>
#include <cmath>
#include <algorithm>

//OpenGL includes
#include <GL/glew.h>

//Eigen includes
#include <Eigen/Eigen>
#include <Eigen/Core>
#include <Eigen/Geometry>
//#include <Eigen/CholmodSupport>

using namespace Eigen;

#define M_RAD_RATIO 0.01745329251
#define TO_RAD(deg) deg * M_RAD_RATIO

typedef unsigned int uint;

#endif /* UTILITIES_H */