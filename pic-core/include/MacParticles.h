#pragma once
#ifndef MAC_PARTICLES_H
#define MAC_PARTICLES_H

#include "MathUtils.h"


namespace pic{

//RefVector3d & RefVector3d::operator=(Vector3d const & other){
//    if (this != &other) { // self-assignment check expected
//    	this[0] = other[0];
//    	this[1] = other[1];
//    	this[2] = other[2];
//    }
//    return *this;
//}

struct MacParticles{

	MacParticles(double const _mass, std::vector<Vector3d> _positions, std::vector<Vector3d> _velocities) :
	size(_positions.size()),
	mass(_mass),
	positions(_positions),
	velocities(_velocities){
		assert(_positions.size() == _velocities.size());
	}

	MacParticles(double const _mass, std::vector<double> _positions, std::vector<double> _velocities) :
	size(_positions.size()),
	mass(_mass),
	positionsTest(_positions),
	velocitiesTest(_velocities){
		assert(_positions.size() == _velocities.size());
	}

	RefVector3d getVec3Ref(uint i){
		return RefVector3d(&positionsTest[i * 3]);
	}
	Vector3d getVec3Pos(uint i){
		return {positionsTest[i * 3], positionsTest[i * 3 + 1], positionsTest[i * 3 + 1]};
	}

	void setVec3Pos(uint i, Vector3d const & newPos){
		positionsTest[i * 3] = newPos[0];
		positionsTest[i * 3 + 1] = newPos[1];
		positionsTest[i * 3 + 1] = newPos[2];
	}
	//Matrix<double &, 1, 3> getVec3Ref(uint i){
	//	return {& positionsTest[i * 3], & positionsTest[i * 3 + 1], & positionsTest[i * 3 + 1]}
	//}

	int size;
	double mass;
	std::vector<double> positionsTest;
	std::vector<double> velocitiesTest;

	std::vector<Vector3d> positions;
	std::vector<Vector3d> velocities;
};
}

#endif /* MAC_PARTICLES_H */