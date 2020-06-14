#pragma once
#ifndef PARTICLES_H
#define PARTICLES_H

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

struct Particles{

	Particles(double const _mass,
		std::vector<Vector3d> _positions,
		std::vector<Vector3d> _velocities) :
	num(_positions.size()),
	mass(_mass){
		assert(_positions.size() == _velocities.size());
		positions.resize(_positions.size() * 3);
		velocities.resize(_velocities.size() * 3);
		for(uint i = 0 ; i < _positions.size() ; i++ ){
			getPos(i) = _positions[i];
			getVel(i) = _velocities[i];
		}
	}

	Particles(double const _mass,
		std::vector<double> _positions,
		std::vector<double> _velocities) :
	num(_positions.size() / 3),
	mass(_mass),
	positions(_positions),
	velocities(_velocities){
		assert(_positions.size() == _velocities.size());
		assert(_positions.size() % 3 == 0);
		assert(_velocities.size() % 3 == 0);
	}

	Vector3dRef getPos(uint i){
		return Vector3dRef(&positions[i * 3]);
	}
	Vector3d getPos(uint i) const{
		return Vector3d{positions[i * 3], positions[i * 3 + 1], positions[i * 3 + 2]};
	}


	Vector3dRef getPrevPos(uint i){
		return Vector3dRef(&previousPositions[i * 3]);
	}
	Vector3d getPrevPos(uint i) const{
		return Vector3d{previousPositions[i * 3], previousPositions[i * 3 + 1], previousPositions[i * 3 + 2]};
	}


	Vector3dRef getVel(uint i){
		return Vector3dRef(&velocities[i * 3]);
	}
	Vector3d getVel(uint i) const{
		return Vector3d{velocities[i * 3], velocities[i * 3 + 1], velocities[i * 3 + 2]};
	}

	uint num;
	double mass;
	std::vector<double> positions;
	std::vector<double> previousPositions;
	std::vector<double> velocities;
};

struct AffineParticles : Particles{
	std::vector<double> affine;
};

//struct FlipParticles : public Particles{
//	
//};

}

#endif /* PARTICLES_H */