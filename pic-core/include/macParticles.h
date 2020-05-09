#pragma once
#ifndef MAC_PARTICLES_H
#define MAC_PARTICLES_H

#include "mathUtils.h"


namespace pic{

	struct MacParticles{

		MacParticles(double const _mass, std::vector<Vector3d> _pos, std::vector<Vector3d> _vel) :
		size(_pos.size()),
		mass(_mass),
		pos(_pos),
		vel(_vel)
		{
			assert(_pos.size() == _vel.size());
			//if (_pos.size() != _vel.size()){
			//	std::string error_msg = "Number of elements in pos and vel is not equal ";
			//	error_msg += _pos.size();
			//	error_msg += " != ";
			//	error_msg += _vel.size();
			//	throw std::runtime_error(error_msg);
			//}
		}

		int size;
		double mass;
		std::vector<Vector3d> positions;
		std::vector<Vector3d> velocities;
	};
}

#endif /* MAC_PARTICLES_H */