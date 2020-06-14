#pragma once
#ifndef COLLISIONS_H
#define COLLISIONS_H

#include "MathUtils.h"

namespace pic{
	bool boxCollisionBounce(Vector3d const & c000,
					  Vector3d const & c111,
					  Vector3d & prev, 
					  Vector3d & vel);
	bool boxCollision(Vector3d const & c000,
					  Vector3d const & c111,
					  Vector3d & prev, 
					  Vector3d & vel);
}

#endif /* COLLISIONS_H */