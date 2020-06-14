#include "Collisions.h"

namespace pic {

inline Vector3d projectPlaneCollision(Vector3d const & planePos,
	Vector3d const & planeNorm, Vector3d & pos, Vector3d & vel){

	Vector3d rayDir = vel;
	double d = planePos.dot(-planeNorm);

	double t = -(d + (pos[0] * planeNorm[0] + pos[1] * planeNorm[1] + pos[2] * planeNorm[2])) 
		/ ((rayDir[0] * planeNorm[0] + rayDir[1] * planeNorm[1] + rayDir[2] * planeNorm[2]));
	//double t = -(d + (pos * planeNorm).sum()) / ((rayDir * planeNorm).sum());
	vel = vel - 2 * (vel.dot(planeNorm)) * planeNorm;
	pos = pos + t * rayDir;
}

bool boxCollisionBounce(Vector3d const & c000, 	//Minimum corner
 						Vector3d const & c111, 	//Maximum corner
 						Vector3d & prev, 	   	//position
 						Vector3d & vel){		//velocity

	Vector3d next = prev + vel;
/*
	bool prio_x = false;
	bool prio_y = false;
	bool prio_z = false;

	double distToFace_x = std::min(c111[0] - prev[0], c000[0] - prev[0]);
	double distToFace_y = std::min(c111[1] - prev[1], c000[1] - prev[1]);
	double distToFace_z = std::min(c111[2] - prev[2], c000[2] - prev[2]);

	prio_x = distToFace_x ?
*/

	//x collisions
	if(c111[1] > prev[1] && prev[1] > c000[1] && 
		c111[2] > prev[2] && prev[2] > c000[2]){

		if(prev[0] > c000[0] && next[0] < c000[0]){
			projectPlaneCollision(c000, Vector3d(1.f, 0.f, 0.f), prev, vel);
			return true;
		}
		else if(next[0] > c111[0] && c111[0] > prev[0]){
			projectPlaneCollision(c111, Vector3d(-1.f, 0.f, 0.f), prev, vel);
			return true;
		}
	}

	//y collisions
	if(c111[0] > prev[0] && prev[0] > c000[0] && 
		c111[2] > prev[2] && prev[2] > c000[2]){

		if(prev[1] > c000[1] && next[1] < c000[1]){
			projectPlaneCollision(c000, Vector3d(0.f, 1.f, 0.f), prev, vel);
			return true;
		}
		else if(next[1] > c111[1] && c111[1] > prev[1]){
			projectPlaneCollision(c111, Vector3d(0.f, -1.f, 0.f), prev, vel);
			return true;
		}
	}
	//z collisions
	if(c111[0] > prev[0] && prev[0] > c000[0] && 
		c111[1] > prev[1] && prev[1] > c000[1]){

		if(prev[2] > c000[2] && next[2] < c000[2]){
			projectPlaneCollision(c000, Vector3d(0.f, 0.f, 1.f), prev, vel);
			return true;
		}
		else if(next[2] > c111[2] && c111[2] > prev[2]){
			projectPlaneCollision(c111, Vector3d(0.f, 0.f, -1.f), prev, vel);
			return true;
		}
	}
	return false;
}

/*
bool boxCollision(Vector3d const & c000, //Minimum corner
 				  Vector3d const & c111, //Maximum corner
 				  Vector3d & prev, 	   	 //position
 				  Vector3d & vel){		 //velocity

	Vector3d next = prev + vel;

	//x collisions
	if(c111[1] > prev[1] && prev[1] > c000[1] && 
		c111[2] > prev[2] && prev[2] > c000[2]){

		if(prev[0] > c000[0] && next[0] < c000[0]){
			projectPlaneCollision(c000, Vector3d(1.f, 0.f, 0.f), prev, vel);
			return true;
		}
		else if(next[0] > c111[0] && c111[0] > prev[0]){
			projectPlaneCollision(c111, Vector3d(-1.f, 0.f, 0.f), prev, vel);
			return true;
		}
	}

	//y collisions
	if(c111[0] > prev[0] && prev[0] > c000[0] && 
		c111[2] > prev[2] && prev[2] > c000[2]){

		if(prev[1] > c000[1] && next[1] < c000[1]){
			projectPlaneCollision(c000, Vector3d(0.f, 1.f, 0.f), prev, vel);
			return true;
		}
		else if(next[1] > c111[1] && c111[1] > prev[1]){
			projectPlaneCollision(c111, Vector3d(0.f, -1.f, 0.f), prev, vel);
			return true;
		}
	}
	//z collisions
	if(c111[0] > prev[0] && prev[0] > c000[0] && 
		c111[1] > prev[1] && prev[1] > c000[1]){

		if(prev[2] > c000[2] && next[2] < c000[2]){
			projectPlaneCollision(c000, Vector3d(0.f, 0.f, 1.f), prev, vel);
			return true;
		}
		else if(next[2] > c111[2] && c111[2] > prev[2]){
			projectPlaneCollision(c111, Vector3d(0.f, 0.f, -1.f), prev, vel);
			return true;
		}
	}
	return false;
}
*/

bool boxCollision(Vector3d const & c000, 	//Minimum corner
 				  Vector3d const & c111, 	//Maximum corner
 				  Vector3d & prev, 	   	//position
 				  Vector3d & vel){		//velocity

	Vector3d next = prev + vel;

	//x collisions
	if(c111[1] > prev[1] && prev[1] > c000[1] && 
		c111[2] > prev[2] && prev[2] > c000[2]){

		if((prev[0] > c000[0] && next[0] < c000[0]) 
			|| (next[0] > c111[0] && c111[0] > prev[0])){
			double velDissipate = std::abs(vel[0] * .5f);
			vel[0] = 0;
			vel[1] = vel[1] > 0 ?  vel[1] + velDissipate : vel[1] - velDissipate;
			vel[2] = vel[2] > 0 ?  vel[2] + velDissipate : vel[2] - velDissipate;
			//vel = {0, vel[1] + vel[0] * .5f, vel[2] + vel[0] * .5f};
			return true;
		}
	}

	//y collisions
	if(c111[0] > prev[0] && prev[0] > c000[0] && 
		c111[2] > prev[2] && prev[2] > c000[2]){

		if((prev[1] > c000[1] && next[1] < c000[1]) 
			|| (next[1] > c111[1] && c111[1] > prev[1])){
			double velDissipate = std::abs(vel[1] * .5f);
			vel[0] = vel[0] > 0 ?  vel[0] + velDissipate : vel[0] - velDissipate;
			vel[1] = 0;
			vel[2] = vel[2] > 0 ?  vel[2] + velDissipate : vel[2] - velDissipate;
			//vel = {vel[0] + vel[1] * .5f, 0, vel[2] + vel[1] * .5f};
			return true;
		}
	}
	//z collisions
	if(c111[0] > prev[0] && prev[0] > c000[0] && 
		c111[1] > prev[1] && prev[1] > c000[1]){

		if((prev[2] > c000[2] && next[2] < c000[2])
			|| (next[2] > c111[2] && c111[2] > prev[2])){
			double velDissipate = std::abs(vel[2] * .5f);
			vel[0] = vel[0] > 0 ?  vel[0] + velDissipate : vel[0] - velDissipate;
			vel[1] = vel[1] > 0 ?  vel[1] + velDissipate : vel[1] - velDissipate;
			vel[2] = 0;

			//vel = {vel[0] + vel[2] * .5f, vel[1] + vel[2] * .5f, 0};
			return true;
		}
	}
	return false;
}
}