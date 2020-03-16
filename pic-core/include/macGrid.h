#include "mathUtils.h"

namespace pic{

template<typename T, int u, int v, int w>
struct MacGrid{

	MacGrid(double origin_x, double origin_y, double origin_z
		double opposite_x, double opposite_y, double opposite_z);

	vec3<int> getIndexFromPos() 

	T cellVel_x[u];
	T cellVel_y[v];
	T cellVel_z[w];
}

}