#include "Utilities.h"

std::array<double, 3> blueToRedInterpolation(double scalar){
	std::array<double, 3> rgb{ 0,0,0 };

	scalar = scalar/10.f;

	if (scalar <= 0.5f)
	{
	    scalar *= 2.0f;
	    rgb[0] = scalar + 0.5f;
	    rgb[1] = 1.0f - scalar + 0.5f;
	}
	else
	{
	    scalar = scalar * 2.0f - 1.0f;
	    rgb[1] = scalar + 0.5f;
	    rgb[2] = 1.0f - scalar + 0.5f;
	}
	return rgb;
}
