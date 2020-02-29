#ifndef GRID_H
#define GRID_H

#include "particleSystem.h"

struct cellGridData
{
	float x_fieldPosition, y_fieldPosition, z_fieldPosition;
	float x_fieldSize, y_fieldSize, z_fieldSize;  
	uint i_dimension, j_dimension, k_dimension;

	//Container to access all particle related data
	std::vector<float> cellVertexData;
	std::vector<float> cellVelocityData;
	std::vector<float> cellMassData;
};

std::vector<float> createCellGridVertexData(float x_fieldPosition, float y_fieldPosition, float z_fieldPosition,
										float x_fieldSize, float y_fieldSize, float z_fieldSize,  
										uint i_dimension, uint j_dimension, uint k_dimension);

//std::vector<std::vector<std::vector<float>>> createGridFromDimensions(float x_fieldPosition, float y_fieldPosition, float z_fieldPosition, 
//																	  float x_fieldSize, float y_fieldSize, float z_fieldSize);

#endif /* GRID_H */