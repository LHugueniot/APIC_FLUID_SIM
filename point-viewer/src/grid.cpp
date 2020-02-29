#include "grid.h"

std::vector<float> createCellGridVertexData(float x_fieldPosition, float y_fieldPosition, float z_fieldPosition,
										float x_fieldSize, float y_fieldSize, float z_fieldSize,  
										uint i_dimension, uint j_dimension, uint k_dimension)
{
	std::vector<float> gridVertexData(i_dimension * j_dimension * k_dimension * 3);

	float i_cellSize = x_fieldSize / i_dimension;
	float j_cellSize = y_fieldSize / j_dimension;
	float k_cellSize = z_fieldSize / k_dimension;

	float i_halfCellLen = i_cellSize / 2;
	float j_halfCellLen = j_cellSize / 2;
	float k_halfCellLen = k_cellSize / 2;

	for (size_t i = 0 ; i < i_dimension ; i++ )
	{
		for (size_t j = 0 ; j < j_dimension ; j++ )
		{
			for (size_t k = 0 ; k < k_dimension ; k++ )
			{
				size_t index = (k * (j_dimension * i_dimension) + j * i_dimension + i) * 3 ;

				gridVertexData[index] 	  = i * i_cellSize + i_halfCellLen + x_fieldPosition;
				gridVertexData[index + 1] = j * j_cellSize + j_halfCellLen + y_fieldPosition;
				gridVertexData[index + 2] = k * k_cellSize + k_halfCellLen + z_fieldPosition;
			}
		}
	}

	return gridVertexData;
}

void particleToCellIndex(float x_pos, float y_pos, float z_pos, cellGridData & grid)
{
	//grid
	//i * i_cellSize + i_halfCellLen + x_fieldPosition
}

size_t getIndex(uint i_dimension, uint j_dimension, uint k_dimension)
{
	return (k * (j_dimension * i_dimension) + j * i_dimension + i) * 3;
}

void PICStep(particleSystemData & particleSystem, cellGridData & grid)
{

    for(uint i = 0 ; i < particleSystem.vertexDataSize ; i += 3)
    {
        auto x = particleSystem.vertexData[i];
        auto y = particleSystem.vertexData[i+1];
        auto z = particleSystem.vertexData[i+2];


    }
}
