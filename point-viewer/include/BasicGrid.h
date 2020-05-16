#ifndef BASIC_GRID_H
#define BASIC_GRID_H

#include "Utilities.h"

//=====================================GRID====================================================

struct GridGLData{

	GridGLData(){}
	GridGLData(std::vector<double> & _vertices){
		vertices = & _vertices;
	}

    GLuint verticesSize = 0;
    GLuint verticesBufferObject = 0;
    GLuint verticesArrayObject = 0;

    std::vector<double> * vertices;
};

void generateGridVertexData(std::vector<double> & gridVertices, 
	double squareSize, uint x_gridSize, uint z_gridSize);

void initGridVAO(GridGLData & grid);
void updateGridVAO(GridGLData const & grid);
void drawGrid(GridGLData const & grid, GLuint shaderProgram, Eigen::Matrix4f & cameraMat);

#endif /* BASIC_GRID_H */