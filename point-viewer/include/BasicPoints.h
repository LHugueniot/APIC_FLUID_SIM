#ifndef BASIC_POINTS_H
#define BASIC_POINTS_H

#include "Utilities.h"

//=====================================POINTS==================================================

struct PointGLData{

	PointGLData(){}
	PointGLData(std::vector<double> & _vertices){
		vertices = & _vertices;
	}

    GLuint verticesSize = 0;
    GLuint verticesBufferObject = 0;
    GLuint verticesArrayObject = 0;

    std::vector<double> * vertices;
};

void initPointsVAO(PointGLData & points);
void updatePointsVAO(PointGLData const & points);
void drawPoints(PointGLData const & points, GLuint shaderProgram, Matrix4f & cameraMat);

#endif /* BASIC_POINTS_H */