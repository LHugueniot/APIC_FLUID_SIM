#ifndef POINTS_GL_DATA_H
#define POINTS_GL_DATA_H

#include "Utilities.h"

//=====================================POINTS==================================================

struct PointGLData{

	PointGLData(){}
	PointGLData(std::vector<double> * _vertices, GLuint * _monoColourShader, Vector3f _baseColour = {0.f, 0.f, 1.f}) :
	monoColourShader(_monoColourShader),
	baseColour(_baseColour){
		vertices = _vertices;
	}

 	GLuint * monoColourShader;

    GLuint verticesSize = 0;
    GLuint verticesBufferObject = 0;
    GLuint verticesArrayObject = 0;

    std::vector<double> * vertices;

    Vector3f baseColour;
};

void initPointsVAO(PointGLData & glData);
void updatePointsVAO(PointGLData const & glData);
void drawPoints(PointGLData const & glData, Matrix4f & cameraMat);

#endif /* POINTS_GL_DATA_H */