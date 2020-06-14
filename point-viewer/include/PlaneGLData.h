#ifndef POINT_GL_DATA_H
#define POINT_GL_DATA_H

#include "Utilities.h"

//=====================================GRID====================================================

struct PlaneGLData{

	PlaneGLData(){}
	PlaneGLData(std::vector<float> * _vertices, GLuint * _monoColourShader, Vector3f _baseColour = {1.f, 1.f, 1.f}) :
	monoColourShader(_monoColourShader),
	baseColour(_baseColour){
		vertices = _vertices;
	}

    GLuint * monoColourShader;

    GLuint verticesSize = 0;
    GLuint verticesBufferObject = 0;
    GLuint verticesArrayObject = 0;

    std::vector<float> * vertices;

    Vector3f baseColour;
};

void generatePlaneVertexData(std::vector<float> & gridVertices, 
	float squareSize, uint x_gridSize, uint z_gridSize);
void generateTile(std::vector<float> & gridVertices);
void generateLine(std::vector<float> & gridVertices);

void initPlaneVAO(PlaneGLData & glData);
void updatePlaneVAO(PlaneGLData const & glData);
void drawPlane(PlaneGLData const & glData, Eigen::Matrix4f & cameraMat);

#endif /* POINT_GL_DATA_H */