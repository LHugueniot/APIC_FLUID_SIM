#ifndef BASIC_GEOM_H
#define BASIC_GEOM_H

#include "utilities.h"

//=====================================POINTS======================================================

struct pointData
{
	std::vector<float> vertexData;

    GLuint vertexDataSize = 0;
    GLuint vertexBufferObject = 0;
    GLuint vertexArrayObject = 0;
};

std::vector<float> flattenPointCoorAttr(std::vector<float> positions_x, std::vector<float> positions_y, std::vector<float> positions_z);
void populatePoints(pointData & points, std::vector<float> positions_x, std::vector<float> positions_y,std::vector<float> positions_z);
void updatePointsVAO(const pointData & points);
void drawPoints(const pointData & points, GLuint shaderProgram, glm::mat4 cameraMat);

//==========================================GRID===================================================

struct gridData
{
    std::vector<float> vertexData;

    GLuint vertexDataSize = 0;
    GLuint vertexBufferObject = 0;
    GLuint vertexArrayObject = 0;
};

std::vector<float> generateGridVertexData(float squareSize, uint x_gridSize, uint z_gridSize);
void populateGrid(gridData & grid, std::vector<float> vertexData);
//void updateGridVAO(const particleSystemData & particleSystem);
void drawGrid(const gridData & grid, GLuint shaderProgram, glm::mat4 cameraMat);


GLuint createBasicGeometryTri();
void drawBasicGeometryTri(const gridData & grid, GLuint shaderProgram, glm::mat4 cameraMat);

#endif /* BASIC_GEOM_H */