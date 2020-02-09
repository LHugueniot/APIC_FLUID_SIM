#ifndef BASIC_GEOM_H
#define BASIC_GEOM_H

#include "utilities.h"


struct gridData
{
    //Container to access all particle related data
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