#include "BasicGrid.h"

void generateGridVertexData(std::vector<double> & gridVertices,
    double squareSize, uint x_gridSize, uint z_gridSize){

    uint totalVertices = (x_gridSize + 1) * 6 + (z_gridSize + 1) * 6;
    if(totalVertices != gridVertices.size())
        gridVertices.resize(totalVertices);

    for (uint i = 0 ; i < totalVertices / 2 ; i += 6){
        uint row_col = i / 6;
        gridVertices[i]     = row_col * squareSize;
        gridVertices[i + 1] = 0;
        gridVertices[i + 2] = x_gridSize * squareSize;
        gridVertices[i + 3] = row_col * squareSize;
        gridVertices[i + 4] = 0;
        gridVertices[i + 5] = 0;
    }
    for (uint i = totalVertices /2 ; i < totalVertices ; i += 6) {
        uint row_col = (i - totalVertices /2 )/ 6;

        gridVertices[i]     = z_gridSize * squareSize;
        gridVertices[i + 1] = 0;
        gridVertices[i + 2] = row_col * squareSize;
        gridVertices[i + 3] = 0;
        gridVertices[i + 4] = 0;
        gridVertices[i + 5] = row_col * squareSize;
    }
}

void initGridVAO(GridGLData & grid){

    grid.verticesSize = grid.vertices->size();
    glGenBuffers(1, &grid.verticesBufferObject);
    glBindBuffer(GL_ARRAY_BUFFER, grid.verticesBufferObject);
    glBufferData(GL_ARRAY_BUFFER,
        sizeof(GLdouble)  * grid.verticesSize,
        &grid.vertices[0],
        GL_STATIC_DRAW);

    glGenVertexArrays(1, &grid.verticesArrayObject);
    glBindVertexArray(grid.verticesArrayObject);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, grid.verticesBufferObject);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
}

void updateGridVAO(GridGLData const & grid){
    glBindBuffer(GL_ARRAY_BUFFER, grid.verticesBufferObject);
    glBufferData(GL_ARRAY_BUFFER,
        sizeof(GLdouble) * grid.verticesSize,
        &grid.vertices[0], GL_STATIC_DRAW);
}

void drawGrid(GridGLData const & grid, GLuint shaderProgram, Eigen::Matrix4f & cameraMat){

    GLuint mvpID = glGetUniformLocation(shaderProgram, "MVP");
    glUniformMatrix4fv(mvpID, 1, GL_FALSE, cameraMat.data());

    glUseProgram(shaderProgram);
    glBindVertexArray(grid.verticesArrayObject);
    glDrawArrays(GL_LINES, 0, grid.verticesSize/3);
}
