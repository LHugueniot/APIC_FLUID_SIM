#include "basicGeom.h"

//=====================================POINTS======================================================

std::vector<float> flattenPointCoorAttr(
    std::vector<float> positions_x,
    std::vector<float> positions_y,
    std::vector<float> positions_z)
{
    auto vertexData = std::vector<float>(positions_x.size() + positions_y.size() + positions_z.size());
    for (int i = 0 ; i < positions_x.size() ; i++)
    {
        auto index_stride = i * 3;
        vertexData[i] = positions_x[i];
        vertexData[i + 1] = positions_y[i];
        vertexData[i + 2] = positions_z[i];
        std::cout<<"x: "<<vertexData[i]<<" y: "<<vertexData[i + 1]<<" z: "<<vertexData[i + 2]<<std::endl;
    }
    return vertexData;
}

void populatePoints(pointData & points,     
    std::vector<float> positions_x,
    std::vector<float> positions_y,
    std::vector<float> positions_z)
{
    points.vertexData =
    flattenPointCoorAttr(positions_x,
                        positions_y,
                        positions_z);

    points.vertexDataSize = points.vertexData.size();

    glGenBuffers(1, &points.vertexBufferObject);

    glBindBuffer(GL_ARRAY_BUFFER, points.vertexBufferObject);

    glBufferData(GL_ARRAY_BUFFER,
        sizeof(GLfloat)  * points.vertexDataSize,
        &points.vertexData[0], GL_STATIC_DRAW);

    glGenVertexArrays(1, &points.vertexArrayObject);
    glBindVertexArray(points.vertexArrayObject);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, points.vertexBufferObject);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
}

void updatePointsVAO(const pointData & points)
{
    glBindBuffer(GL_ARRAY_BUFFER, points.vertexBufferObject);
    glBufferData(GL_ARRAY_BUFFER,
        sizeof(GLfloat) * points.vertexDataSize,
        &points.vertexData[0], GL_STATIC_DRAW);
}

void drawPoints(const pointData & points, GLuint shaderProgram, glm::mat4 cameraMat)
{
    GLuint mvpID = glGetUniformLocation(shaderProgram, "MVP");
    glUniformMatrix4fv(mvpID, 1, GL_FALSE, &cameraMat[0][0]);

    glUseProgram(shaderProgram);
    glBindVertexArray(points.vertexArrayObject);
    glDrawArrays(GL_POINTS, 0, points.vertexDataSize/3);
}

//==========================================GRID===================================================

std::vector<float> generateGridVertexData(float squareSize, uint x_gridSize, uint z_gridSize)
{

    std::vector<float> vertexData((x_gridSize + 1) * 6 + (z_gridSize + 1) * 6);

    for (uint i = 0 ; i < vertexData.size() / 2 ; i += 6)
    {
        uint row_col = i / 6;

        vertexData[i]     = row_col * squareSize;
        vertexData[i + 1] = 0;
        vertexData[i + 2] = x_gridSize * squareSize;
        vertexData[i + 3] = row_col * squareSize;
        vertexData[i + 4] = 0;
        vertexData[i + 5] = 0;
    }

    for (uint i = vertexData.size() /2 ; i < vertexData.size() ; i += 6)
    {
        uint row_col = (i - vertexData.size() /2 )/ 6;

        vertexData[i]     = z_gridSize * squareSize;
        vertexData[i + 1] = 0;
        vertexData[i + 2] = row_col * squareSize;
        vertexData[i + 3] = 0;
        vertexData[i + 4] = 0;
        vertexData[i + 5] = row_col * squareSize;
    }

    return vertexData;
}

void populateGrid(gridData & grid, std::vector<float> vertexData)
{
    grid.vertexData = vertexData;
    grid.vertexDataSize = vertexData.size();

    glGenBuffers(1, &grid.vertexBufferObject);

    glBindBuffer(GL_ARRAY_BUFFER, grid.vertexBufferObject);

    glBufferData(GL_ARRAY_BUFFER,
        sizeof(GLfloat)  * grid.vertexDataSize,
        &grid.vertexData[0], GL_STATIC_DRAW);

    glGenVertexArrays(1, &grid.vertexArrayObject);
    glBindVertexArray(grid.vertexArrayObject);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, grid.vertexBufferObject);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
}


void drawGrid(const gridData & grid, GLuint shaderProgram, glm::mat4 cameraMat)
{
    GLuint mvpID = glGetUniformLocation(shaderProgram, "MVP");
    glUniformMatrix4fv(mvpID, 1, GL_FALSE, &cameraMat[0][0]);

    glUseProgram(shaderProgram);
    glBindVertexArray(grid.vertexArrayObject);
    glDrawArrays(GL_LINES, 0, grid.vertexDataSize/3);
}
 

static const GLfloat basicTriVertexData[] = {
    -1.0f, -1.0f, 0.0f,
     1.0f, -1.0f, 0.0f,
     0.0f,  1.0f, 0.0f,
};

GLuint createBasicGeometryTri()
{
    GLuint vertexBufferObject = 0;

    glGenBuffers(1, &vertexBufferObject);

    glBindBuffer(GL_ARRAY_BUFFER, vertexBufferObject);

    glBufferData(GL_ARRAY_BUFFER,
        sizeof(basicTriVertexData),
        basicTriVertexData, GL_STATIC_DRAW);

    GLuint vertexArrayObject = 0;
    glGenVertexArrays(1, &vertexArrayObject);
    glBindVertexArray(vertexArrayObject);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, vertexBufferObject);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    return vertexBufferObject;
}

void drawBasicGeometryTri(GLuint vertexArrayObject, GLuint shaderProgram, glm::mat4 cameraMat)
{
    GLuint mvpID = glGetUniformLocation(shaderProgram, "MVP");
    glUniformMatrix4fv(mvpID, 1, GL_FALSE, &cameraMat[0][0]);

    glUseProgram(shaderProgram);
    glBindVertexArray(vertexArrayObject);
    glDrawArrays(GL_TRIANGLES, 0, 3);
}