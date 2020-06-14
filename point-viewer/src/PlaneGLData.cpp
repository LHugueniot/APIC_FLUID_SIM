#include "PlaneGLData.h"

void generateLine(std::vector<float> & planeVertices){
    planeVertices.resize(3 * 2);
    planeVertices[0] = 0;
    planeVertices[1] = 0;
    planeVertices[2] = 0;
    planeVertices[3] = 10;
    planeVertices[4] = 0;
    planeVertices[5] = 0;
}

void generateTile(std::vector<float> & planeVertices){
    planeVertices.resize(3 * 8);
    planeVertices[0] = 0;
    planeVertices[1] = 0;
    planeVertices[2] = 0;
    planeVertices[3] = 10;
    planeVertices[4] = 0;
    planeVertices[5] = 0;

    planeVertices[6] = 10;
    planeVertices[7] = 0;
    planeVertices[8] = 0;
    planeVertices[9] = 10;
    planeVertices[10] = 0;
    planeVertices[11] = 10;

    planeVertices[12] = 10;
    planeVertices[13] = 0;
    planeVertices[14] = 10;
    planeVertices[15] = 0;
    planeVertices[16] = 0;
    planeVertices[17] = 10;

    planeVertices[18] = 0;
    planeVertices[19] = 0;
    planeVertices[20] = 10;
    planeVertices[21] = 0;
    planeVertices[22] = 0;
    planeVertices[23] = 0;
}

void generatePlaneVertexData(std::vector<float> & planeVertices,
    float squareSize, uint planeSize_x, uint planeSize_z){

    uint totalVertices = (planeSize_x + 1) * 6 + (planeSize_z + 1) * 6;
    if(totalVertices != planeVertices.size())
        planeVertices.resize(totalVertices);

    for (uint i = 0 ; i < totalVertices / 2 ; i += 6){
        uint row_col = i / 6;
        planeVertices[i]     = row_col * squareSize;
        planeVertices[i + 1] = 0;
        planeVertices[i + 2] = planeSize_x * squareSize;
        planeVertices[i + 3] = row_col * squareSize;
        planeVertices[i + 4] = 0;
        planeVertices[i + 5] = 0;
    }
    for (uint i = totalVertices / 2 ; i < totalVertices ; i += 6) {
        uint row_col = (i - totalVertices / 2 )/ 6;

        planeVertices[i]     = planeSize_z * squareSize;
        planeVertices[i + 1] = 0;
        planeVertices[i + 2] = row_col * squareSize;
        planeVertices[i + 3] = 0;
        planeVertices[i + 4] = 0;
        planeVertices[i + 5] = row_col * squareSize;
    }
}

void initPlaneVAO(PlaneGLData & glData){

    glData.verticesSize = glData.vertices->size();
    glGenBuffers(1, &glData.verticesBufferObject);
    updatePlaneVAO(glData);

    glGenVertexArrays(1, &glData.verticesArrayObject);
    glBindVertexArray(glData.verticesArrayObject);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, glData.verticesBufferObject);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
}

void updatePlaneVAO(PlaneGLData const & glData){

    glBindBuffer(GL_ARRAY_BUFFER, glData.verticesBufferObject);
    glBufferData(GL_ARRAY_BUFFER,
        sizeof(GLfloat) * glData.verticesSize,
        glData.vertices->data(), GL_STATIC_DRAW);
}

void drawPlane(PlaneGLData const & glData, Eigen::Matrix4f & cameraMat){

    glUseProgram(*glData.monoColourShader);
    GLuint mvpID = glGetUniformLocation(*glData.monoColourShader, "MVP");
    glUniformMatrix4fv(mvpID, 1, GL_FALSE, cameraMat.data());

    GLuint baseColID = glGetUniformLocation(*glData.monoColourShader, "base_colour");
    glUniform3fv(baseColID, 1, glData.baseColour.data());

    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, glData.verticesBufferObject);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    glDrawArrays(GL_LINES, 0, glData.verticesSize/3);
    glDisableVertexAttribArray(0);
}
