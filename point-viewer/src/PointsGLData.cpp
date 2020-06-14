#include "PointsGLData.h"

//=====================================POINTS======================================================

void initPointsVAO(PointGLData & glData){
    glData.verticesSize = glData.vertices->size();
    std::cout<<glData.verticesSize<<std::endl;
    glGenBuffers(1, &glData.verticesBufferObject);
    updatePointsVAO(glData);

    glGenVertexArrays(1, &glData.verticesArrayObject);
    glBindVertexArray(glData.verticesArrayObject);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, glData.verticesBufferObject);
    glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0, NULL);
}

void updatePointsVAO(PointGLData const & glData){
    glBindBuffer(GL_ARRAY_BUFFER, glData.verticesBufferObject);
    glBufferData(GL_ARRAY_BUFFER,
        sizeof(GLdouble) * glData.verticesSize,
        glData.vertices->data(), GL_STATIC_DRAW);
}

void drawPoints(PointGLData const & glData, Eigen::Matrix4f & cameraMat){

    float defaultPointSize;
    glGetFloatv(GL_POINT_SIZE , &defaultPointSize);

    glUseProgram(*glData.monoColourShader);

    glPointSize(3);
    GLuint mvpID = glGetUniformLocation(*glData.monoColourShader, "MVP");
    glUniformMatrix4fv(mvpID, 1, GL_FALSE, cameraMat.data());

    GLuint baseColID = glGetUniformLocation(*glData.monoColourShader, "base_colour");
    glUniform3fv(baseColID, 1, glData.baseColour.data());

    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, glData.verticesBufferObject);
    glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0, NULL);

    glDrawArrays(GL_POINTS, 0, glData.verticesSize/3);

    glDisableVertexAttribArray(0);

    glPointSize(defaultPointSize);
}
