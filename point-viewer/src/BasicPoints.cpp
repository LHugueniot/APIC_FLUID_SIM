#include "BasicPoints.h"

//=====================================POINTS======================================================

void initPointsVAO(PointGLData & points){

    points.verticesSize = points.vertices->size();
    glGenBuffers(1, &points.verticesBufferObject);
    glBindBuffer(GL_ARRAY_BUFFER, points.verticesBufferObject);
    glBufferData(GL_ARRAY_BUFFER,
        sizeof(GLdouble)  * points.verticesSize,
        &points.vertices[0],
        GL_STATIC_DRAW);

    glGenVertexArrays(1, &points.verticesArrayObject);
    glBindVertexArray(points.verticesArrayObject);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, points.verticesBufferObject);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
}

void updatePointsVAO(PointGLData const & points){
    glBindBuffer(GL_ARRAY_BUFFER, points.verticesBufferObject);
    glBufferData(GL_ARRAY_BUFFER,
        sizeof(GLdouble) * points.verticesSize,
        &points.vertices[0], GL_STATIC_DRAW);
}

void drawPoints(PointGLData const & points, GLuint shaderProgram, Eigen::Matrix4f & cameraMat){

    GLuint mvpID = glGetUniformLocation(shaderProgram, "MVP");
    glUniformMatrix4fv(mvpID, 1, GL_FALSE, cameraMat.data());

    glUseProgram(shaderProgram);
    glBindVertexArray(points.verticesArrayObject);
    glDrawArrays(GL_POINTS, 0, points.verticesSize/3);
}
