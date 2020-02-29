#include "particleSystem.h"

std::vector<float> generateRandomParticleData(float x_fieldPosition, float y_fieldPosition, float z_fieldPosition, 
                                                   float x_fieldSize, float y_fieldSize, float z_fieldSize, 
                                                   size_t particleSystemSize)
{
    std::vector<float> particleSystemVertexData(particleSystemSize * 3);

    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_real_distribution<> x_dist(x_fieldPosition, x_fieldPosition + x_fieldSize);
    std::uniform_real_distribution<> y_dist(y_fieldPosition, y_fieldPosition + y_fieldSize);
    std::uniform_real_distribution<> z_dist(z_fieldPosition, z_fieldPosition + z_fieldSize);

    for (size_t i = 0 ; i < particleSystemSize * 3 ; i += 3)
    {
        particleSystemVertexData[i] = x_dist(eng);
        particleSystemVertexData[i + 1] = y_dist(eng);
        particleSystemVertexData[i + 2] = z_dist(eng);
    }
    return particleSystemVertexData;
}

void populateParticleSystem(particleSystemData & particleSystem, std::vector<float> vertexData)
{
    particleSystem.vertexData = vertexData;
    particleSystem.vertexDataSize = particleSystem.vertexData.size();
    particleSystem.velocityData = std::vector<float>(particleSystem.vertexData.size());
    particleSystem.massData = std::vector<float>(particleSystem.vertexData.size()/3, 1);

    glGenBuffers(1, &particleSystem.vertexBufferObject);

    glBindBuffer(GL_ARRAY_BUFFER, particleSystem.vertexBufferObject);

    glBufferData(GL_ARRAY_BUFFER,
        sizeof(GLfloat) * particleSystem.vertexDataSize,
        &particleSystem.vertexData[0], GL_STATIC_DRAW);

    glGenVertexArrays(1, &particleSystem.vertexArrayObject);
    glBindVertexArray(particleSystem.vertexArrayObject);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, particleSystem.vertexBufferObject);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
}

void randomiseParticleStep(particleSystemData & particleSystem)
{
    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_real_distribution<> x_dist(0, 0.01);
    std::uniform_real_distribution<> y_dist(0, 0.01);
    std::uniform_real_distribution<> z_dist(0, 0.01);

    for(uint i = 0 ; i < particleSystem.vertexDataSize ; i += 3)
    {
        particleSystem.vertexData[i] += x_dist(eng);
        particleSystem.vertexData[i+1] += y_dist(eng);
        particleSystem.vertexData[i+2] += z_dist(eng);
    }
}

void updateParticleSystemVAO(const particleSystemData & particleSystem)
{
    glBindBuffer(GL_ARRAY_BUFFER, particleSystem.vertexBufferObject);
    glBufferData(GL_ARRAY_BUFFER,
        sizeof(GLfloat) * particleSystem.vertexDataSize,
        &particleSystem.vertexData[0], GL_STATIC_DRAW);
}

void drawParticleSystem(const particleSystemData & particleSystem, GLuint shaderProgram, glm::mat4 cameraMat)
{
    GLuint mvpID = glGetUniformLocation(shaderProgram, "MVP");
    glUniformMatrix4fv(mvpID, 1, GL_FALSE, &cameraMat[0][0]);

    glUseProgram(shaderProgram);
    glBindVertexArray(particleSystem.vertexArrayObject);
    glDrawArrays(GL_POINTS, 0, particleSystem.vertexDataSize/3);
}
