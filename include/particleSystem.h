#ifndef PARTICLE_SYSTEM_H
#define PARTICLE_SYSTEM_H

#include "utilities.h"

struct particleSystemData
{
	//Container to access all particle related data
	std::vector<float> vertexData;
	std::vector<float> velocityData;
	std::vector<float> massData;

	GLuint vertexDataSize = 0;
	GLuint vertexBufferObject = 0;
	GLuint vertexArrayObject = 0;
};

std::vector<float> generateRandomParticleData(float x_fieldPosition, float y_fieldPosition, float z_fieldPosition, 
												   float x_fieldSize, float y_fieldSize, float z_fieldSize, 
												   size_t particleSystemSize);
void populateParticleSystem(particleSystemData & particleSystem, std::vector<float> vertexData);
void randomiseParticleStep(particleSystemData & particleSystem);
void updateParticleSystemVAO(const particleSystemData & particleSystem);
void drawParticleSystem(const particleSystemData & particleSystem, GLuint shaderProgram, glm::mat4 cameraMat);

#endif /* PARTICLE_SYSTEM_H */