#ifndef CAMERA_H
#define CAMERA_H

#include "utilities.h"

enum cameraActions
{
	ORBIT_LEFT,
	ORBIT_RIGHT,
	ORBIT_UP,
	ORBIT_DOWN,
	ZOOM_IN,
	ZOOM_OUT
};

struct cameraData
{
	cameraData(float windowWidth, float windowHeight);
	
	glm::mat4 viewMat;
	glm::mat4 projectionMat;

	double yaw, pitch, zoom;
	
	glm::dvec3 target, eye, transformedEye;
};

void rotateCamera(cameraData& camera, double rotateAngle);
void pitchCamera(cameraData& camera, double pitchAngle);
void zoomCamera(cameraData& camera, double zoomAmount);

void updateCamera(cameraData& camera);
void moveCamera(cameraData& camera, cameraActions action);

#endif /* CAMERA_H */
