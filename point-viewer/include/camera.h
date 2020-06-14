#ifndef CAMERA_H
#define CAMERA_H

#include "Utilities.h"

struct Camera{

	enum Actions{
		ORBIT_LEFT,
		ORBIT_RIGHT,
		ORBIT_UP,
		ORBIT_DOWN,
		ZOOM_IN,
		ZOOM_OUT
	};

	Camera(float _windowWidth, float _windowHeight, 
		float _fov = TO_RAD(50), //In rads 
		float _far = 200.f, float _near = 1.f,
		float _rotationSpeed = 0.05f,
		float _zoomSpeed = 0.05f);
	
	Matrix4f viewMat;

	Matrix4f projMat;

	//Perspective Parameters
	float windowWidth, windowHeight, fov, far, near;

	//Trackball parameters
	float yaw, pitch, zoom;

	float rotationSpeed, zoomSpeed;

	Vector3f target, eye, transformedEye;
};

void setProjMat(Eigen::Matrix4f & projMat,
	float windowWidth, 
    float windowHeight, 
    float fov, 
    float far, float near);
void updateProjMat(Camera & camera);

void setLookAt(Matrix4f & viewMat,
	Vector3f const & position,
	Vector3f const & target,
	Vector3f const & up);
void updateLookAt(Camera & camera);

void rotateCamera(Camera& camera, float rotateAngle);
void pitchCamera(Camera& camera, float pitchAngle);
void zoomCamera(Camera& camera, float zoomAmount);
void moveCamera(Camera& camera, Camera::Actions action);

void updateCamera(Camera& camera);

#endif /* CAMERA_H */
