#include "camera.h"

cameraData::cameraData(float windowWidth, float windowHeight)
{
    eye = glm::dvec3(0.f, 0.f, 30.f);
    target = glm::dvec3(0, 0, 0);
    transformedEye = eye;

    yaw = 0.0;
    pitch = 0.0;
    zoom = 1.0;

    projectionMat = glm::perspective(glm::radians(45.0f), (float)windowWidth / (float)windowHeight, 0.1f, 100.0f);
}

static const float rotationSpeed = 0.05f;
static const float zoomSpeed = 0.05f;

void rotateCamera(cameraData& camera, double rotateAngle)
{
    camera.yaw -= rotateAngle;

    if (camera.yaw > glm::pi<double>())
        camera.yaw -= 2.0 * glm::pi<double>();
    else if (camera.yaw < -glm::pi<double>())
        camera.yaw += 2.0 * glm::pi<double>();
}

void pitchCamera(cameraData& camera, double pitchAngle)
{
    camera.pitch = glm::clamp(camera.pitch + pitchAngle,
                       -glm::pi<double>() * 0.5,
                       glm::pi<double>() * 0.5);
}

void zoomCamera(cameraData& camera, double zoomAmount)
{
    camera.zoom = glm::clamp(camera.zoom + zoomAmount, 0.0, 10.0);
}

void moveCamera(cameraData& camera, cameraActions action)
{    
    switch (action)
    {
        case ORBIT_LEFT:
            rotateCamera(camera, rotationSpeed);
            break;
        case ORBIT_RIGHT:
            rotateCamera(camera, -rotationSpeed);
            break;
        case ORBIT_UP:
            pitchCamera(camera, -rotationSpeed);
            break;
        case ORBIT_DOWN:
            pitchCamera(camera, rotationSpeed);
            break;
        case ZOOM_IN:
            zoomCamera(camera, zoomSpeed);
            break;
        case ZOOM_OUT:
            zoomCamera(camera, -zoomSpeed);
            break;
    }
}

void updateCamera(cameraData& camera)
{
    glm::dmat3 R_yaw = glm::mat3_cast(glm::angleAxis(camera.yaw, glm::dvec3(0.0, 1.0, 0.0)));
    glm::dmat3 R_pitch = glm::mat3_cast(glm::angleAxis(camera.pitch, glm::dvec3(1.0, 0.0, 0.0)));
    camera.transformedEye = (R_yaw * R_pitch * (camera.zoom * (camera.eye-camera.target))) + camera.target;
    camera.viewMat = glm::lookAt(glm::vec3(camera.transformedEye), glm::vec3(camera.target), glm::vec3(0.0f,1.0f,0.0f));
}

/*
void moveCamera(cameraData& camera, cameraActions action)
{

	
    float rotationSpeed = 0.05f;

	switch (action)
	{
		case ORBIT_LEFT:
			camera.translationMat = glm::rotate(glm::mat4(1.0f), rotationSpeed, glm::vec3(0.0f, 1.0f, 0.0f)) * camera.translationMat; 
			break;
		case ORBIT_RIGHT:
			camera.translationMat = glm::rotate(glm::mat4(1.0f), -rotationSpeed, glm::vec3(0.0f, 1.0f, 0.0f)) * camera.translationMat; 
			break;
		case PAN_RIGHT:
			camera.pivotPointMat = glm::translate(glm::mat4(1.0f), glm::vec3(0.1f, 0.0f,0.0f)) * camera.pivotPointMat;
			break;
		case PAN_LEFT:
			camera.pivotPointMat = glm::translate(glm::mat4(1.0f), glm::vec3(-0.1f, 0.0f,0.0f)) * camera.pivotPointMat;
			break;
		case ZOOM_IN:
			camera.translationMat = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f,0.1f)) * camera.translationMat;
			break;
		case ZOOM_OUT:
			camera.translationMat = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f,-0.1f)) * camera.translationMat;
			break;
		case FORWARD:
			camera.pivotPointMat = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f,0.1f)) * camera.pivotPointMat;
			break;
		case BACK:
			camera.pivotPointMat = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f,-0.1f)) * camera.pivotPointMat;
			break;
	}

    glm::dmat3 R_yaw = glm::mat3_cast(glm::angleAxis(camera.yaw, glm::dvec3(0.0, 1.0, 0.0)));
    glm::dmat3 R_pitch = glm::mat3_cast(glm::angleAxis(m_pitch, glm::dvec3(1.0, 0.0, 0.0)));
    m_transformedEye = (R_yaw * R_pitch * (m_zoom * (m_eye-m_target))) + m_target;
    m_V = glm::lookAt(glm::vec3(m_transformedEye), glm::vec3(m_target), glm::vec3(0.0f,1.0f,0.0f));

	camera.viewMat = glm::inverse(camera.translationMat * camera.pivotPointMat );
}
*/

//glm::mat4 defaultcameraMatrix(float width, float height)
//{
//	glm::mat4 Projection = glm::perspective(glm::radians(45.0f), (float) width / (float)height, 0.1f, 100.0f);
//	  
//	// Or, for an ortho camera :
//	//glm::mat4 Projection = glm::ortho(-10.0f,10.0f,-10.0f,10.0f,0.0f,100.0f); // In world coordinates
//	  
//	// cameraMatrix matrix
//	glm::mat4 View = glm::lookAt(
//	    glm::vec3(10,10,10), // cameraMatrix is at (4,3,3), in World Space
//	    glm::vec3(0,0,0), // and looks at the origin
//	    glm::vec3(0,1,0)  // Head is up (set to 0,-1,0 to look upside-down)
//	    );
//	  
//	// Model matrix : an identity matrix (model will be at the origin)
//	glm::mat4 Model = glm::mat4(1.0f);
//
//	//glm::mat4 Projection = glm::perspective(glm::radians(45.0f), 4.0f / 3.0f, 0.1f, 100.f);
//	//glm::mat4 View = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, -Translate));
//	//View = glm::rotate(View, Rotate.y, glm::vec3(-1.0f, 0.0f, 0.0f));
//	//View = glm::rotate(View, Rotate.x, glm::vec3(0.0f, 1.0f, 0.0f));
//	//glm::mat4 Model = glm::scale(glm::mat4(1.0f), glm::vec3(0.5f));
//	return Projection * View * Model;
//}