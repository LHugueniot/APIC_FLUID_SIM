#include "Camera.h"

Camera::Camera(float _windowWidth, float _windowHeight, 
    float _fieldOfView, float _far, float _near,
    float _rotationSpeed, float _zoomSpeed) :
        windowWidth(_windowWidth),
        windowHeight(_windowHeight),
        fieldOfView(_fieldOfView),
        far(_far), near(_near),
        rotationSpeed(_rotationSpeed),
        zoomSpeed(_zoomSpeed),
        eye(0.f, 0.f, 30.f),
        target(0, 0, 0){

    transformedEye = eye;
    yaw = 0.f;
    pitch = 0.f;
    zoom = 1.f;
    updateProjMat(*this);
}

void setProjMat(Eigen::Matrix4f & projMat, float windowWidth, 
    float windowHeight, float fieldOfView, float far, float near){

    projMat.setIdentity();
    float aspect = float(windowWidth)/float(windowHeight);
    float theta = fieldOfView*0.5;
    float range = far - near;
    float invtan = 1./tan(theta);

    projMat(0,0) = invtan / aspect;
    projMat(1,1) = invtan;
    projMat(2,2) = -(near + far) / range;
    projMat(3,2) = -1;
    projMat(2,3) = -2 * near * far / range;
    projMat(3,3) = 0;
}

void updateProjMat(Camera & camera){
    setProjMat(camera.projMat, 
        camera.windowWidth,
        camera.windowHeight,
        camera.fieldOfView,
        camera.far,
        camera.near);
}

void rotateCamera(Camera& camera, float rotateAngle)
{
    camera.yaw -= rotateAngle;

    if (camera.yaw > M_PI)
        camera.yaw -= 2.0 * M_PI;
    else if (camera.yaw < -M_PI)
        camera.yaw += 2.0 * M_PI;
}

void pitchCamera(Camera& camera, float pitchAngle)
{
    camera.pitch = std::clamp(camera.pitch + pitchAngle, - .5f * (float)M_PI, .5f * (float)M_PI);
}

void zoomCamera(Camera& camera, float zoomAmount)
{
    camera.zoom = std::clamp(camera.zoom + zoomAmount, 0.f, 10.f);
}

void moveCamera(Camera& camera, Camera::Actions action)
{    
    switch (action)
    {
        case Camera::ORBIT_LEFT:
            rotateCamera(camera, camera.rotationSpeed);
            break;
        case Camera::ORBIT_RIGHT:
            rotateCamera(camera, -camera.rotationSpeed);
            break;
        case Camera::ORBIT_UP:
            pitchCamera(camera, -camera.rotationSpeed);
            break;
        case Camera::ORBIT_DOWN:
            pitchCamera(camera, camera.rotationSpeed);
            break;
        case Camera::ZOOM_IN:
            zoomCamera(camera, camera.zoomSpeed);
            break;
        case Camera::ZOOM_OUT:
            zoomCamera(camera, -camera.zoomSpeed);
            break;
    }
}

void setLookAt(Matrix4f & viewMat, Vector3f const & position, Vector3f const & target, Vector3f const & up){

    Matrix3f R;
    R.col(2) = (position-target).normalized();
    R.col(0) = up.cross(R.col(2)).normalized();
    R.col(1) = R.col(2).cross(R.col(0));
    viewMat.topLeftCorner<3,3>() = R.transpose();
    viewMat.topRightCorner<3,1>() = -R.transpose() * position;
    viewMat(3,3) = 1.0f;
}

void updateLookAt(Camera & camera){
    setLookAt(camera.viewMat, camera.eye, camera.target, Vector3f(0, 1, 0));
}

void updateCamera(Camera& camera){

    Matrix3f R_yaw; 
    R_yaw = AngleAxisf(camera.yaw, Vector3f::UnitY());
    Matrix3f R_pitch;
    R_pitch = AngleAxisf(camera.yaw, Vector3f::UnitX());

    camera.transformedEye = (R_yaw * R_pitch * (camera.zoom * (camera.eye - camera.target))) + camera.target;
    updateLookAt(camera);
}

/*

void updateCamera(Camera& camera){

    glm::dmat3 R_yaw = glm::mat3_cast(glm::angleAxis(camera.yaw, Vector3f(0.0, 1.0, 0.0)));
    glm::dmat3 R_pitch = glm::mat3_cast(glm::angleAxis(camera.pitch, Vector3f(1.0, 0.0, 0.0)));
    camera.transformedEye = (R_yaw * R_pitch * (camera.zoom * (camera.eye-camera.target))) + camera.target;
    camera.viewMat = glm::lookAt(glm::vec3(camera.transformedEye), glm::vec3(camera.target), glm::vec3(0.0f,1.0f,0.0f));
}

void moveCamera(Camera& camera, cameraActions action)
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

    glm::dmat3 R_yaw = glm::mat3_cast(glm::angleAxis(camera.yaw, Vector3f(0.0, 1.0, 0.0)));
    glm::dmat3 R_pitch = glm::mat3_cast(glm::angleAxis(m_pitch, Vector3f(1.0, 0.0, 0.0)));
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