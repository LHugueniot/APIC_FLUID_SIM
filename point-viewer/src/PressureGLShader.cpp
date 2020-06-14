#include "PressureGLShader.h"

GLuint compilePressureShaderProgram(){
	return compileShaderProgram(vertexPressureSource, fragmentPressureSource);
}
