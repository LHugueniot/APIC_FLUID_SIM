#include "MonoColourGLShader.h"

GLuint compileMonoColourShaderProgram(){
	return compileShaderProgram(vertexMonoColourSource, fragmentMonoColourSource);
}