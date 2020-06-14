#include "MultiColourGLShader.h"

GLuint compileMultiColourShaderProgram(){
	return compileShaderProgram(vertexMultiColourSource, fragmentMultiColourSource);
}
