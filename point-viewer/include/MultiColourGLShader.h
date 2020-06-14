#ifndef MULTI_COLOUR_GL_SHADER_H
#define MULTI_COLOUR_GL_SHADER_H

#include "GLShader.h"

static const char* vertexMultiColourSource =
"#version 330\n"
"uniform mat4 MVP;\n"
"layout(location = 0) in vec3 vertex_position;\n"
"layout(location = 1) in vec3 vertex_color;\n"
"out vec4 frag_colour;\n"
"void main() {\n"
"	frag_colour = vec4(vertex_color, 1.0);\n"
"  	gl_Position = MVP * vec4(vertex_position, 1.0);\n"
"}\n";

static const char *fragmentMultiColourSource =
"#version 330\n"
"in vec4 frag_colour;\n"
"out vec4 colour;\n"
"void main() {\n"
"  colour = frag_colour;\n"
//"  colour = vec4(1.0, 1.0, 1.0, 1.0);"
"}\n";


GLuint compileMultiColourShaderProgram();

#endif /* MULTI_COLOUR_GL_SHADER_H */