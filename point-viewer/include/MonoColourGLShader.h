#ifndef MONO_COLOUR_GL_SHADER_H
#define MONO_COLOUR_GL_SHADER_H

#include "GLShader.h"

static const char* vertexMonoColourSource =
"#version 330\n"
"layout(location = 0) in vec3 vertex_position;\n"
"uniform vec3 base_colour;\n"
"uniform mat4 MVP;\n"
"out vec3 frag_colour;\n"
"void main() {\n"
"	frag_colour = base_colour;\n"
"  	gl_Position = MVP * vec4(vertex_position, 1.0);\n"
"}\n";

static const char *fragmentMonoColourSource =
"#version 330\n"
"in vec3 frag_colour;\n"
"out vec4 colour;\n"
"void main() {\n"
"  colour = vec4(frag_colour, 1.0);\n"
"}\n";


GLuint compileMonoColourShaderProgram();

#endif /* MONO_COLOUR_GL_SHADER_H */