#ifndef PRESSURE_GL_SHADER_H
#define PRESSURE_GL_SHADER_H

#include "GLShader.h"

static const char* vertexPressureSource =
"#version 330\n"
"uniform mat4 MVP;\n"
"layout(location = 0) in vec3 vertex_position;\n"
"layout(location = 1) in float pressure;\n"
"out vec4 frag_colour;\n"
"uniform float max_pressure;\n"
"void main(){\n"
" 	float new_pressure = 0.75/(-mod(pressure, max_pressure)/max_pressure-0.5) + 1.5;\n"
"	frag_colour = vec4(0.0, 0.0, 0.0, 0.5);\n"
"   if(new_pressure > 0.001){\n"
"		if (new_pressure > 0.5){\n"
"			new_pressure = new_pressure * 2.0 - 1.0;\n"
"			frag_colour[0] = (new_pressure) + 0.5;\n"
"			frag_colour[1] = (1.0f - new_pressure) + 0.5;\n"
"		}\n"
"		else{\n"
"			new_pressure *= 2.0;\n"
"			frag_colour[1] = (new_pressure) + 0.5;\n"
"			frag_colour[2] = (1.0 - new_pressure) + 0.5;\n"
"		}\n"
"  	gl_Position = MVP * vec4(vertex_position, 1.0);\n"
"	}else{\n"
"		gl_Position = vec4(-1, -1, -1, 1);\n"
"	}\n"
"}";

static const char *fragmentPressureSource =
"#version 330\n\n"
"in vec4 frag_colour;\n"
"out vec4 colour;\n"
"void main(){\n"
"  colour = frag_colour;\n"
"}";


GLuint compilePressureShaderProgram();

#endif /* PRESSURE_GL_SHADER_H */