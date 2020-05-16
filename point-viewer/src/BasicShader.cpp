#include "BasicShader.h"

static const char* vertexSource =
"#version 330\n"
"in vec3 position;\n"
"uniform mat4 MVP;\n"
"void main() {\n"
"  gl_Position = MVP * vec4(position, 1.0);\n"
"}\n";

static const char *fragmentSource =
"#version 330\n"
"out vec4 frag_colour;\n"
"void main() {\n"
"  frag_colour = vec4(1.0, 1.0, 1.0, 1.0);\n"
"}\n";

GLuint compileBasicShaderProgram()
{
	// Create an empty vertex shader handle
	GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);

	// Send the vertex shader source code to GL
	glShaderSource(vertexShader, 1, &vertexSource, 0);

	// Compile the vertex shader
	glCompileShader(vertexShader);

	GLint isCompiled = 0;

	glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &isCompiled);

	if(isCompiled == GL_FALSE)
	{
		GLint maxLength = 0;
		glGetShaderiv(vertexShader, GL_INFO_LOG_LENGTH, &maxLength);

		// The maxLength includes the NULL character
		char infoLog[maxLength];
		glGetShaderInfoLog(vertexShader, maxLength, &maxLength, &infoLog[0]);

		// We don't need the shader anymore.
		glDeleteShader(vertexShader);

		// Use the infoLog as you see fit.
		std::cout<<infoLog<<std::endl;
		
		// In this simple program, we'll just leave
		return 0;
	}

	GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	
	// Send the fragment shader source code to GL
	glShaderSource(fragmentShader, 1, &fragmentSource, 0);
	
	// Compile the fragment shader
	glCompileShader(fragmentShader);
	
	glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &isCompiled);

	if (isCompiled == GL_FALSE)
	{
		GLint maxLength = 0;
		glGetShaderiv(fragmentShader, GL_INFO_LOG_LENGTH, &maxLength);

		// The maxLength includes the NULL character
		char infoLog[maxLength];
		//std::vector<GLchar> infoLog(maxLength);
		glGetShaderInfoLog(fragmentShader, maxLength, &maxLength, &infoLog[0]);

		// We don't need the shader anymore.
		glDeleteShader(fragmentShader);
		// Either of them. Don't leak shaders.
		glDeleteShader(vertexShader);

		// Use the infoLog as you see fit.
		std::cout<<infoLog<<std::endl;
		
		// In this simple program, we'll just leave
		return 0;
	}

	// Vertex and fragment shaders are successfully compiled.
	// Now time to link them together into a program.
	// Get a program object.
	GLuint program = glCreateProgram();
	
	// Attach our shaders to our program
	glAttachShader(program, vertexShader);
	glAttachShader(program, fragmentShader);
	
	// Link our program
	glLinkProgram(program);
	
	// Note the different functions here: glGetProgram* instead of glGetShader*.
	GLint isLinked = 0;
	glGetProgramiv(program, GL_LINK_STATUS, (int *)&isLinked);
	if (isLinked == GL_FALSE)
	{
		GLint maxLength = 0;
		glGetProgramiv(program, GL_INFO_LOG_LENGTH, &maxLength);
	
		// The maxLength includes the NULL character
		std::vector<GLchar> infoLog(maxLength);
		glGetProgramInfoLog(program, maxLength, &maxLength, &infoLog[0]);
		
		// We don't need the program anymore.
		glDeleteProgram(program);
		// Don't leak shaders either.
		glDeleteShader(vertexShader);
		glDeleteShader(fragmentShader);
	
		// Use the infoLog as you see fit.
		
		// In this simple program, we'll just leave
		return 0;
	}
	// Always detach shaders after a successful link.
	glDetachShader(program, vertexShader);
	glDetachShader(program, fragmentShader);

	return program;
}


//static const GLchar *vertexSource =
//	"#version 330\n"
//    "attribute highp vec4 posAttr;\n"
//    "attribute lowp vec4 colAttr;\n"
//    "varying lowp vec4 col;\n"
//    "uniform highp mat4 matrix;\n"
//    "void main() {\n"
//    "   col = colAttr;\n"
//    "   gl_Position = matrix * posAttr;\n"
//    "}\n";
//
//static const char *fragmentSource =
//	"#version 330\n"
//    "varying lowp vec4 col;\n"
//    "void main() {\n"
//    "   gl_FragColor = col;\n"
//    "}\n";
