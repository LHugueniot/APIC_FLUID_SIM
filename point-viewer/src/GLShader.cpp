#include "GLShader.h"

GLuint compileShaderProgram(std::string const & vertexSource, std::string const & fragmentSource){

	// Create an empty vertex shader handle
	GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);

	const GLchar *vertexSource_c_str = (const GLchar *)vertexSource.c_str();
	// Send the vertex shader source code to GL
	glShaderSource(vertexShader, 1, &vertexSource_c_str, 0);

	// Compile the vertex shader
	glCompileShader(vertexShader);

	GLint isCompiled = 0;

	glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &isCompiled);

	if(isCompiled == GL_FALSE){
		std::cout<<"Vertex shader failed to compile."<<std::endl;
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
	
	const GLchar * fragmentSource_c_str = (const GLchar *)fragmentSource.c_str();

	// Send the fragment shader source code to GL
	glShaderSource(fragmentShader, 1, &fragmentSource_c_str, 0);
	
	// Compile the fragment shader
	glCompileShader(fragmentShader);
	
	glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &isCompiled);

	if (isCompiled == GL_FALSE){
		std::cout<<"Fragment shader failed to compile."<<std::endl;
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
	if (isLinked == GL_FALSE){
		std::cout<<"Shader failed to link."<<std::endl;
		GLint maxLength = 0;
		glGetProgramiv(program, GL_INFO_LOG_LENGTH, &maxLength);
	
		// The maxLength includes the NULL character
		char infoLog[maxLength];
		glGetProgramInfoLog(program, maxLength, &maxLength, &infoLog[0]);
		
		// We don't need the program anymore.
		glDeleteProgram(program);
		// Don't leak shaders either.
		glDeleteShader(vertexShader);
		glDeleteShader(fragmentShader);
	
		// Use the infoLog as you see fit.
		std::cout<<infoLog<<std::endl;
		
		// In this simple program, we'll just leave
		return 0;
	}
	// Always detach shaders after a successful link.
	glDetachShader(program, vertexShader);
	glDetachShader(program, fragmentShader);

	return program;
}