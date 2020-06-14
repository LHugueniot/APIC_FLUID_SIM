#pragma once
#ifndef EXPORTER_H
#define EXPORTER_H

#include "Particles.h"

#include <iostream>
#include <fstream> 
#include <string>
#include <regex>

namespace pic{


static std::string const fileExtension = ".pic";

static std::regex 
POINT_REGEX("^p (\\d{1,6}\\.\\d{1,6}) (\\d{1,6}\\.\\d{1,6}) (\\d{1,6}\\.\\d{1,6})");
static std::regex 
FRAME_REGEX("^frame (\\d*)");
static std::regex 
FPS_REGEX("^fps (\\d+\\.\\d*){0,1}");

class ParticleCacheO{

public:
	ParticleCacheO(std::string _outputFilePath, uint _fps) :
    outputFilePath(_outputFilePath),
    fps(_fps){
    	assert(outputFilePath.length());
    	assert(Initialize());
    }


	~ParticleCacheO(){
		pointFile.close();
	}

	bool writePositions( uint frame, std::vector<Vector3d> const & pointPositions);

private:

	bool Initialize();

	std::string outputFilePath;
	uint totalFrames;
	uint fps;

	std::ofstream pointFile;
};

class ParticleCacheI{

public:
	ParticleCacheI(std::string _inputFilePath) :
    _inputFilePath(_inputFilePath){
    	assert(_inputFilePath.length());
    	assert(Initialize());
    }


	~ParticleCacheI(){
		for(int frame = 0 ; frame < totalFrames ; frame++)
			pointFile.close();
	}

	bool readPositions( uint frame, std::vector<Vector3d> & pointPositions);

private:

	bool Initialize();

	std::string _inputFilePath;
	uint totalFrames;
	uint fps;

	std::ifstream pointFile;
};

}


#endif /* EXPORTER_H */