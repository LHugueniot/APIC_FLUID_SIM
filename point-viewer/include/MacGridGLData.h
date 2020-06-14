#ifndef MAC_GRID_GL_DATA_H
#define MAC_GRID_GL_DATA_H

#include "MacGrid.h"
#include "Utilities.h"
#include "GLError.h"

//=====================================GRID====================================================

struct MacGridGLData{

	MacGridGLData(){}
	MacGridGLData(pic::MacGrid * _grid,
        GLuint * _monoColourShader,
        GLuint * _pressureShader,
        double _vecScale = 1.f,
        Vector3f _baseColour = {1.f, 1.f, 1.f}, 
        Vector3f _cellFaceVelColour = {1.f, 0.f, 0.f}) :
    grid(_grid),
    monoColourShader(_monoColourShader),
    pressureShader(_pressureShader),
	baseColour(_baseColour),
    cellFaceVelColour(_cellFaceVelColour),
	vecScale(_vecScale)
    {
        std::cout<<baseColour<<std::endl;
        std::cout<<cellFaceVelColour<<std::endl;
        //throw 1;
	}

    pic::MacGrid * grid;

    GLuint * monoColourShader;
    GLuint * pressureShader;
    double maxPressure = 10;

    Vector3f baseColour;

	double vecScale;

    //bool showCellVels = true;

    bool showCellFaceVels = false;
    GLuint cellFaceVelVerticesSize = 0;
    GLuint cellFaceVelVerticesBufferObject = 0;
    GLuint cellFaceVelVerticesArrayObject = 0;
    std::vector<double> cellFaceVelVertices;

    Vector3f cellFaceVelColour;

    bool showCellEdges = false;
    GLuint cellEdgesVerticesSize = 0;
    GLuint cellEdgesVerticesBufferObject = 0;
    GLuint cellEdgesVerticesArrayObject = 0;
    std::vector<double> cellEdgeVertices;

    bool showBoundary = true;
    GLuint boundEdgesVerticesSize = 0;
    GLuint boundEdgesVerticesBufferObject = 0;
    GLuint boundEdgesVerticesArrayObject = 0;
    std::vector<double> boundEdgesVertices;

    bool showPressures = false;
    GLuint pressuresVerticesSize = 0;
    GLuint pressuresVerticesBufferObject = 0;
    GLuint pressuresVerticesArrayObject = 0;
    std::vector<double> pressuresVertices;

    GLuint pressuresValsSize = 0;
    GLuint pressuresValsBufferObject = 0;
    GLuint pressuresValsArrayObject = 0;
    std::vector<double> pressuresVals;
};

//=====================================CELL FACE VEL========================================

void initCellFaceVelVertices(MacGridGLData & glData);
void initCellFaceVelVAO(MacGridGLData & glData);
void updateCellFaceVelVAO(MacGridGLData & glData);

//=====================================CELL VERTICES========================================

void initCellEdgesVAO(MacGridGLData & glData);
void updateCellEdgesVAO(MacGridGLData & glData);

//=====================================BOUND EDGE VERTICES==================================

void initBoundEdgesVAO(MacGridGLData & glData);
void updateBoundEdgesVAO(MacGridGLData & glData);

//=====================================PRESSURE VERTICES====================================

void updatePressureColours(MacGridGLData & glData);
void initPressureAttrs(MacGridGLData & glData);
void updatePressureVerticesVAO(MacGridGLData & glData);
void updatePressureColoursVAO(MacGridGLData & glData);
void initPressureVAO(MacGridGLData & glData);

//=====================================CORE FUNCS===========================================

void initMacGridVAO(MacGridGLData & glData);
void updateMacGridVAO(MacGridGLData const & glData);
void drawMacGrid(MacGridGLData const & glData, Eigen::Matrix4f & cameraMat);

#endif /* MAC_GRID_GL_DATA_H */