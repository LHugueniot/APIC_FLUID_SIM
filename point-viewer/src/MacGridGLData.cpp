#include "MacGridGLData.h"

//=====================================CELL FACE VEL========================================

void initCellFaceVelVertices(MacGridGLData & glData){

    double vecScale = glData.vecScale;

    glData.cellFaceVelVertices.resize(
        glData.grid->cellFaceVel_u.size() * 3 * 2 + 
        glData.grid->cellFaceVel_v.size() * 3 * 2 +
        glData.grid->cellFaceVel_w.size() * 3 * 2);
    glData.cellFaceVelVerticesSize = glData.cellFaceVelVertices.size();

    int cellFaceVel_u_size = glData.grid->cellFaceVel_u.size();
    int cellFaceVel_v_size = glData.grid->cellFaceVel_v.size();
    int cellFaceVel_w_size = glData.grid->cellFaceVel_w.size();

    for (uint i = 0 ; i < glData.grid->cellNum_i ; i ++)
        for (uint j = 0 ; j < glData.grid->cellNum_j ; j++)
            for (uint k = 0 ; k < glData.grid->cellNum_k ; k ++){

                if(glData.grid->isValidFace<pic::minU>(i, j, k)){
                    auto pos = glData.grid->cellFacePos<pic::minU>(i, j, k);
                    double vel_u = glData.grid->cellFaceVel<pic::minU>(i, j, k);

                    auto g_idx_u = glData.grid->cellFaceIdx<pic::minU>(i, j, k);
                    int base_idx = g_idx_u * 6;

                    glData.cellFaceVelVertices[base_idx] = pos[0];
                    glData.cellFaceVelVertices[base_idx + 1] = pos[1];
                    glData.cellFaceVelVertices[base_idx + 2] = pos[2];

                    glData.cellFaceVelVertices[base_idx + 3] = vel_u * vecScale + pos[0];
                    glData.cellFaceVelVertices[base_idx + 4] = pos[1];
                    glData.cellFaceVelVertices[base_idx + 5] = pos[2];
                }
                if(glData.grid->isValidFace<pic::minV>(i, j, k)){
                    auto pos = glData.grid->cellFacePos<pic::minV>(i, j, k);
                    double vel_v = glData.grid->cellFaceVel<pic::minV>(i, j, k);

                    auto g_idx_v = glData.grid->cellFaceIdx<pic::minV>(i, j, k);
                    int base_idx = (g_idx_v + cellFaceVel_u_size) * 6;

                    glData.cellFaceVelVertices[base_idx] = pos[0];
                    glData.cellFaceVelVertices[base_idx + 1] = pos[1]; 
                    glData.cellFaceVelVertices[base_idx + 2] = pos[2];

                    glData.cellFaceVelVertices[base_idx + 3] = pos[0];
                    glData.cellFaceVelVertices[base_idx + 4] = vel_v * vecScale + pos[1];
                    glData.cellFaceVelVertices[base_idx + 5] = pos[2];
                }
                if(glData.grid->isValidFace<pic::minW>(i, j, k)){
                    auto pos = glData.grid->cellFacePos<pic::minW>(i, j, k);
                    double vel_w = glData.grid->cellFaceVel<pic::minW>(i, j, k);

                    auto g_idx_w = glData.grid->cellFaceIdx<pic::minW>(i, j, k);
                    int base_idx = (g_idx_w + cellFaceVel_u_size + cellFaceVel_v_size) * 6;

                    glData.cellFaceVelVertices[base_idx] = pos[0];
                    glData.cellFaceVelVertices[base_idx + 1] = pos[1];
                    glData.cellFaceVelVertices[base_idx + 2] = pos[2];

                    glData.cellFaceVelVertices[base_idx + 3] = pos[0];
                    glData.cellFaceVelVertices[base_idx + 4] = pos[1];
                    glData.cellFaceVelVertices[base_idx + 5] = vel_w * vecScale + pos[2];
                }
            }
}

void updateCellFaceVelVertices(MacGridGLData & glData){

    double vecScale = glData.vecScale;

    int cellFaceVel_u_size = glData.grid->cellFaceVel_u.size();
    int cellFaceVel_v_size = glData.grid->cellFaceVel_v.size();
    int cellFaceVel_w_size = glData.grid->cellFaceVel_w.size();

    for (uint i = 0 ; i < glData.grid->cellNum_i ; i ++)
        for (uint j = 0 ; j < glData.grid->cellNum_j ; j++)
            for (uint k = 0 ; k < glData.grid->cellNum_k ; k ++){

                if(glData.grid->isValidFace<pic::minU>(i, j, k)){
                    double vel_u = glData.grid->cellFaceVel<pic::minU>(i, j, k);
                    auto g_idx_u = glData.grid->cellFaceIdx<pic::minU>(i, j, k);
                    int base_idx = g_idx_u * 6;

                    glData.cellFaceVelVertices[base_idx + 3] = vel_u * vecScale + glData.cellFaceVelVertices[base_idx];
                }
                if(glData.grid->isValidFace<pic::minV>(i, j, k)){
                    double vel_v = glData.grid->cellFaceVel<pic::minV>(i, j, k);
                    auto g_idx_v = glData.grid->cellFaceIdx<pic::minV>(i, j, k);
                    int base_idx = (g_idx_v + cellFaceVel_u_size) * 6;

                    glData.cellFaceVelVertices[base_idx + 4] = vel_v * vecScale + glData.cellFaceVelVertices[base_idx + 1];
                }
                if(glData.grid->isValidFace<pic::minW>(i, j, k)){
                    double vel_w = glData.grid->cellFaceVel<pic::minW>(i, j, k);
                    auto g_idx_w = glData.grid->cellFaceIdx<pic::minW>(i, j, k);
                    int base_idx = (g_idx_w + cellFaceVel_u_size + cellFaceVel_v_size) * 6;

                    glData.cellFaceVelVertices[base_idx + 5] = vel_w * vecScale + glData.cellFaceVelVertices[base_idx + 2];
                }
            }
}

void initCellFaceVelVAO(MacGridGLData & glData){

    glGenBuffers(1, &glData.cellFaceVelVerticesBufferObject);

    initCellFaceVelVertices(glData);

    updateCellFaceVelVAO(glData);


    glGenVertexArrays(1, &glData.cellFaceVelVerticesArrayObject);
    glBindVertexArray(glData.cellFaceVelVerticesArrayObject);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, glData.cellFaceVelVerticesBufferObject);
    glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0, NULL);
    glDisableVertexAttribArray(0);
    DEBUG_VAR("initCellFaceVelVAO successfull!");
}

void updateCellFaceVelVAO(MacGridGLData & glData){

    updateCellFaceVelVertices(glData);

    glBindBuffer(GL_ARRAY_BUFFER, glData.cellFaceVelVerticesBufferObject);
    glBufferData(GL_ARRAY_BUFFER,
        sizeof(GLdouble) * glData.cellFaceVelVerticesSize,
        &glData.cellFaceVelVertices.data()[0], GL_STATIC_DRAW);
}

static const double lineCubeVBData[] = {
    //Back face edges
    0.f, 0.f, 0.f,
    1.f, 0.f, 0.f, 
    1.f, 0.f, 0.f,
    1.f, 1.f, 0.f, 
    1.f, 1.f, 0.f,
    0.f, 1.f, 0.f, 
    0.f, 1.f, 0.f,
    0.f, 0.f, 0.f, 

    //Middle edges
    0.f, 0.f, 0.f,
    0.f, 0.f, 1.f, 
    0.f, 1.f, 0.f,
    0.f, 1.f, 1.f, 
    1.f, 0.f, 0.f,
    1.f, 0.f, 1.f, 
    1.f, 1.f, 0.f,
    1.f, 1.f, 1.f, 

    //Front edges
    0.f, 0.f, 1.f,
    1.f, 0.f, 1.f, 
    1.f, 0.f, 1.f,
    1.f, 1.f, 1.f, 
    1.f, 1.f, 1.f,
    0.f, 1.f, 1.f, 
    0.f, 1.f, 1.f,
    0.f, 0.f, 1.f
};
//=====================================CELL VERTICES========================================

void initCellEdgeVertices(MacGridGLData & glData){

    int lineNum_i = glData.grid->cellNum_i + 1;
    int lineNum_j = glData.grid->cellNum_j + 1;
    int lineNum_k = glData.grid->cellNum_k + 1;

    glData.cellEdgeVertices.resize(
        lineNum_i * lineNum_j * 6 + 
        lineNum_i * lineNum_k * 6 +
        lineNum_j * lineNum_k * 6);
    glData.cellEdgesVerticesSize = glData.cellEdgeVertices.size();

    std::vector<bool> accessed(glData.cellEdgeVertices.size(), false);

    double lineStart_i = glData.grid->origin[0];
    double lineStart_j = glData.grid->origin[1];
    double lineStart_k = glData.grid->origin[2];

    double lineEnd_i = glData.grid->cellNum_i * glData.grid->cellSize + glData.grid->origin[0];
    double lineEnd_j = glData.grid->cellNum_j * glData.grid->cellSize + glData.grid->origin[1];
    double lineEnd_k = glData.grid->cellNum_k * glData.grid->cellSize + glData.grid->origin[2];

    int dimStride = 0;
    DEBUG_VAR(dimStride);
    for (uint i = 0 ; i < lineNum_i ; i ++)
        for (uint j = 0 ; j < lineNum_j ; j++){
                uint strideIdx = (i * lineNum_j + j )* 6;
                glData.cellEdgeVertices.at(strideIdx)     = i * glData.grid->cellSize + glData.grid->origin[0];
                glData.cellEdgeVertices.at(strideIdx + 1) = j * glData.grid->cellSize + glData.grid->origin[1];
                glData.cellEdgeVertices.at(strideIdx + 2) = lineStart_k;

                glData.cellEdgeVertices.at(strideIdx + 3) = glData.cellEdgeVertices[strideIdx];
                glData.cellEdgeVertices.at(strideIdx + 4) = glData.cellEdgeVertices[strideIdx + 1];
                glData.cellEdgeVertices.at(strideIdx + 5) = lineEnd_k;
//accessed.at(strideIdx)     = true;
//accessed.at(strideIdx + 1) = true;
//accessed.at(strideIdx + 2) = true;
//accessed.at(strideIdx + 3) = true;
//accessed.at(strideIdx + 4) = true;
//accessed.at(strideIdx + 5) = true;
            }

    dimStride += lineNum_i * lineNum_j * 6;
DEBUG_VAR(dimStride);
    for (uint i = 0 ; i < lineNum_i ; i ++)
        for (uint k = 0 ; k < lineNum_k ; k++){
                uint strideIdx = dimStride + (i * lineNum_k + k) * 6;
                glData.cellEdgeVertices.at(strideIdx)     = i * glData.grid->cellSize + glData.grid->origin[0];
                glData.cellEdgeVertices.at(strideIdx + 1) = lineStart_j;
                glData.cellEdgeVertices.at(strideIdx + 2) = k * glData.grid->cellSize + glData.grid->origin[2];

                glData.cellEdgeVertices.at(strideIdx + 3) = glData.cellEdgeVertices[strideIdx];
                glData.cellEdgeVertices.at(strideIdx + 4) = lineEnd_j;
                glData.cellEdgeVertices.at(strideIdx + 5) = glData.cellEdgeVertices[strideIdx + 2];
            }

    dimStride += lineNum_i * lineNum_k * 6;
DEBUG_VAR(dimStride);
    for (uint j = 0 ; j < lineNum_j ; j++)
        for (uint k = 0 ; k < lineNum_k ; k++){
            uint strideIdx = dimStride + (j * lineNum_k + k) * 6;
            glData.cellEdgeVertices.at(strideIdx)     = lineStart_i;
            glData.cellEdgeVertices.at(strideIdx + 1) = j * glData.grid->cellSize + glData.grid->origin[1];
            glData.cellEdgeVertices.at(strideIdx + 2) = k * glData.grid->cellSize + glData.grid->origin[2];

            glData.cellEdgeVertices.at(strideIdx + 3) = lineEnd_i;
            glData.cellEdgeVertices.at(strideIdx + 4) = glData.cellEdgeVertices[strideIdx + 1];
            glData.cellEdgeVertices.at(strideIdx + 5) = glData.cellEdgeVertices[strideIdx + 2];
        }
}

void initCellEdgesVAO(MacGridGLData & glData){

    glGenBuffers(1, &glData.cellEdgesVerticesBufferObject);

    initCellEdgeVertices(glData);

    updateCellEdgesVAO(glData);

    glGenVertexArrays(1, &glData.cellEdgesVerticesArrayObject);
    glBindVertexArray(glData.cellEdgesVerticesArrayObject);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, glData.cellEdgesVerticesBufferObject);
    glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0, NULL);
    glDisableVertexAttribArray(0);
    DEBUG_VAR("initCellEdgesVAO successfull!");
}
void updateCellEdgesVAO(MacGridGLData & glData){

    glBindBuffer(GL_ARRAY_BUFFER, glData.cellEdgesVerticesBufferObject);
    glBufferData(GL_ARRAY_BUFFER,
        sizeof(GLdouble) * glData.cellEdgesVerticesSize,
        glData.cellEdgeVertices.data(), GL_STATIC_DRAW);
}

//=====================================BOUND EDGE VERTICES==================================

void initBoundEdgesVertices(MacGridGLData & glData){

    glData.boundEdgesVertices.resize(72);
    glData.boundEdgesVerticesSize = glData.boundEdgesVertices.size();


    double lineStart_i = glData.grid->origin[0] + glData.grid->cellSize;
    double lineStart_j = glData.grid->origin[1] + glData.grid->cellSize;
    double lineStart_k = glData.grid->origin[2] + glData.grid->cellSize;

    double lineEnd_i = (glData.grid->cellNum_i - 1) * glData.grid->cellSize + glData.grid->origin[0];
    double lineEnd_j = (glData.grid->cellNum_j - 1) * glData.grid->cellSize + glData.grid->origin[1];
    double lineEnd_k = (glData.grid->cellNum_k - 1) * glData.grid->cellSize + glData.grid->origin[2];

    for(int i = 0; i < 72/3 ; i++)
    {
        glData.boundEdgesVertices[i * 3 + 0] = lineCubeVBData[i * 3 + 0] ? lineEnd_i : lineStart_i;
        glData.boundEdgesVertices[i * 3 + 1] = lineCubeVBData[i * 3 + 1] ? lineEnd_j : lineStart_j;
        glData.boundEdgesVertices[i * 3 + 2] = lineCubeVBData[i * 3 + 2] ? lineEnd_k : lineStart_k;
    }

}

void initBoundEdgesVAO(MacGridGLData & glData){

    glGenBuffers(1, &glData.boundEdgesVerticesBufferObject);

    initBoundEdgesVertices(glData);

    updateBoundEdgesVAO(glData);

    glGenVertexArrays(1, &glData.boundEdgesVerticesArrayObject);
    glBindVertexArray(glData.boundEdgesVerticesArrayObject);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, glData.boundEdgesVerticesBufferObject);
    glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0, NULL);
    glDisableVertexAttribArray(0);
    DEBUG_VAR("initBoundEdgesVAO successfull!");
}
void updateBoundEdgesVAO(MacGridGLData & glData){

    glBindBuffer(GL_ARRAY_BUFFER, glData.boundEdgesVerticesBufferObject);
    glBufferData(GL_ARRAY_BUFFER,
        sizeof(GLdouble) * glData.boundEdgesVerticesSize,
        glData.boundEdgesVertices.data(), GL_STATIC_DRAW);
}

//=====================================PRESSURE VERTICES====================================

void updatePressureColours(MacGridGLData & glData){

    int gCellNum_i = glData.grid->cellNum_i;
    int gCellNum_j = glData.grid->cellNum_j;
    int gCellNum_k = glData.grid->cellNum_k;
    for (uint i = 0 ; i < gCellNum_i ; i ++)
        for (uint j = 0 ; j < gCellNum_j ; j++)
            for (uint k = 0 ; k < gCellNum_k ; k ++){
                auto idx = glData.grid->cellCenterIdx(i,j,k);
                glData.pressuresVals.at(idx) = glData.grid->cellCenterPressure.at(idx);
                //glData.pressuresVals.at(idx) = glData.grid->cellCenterState[idx];
            }
}

void initPressureAttrs(MacGridGLData & glData){

    int gCellNum_i = glData.grid->cellNum_i;
    int gCellNum_j = glData.grid->cellNum_j;
    int gCellNum_k = glData.grid->cellNum_k;

    glData.pressuresVertices.resize(gCellNum_i * gCellNum_j * gCellNum_k * 3);
    glData.pressuresVerticesSize = glData.pressuresVertices.size();

    glData.pressuresVals.resize(gCellNum_i * gCellNum_j * gCellNum_k);
    glData.pressuresValsSize =  glData.pressuresVals.size();


    for (uint i = 0 ; i < gCellNum_i ; i ++)
        for (uint j = 0 ; j < gCellNum_j ; j++)
            for (uint k = 0 ; k < gCellNum_k ; k ++){
                auto idx = glData.grid->cellCenterIdx(i,j,k);
                auto pos = glData.grid->cellCenterPos(i,j,k);
                glData.pressuresVertices.at(idx * 3 + 0) = pos[0];
                glData.pressuresVertices.at(idx * 3 + 1) = pos[1];
                glData.pressuresVertices.at(idx * 3 + 2) = pos[2];
                glData.pressuresVals.at(idx) = glData.grid->cellCenterPressure.at(idx);
                //glData.pressuresVals.at(idx) = glData.grid->cellCenterState[idx];
            }
}

void updatePressureVerticesVAO(MacGridGLData & glData){

    glBindBuffer(GL_ARRAY_BUFFER, glData.pressuresVerticesBufferObject);
    glBufferData(GL_ARRAY_BUFFER,
        sizeof(GLdouble) * glData.pressuresVerticesSize,
        &glData.pressuresVertices.data()[0], GL_STATIC_DRAW);
}

void updatePressureColoursVAO(MacGridGLData & glData){

    updatePressureColours(glData);
    glBindBuffer(GL_ARRAY_BUFFER, glData.pressuresValsBufferObject);
    glBufferData(GL_ARRAY_BUFFER,
        sizeof(GLdouble) * glData.pressuresValsSize,
        &glData.pressuresVals.data()[0], GL_STATIC_DRAW);
}

void initPressureVAO(MacGridGLData & glData){

    initPressureAttrs(glData);

    glGenBuffers(1, &glData.pressuresVerticesBufferObject);
    updatePressureVerticesVAO(glData);

    glGenBuffers(1, &glData.pressuresValsBufferObject);
    updatePressureColoursVAO(glData);

    DEBUG_VAR("initPressureVAO successfull!");
}


//=====================================CORE FUNCS===========================================

void initMacGridVAO(MacGridGLData & glData){
    DEBUG();
    initCellFaceVelVAO(glData);
    initCellEdgesVAO(glData);
    initBoundEdgesVAO(glData);
    initPressureVAO(glData);
    DEBUG();
}

void drawMacGrid(MacGridGLData const & glData, Eigen::Matrix4f & cameraMat){

    GLuint monoColourMVPID = glGetUniformLocation(*glData.monoColourShader, "MVP");
    glUniformMatrix4fv(monoColourMVPID, 1, GL_FALSE, cameraMat.data());

    GLuint baseColID = glGetUniformLocation(*glData.monoColourShader, "base_colour");

    if(glData.showCellFaceVels){

        glUseProgram(*glData.monoColourShader);
        glUniform3f(baseColID, 
            glData.cellFaceVelColour[0],
            glData.cellFaceVelColour[1],
            glData.cellFaceVelColour[2]);

        glEnableVertexAttribArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, glData.cellFaceVelVerticesBufferObject);
        glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0, NULL);

        glDrawArrays(GL_LINES, 0, glData.cellFaceVelVerticesSize/3);

        glDisableVertexAttribArray(0);
    }
    if(glData.showCellEdges){

        glUseProgram(*glData.monoColourShader);
        glUniform3f(baseColID, 
            glData.baseColour[0],
            glData.baseColour[1],
            glData.baseColour[2]);

        glEnableVertexAttribArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, glData.cellEdgesVerticesBufferObject);
        glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0, NULL);

        glDrawArrays(GL_LINES, 0, glData.cellEdgesVerticesSize/3);

        glDisableVertexAttribArray(0);
    }
    if(glData.showBoundary){

        glUseProgram(*glData.monoColourShader);
        glUniform3f(baseColID, 
            glData.baseColour[0],
            glData.baseColour[1],
            glData.baseColour[2]);

        glEnableVertexAttribArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, glData.boundEdgesVerticesBufferObject);
        glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0, NULL);

        glDrawArrays(GL_LINES, 0, glData.boundEdgesVerticesSize/3);

        glDisableVertexAttribArray(0);
    }
    if(glData.showPressures){
    
        glUseProgram(*glData.pressureShader);

        GLuint multiColourMVPID = glGetUniformLocation(*glData.pressureShader, "MVP");
        glUniformMatrix4fv(multiColourMVPID, 1, GL_FALSE, cameraMat.data());

        GLuint max_pressure_ID = glGetUniformLocation(*glData.pressureShader, "max_pressure");
        glUniform1f(max_pressure_ID, glData.maxPressure);

        float defaultPointSize;
        glGetFloatv(GL_POINT_SIZE, &defaultPointSize);
        
        glPointSize(10);

        glEnableVertexAttribArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, glData.pressuresVerticesBufferObject);
        glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, 0, NULL);

        glEnableVertexAttribArray(1);
        glBindBuffer(GL_ARRAY_BUFFER, glData.pressuresValsBufferObject);
        glVertexAttribPointer(1, 1, GL_DOUBLE, GL_FALSE, 0, NULL);

        glDrawArrays(GL_POINTS, 0, glData.pressuresVerticesSize/3);

        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);

        glPointSize(defaultPointSize);
    }
}
