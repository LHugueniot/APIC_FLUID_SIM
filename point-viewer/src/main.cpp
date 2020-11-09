#include <SDL.h>
#include <PicCore.h>

#include "MonoColourGLShader.h"
#include "MultiColourGLShader.h"
#include "PressureGLShader.h"
#include "PointsGLData.h"
#include "PlaneGLData.h"
#include "MacGridGLData.h"
#include "Camera.h"
#include "GLError.h"

static std::vector<float> particleData = {
    -1.0f, -1.0f, 0.0f,
     1.0f, -1.0f, 0.0f,
     0.0f,  1.0f, 0.0f
};

using namespace pic;
static uint windowWidth = 1920;
static uint windowHeight = 1080;
//static const double TO_RADS = 3.141592654 / 180.0;

struct SDL_state
{
    SDL_Window* window = nullptr;
    SDL_GLContext context = nullptr;
};

void teardown(SDL_state state)
{
    SDL_GL_DeleteContext(state.context);
    SDL_DestroyWindow(state.window);
    SDL_Quit();
}

SDL_state setupSDL()
{
    SDL_state state;

    if (SDL_Init(SDL_INIT_VIDEO)<0 ){
        std::cout<<"Failed to init SDL: "<<SDL_GetError()<<std::endl;
        return state;
    }
    std::atexit(SDL_Quit);

    //SDL_GL_LoadLibrary(0);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, 
               SDL_GL_CONTEXT_PROFILE_CORE);

    //SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8);
    //SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8);
    //SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8);
    //SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    //SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);

    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 2);

    state.window = SDL_CreateWindow("Point Viewer", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                    windowWidth, windowHeight, SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN);
    if(!state.window){
        std::cout<<"Failed to create window: "<<SDL_GetError()<<std::endl;
        return state;
    }

    state.context = SDL_GL_CreateContext(state.window);
    if(!state.context){
        std::cout<<"Failed to create GL context: "<<SDL_GetError()<<std::endl;
        return state;
    }
    return state;
    //SDL_GL_Create_
}

int main(){

    //=====================================OpenGL/SDL SETUP=====================================
    auto state = setupSDL();

    if (GLEW_OK != glewInit())
    {
        // GLEW failed!
        exit(1);
    }

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    glClearColor(0.5f, 0.5f, 0.5f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    SDL_GL_SwapWindow(state.window);

    //Setup Camera
    auto camera = Camera(windowWidth, windowHeight);

    rotateCamera(camera, TO_RAD(-45.f));
    pitchCamera(camera,  TO_RAD(-45.f));

    updateCamera(camera);

    GLuint monoColourShader = compileMonoColourShaderProgram();
    if (monoColourShader == 0){
        std::cout<<"Mono Colour Shader failed to compile."<<std::endl;
        return 1;
    }

    GLuint multiColourShader = compileMultiColourShaderProgram();
    if (multiColourShader == 0){
        std::cout<<"Multi Colour Shader failed to compile."<<std::endl;
        return 1;
    }

    GLuint pressureShader = compilePressureShaderProgram();
    if (pressureShader == 0){
        std::cout<<"Pressure Shader failed to compile."<<std::endl;
        return 1;
    }

    //=====================================PICSIM SETUP=========================================

    double scale = 1.f/3.f;

    //Setup pic sim
    pic::MacGrid grid(Vector3d(0,0,0), 20, 20, 20, scale);
    //pic::FlipMacGrid grid(Vector3d(0,0,0), 20, 20, 20, scale);
    setDefaultCellStates(grid);
    //initializeLaplacianNBRMat(grid);
    //
    std::vector<double> particlePositions = pic::AABCubeUniformParticles(
        Vector3d(2, 4, 4) * scale, Vector3d(14, 16, 16) * scale, .5f * scale);
    //std::vector<double> particlePositions = pic::AABCubeUniformParticles(
    //    Vector3d(4, 4, 4), Vector3d(16, 16, 16), .5f);

    //std::vector<double> particlePositions = pic::AABCubeUniformParticles(Vector3d(2.5 , 7 ,2.5 ), 2, 64);
    //std::vector<double> particlePositions{2, 8 ,2};
    std::vector<double> particleVelocities(particlePositions.size(), 0);

    //for(uint i = 0 ; i < particleVelocities.size()/3 ; i++)
    //    particleVelocities[i * 3] = -9.8f;
    //pic::Particles particles(1.f, particlePositions, particleVelocities);
    pic::AffineParticles particles(1.f, particlePositions, particleVelocities);
    std::cout<<"particles.num: "<<particles.num<<std::endl;
    //throw 1;

    //=====================================OPENGL DATA SETUP====================================

    MacGridGLData MGdata(&grid, &monoColourShader, &pressureShader, 0.1f);
    initMacGridVAO(MGdata);

    double timeStep = 1.f/24.f;

    //Setup point viewer
    PointGLData points(&particles.positions, &monoColourShader);
    initPointsVAO(points);

    //Create center of world grid plain
    std::vector<float> gridPlainVertexData;
    //generateTile(gridPlainVertexData);
    //generateLine(gridPlainVertexData);
    generatePlaneVertexData(gridPlainVertexData, 1, 6, 6);
    PlaneGLData gridPlain(&gridPlainVertexData, &monoColourShader);
    initPlaneVAO(gridPlain);

    //=====================================MAIN LOOP============================================

    int step = 0;

    //initPointsVAO(point);
    bool quit = false;
    bool doStep = false;
    while(!quit)
    {
        check_gl_error();
        SDL_Event event;
        while (SDL_PollEvent(&event) != 0){
            switch (event.type){
                case SDL_QUIT:
                    quit = true;
                    break;
                
                case SDL_KEYDOWN:
                    switch(event.key.keysym.sym){
                        case SDLK_LEFT:
                            moveCamera(camera, Camera::ORBIT_LEFT);
                            break;
                        case SDLK_RIGHT:
                            moveCamera(camera, Camera::ORBIT_RIGHT);
                            break;
                        case SDLK_UP:
                            moveCamera(camera, Camera::ORBIT_UP);
                            break;
                        case SDLK_DOWN:
                            moveCamera(camera, Camera::ORBIT_DOWN);
                            break;
                        case SDLK_SPACE:
                            //doStep = !doStep;
                            pic::advanceStep(particles, grid, timeStep);
                            updatePointsVAO(points);
                            step++;
                            std::cout<<"step: "<<step<<std::endl;
                            break;
                        case SDLK_f:
                            MGdata.showCellFaceVels = !MGdata.showCellFaceVels;
                            break;
                        case SDLK_a:
                            MGdata.showCellEdges = !MGdata.showCellEdges;
                            break;
                        case SDLK_s:
                            MGdata.showBoundary = !MGdata.showBoundary;
                            break;
                        case SDLK_d:
                            MGdata.showPressures = !MGdata.showPressures;
                            break;
                    }
                    break;
                case SDL_MOUSEWHEEL:
                    if(event.wheel.y < 0){ // scroll up
                        moveCamera(camera, Camera::ZOOM_IN);
                    }
                    else if(event.wheel.y > 0){ // scroll down
                        moveCamera(camera, Camera::ZOOM_OUT);
                    }
                    //std::cout<< "Key press detected"<<std::endl;
                    break;
        
                case SDL_KEYUP:
                    //std::cout<< "Key release detected"<<std::endl;
                    break;

            }
            if (event.type == SDL_QUIT){
                quit = true;
            }
        }

        if(doStep){
            pic::advanceStep(particles, grid, timeStep);
            updatePointsVAO(points);
        }

        glClearColor(0.5f, 0.5f, 0.5f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        updateCamera(camera);

        //randomiseParticleStep(particleSystem);
        //updateParticleSystemVAO(particleSystem);
        Matrix4f cameraVP = camera.projMat * camera.viewMat;

        drawPoints(points, cameraVP);

        updatePlaneVAO(gridPlain);
        drawPlane(gridPlain, cameraVP);

        updateCellFaceVelVAO(MGdata);
        updatePressureColoursVAO(MGdata);
        drawMacGrid(MGdata, cameraVP);

        SDL_GL_SwapWindow(state.window);
    }

    teardown(state);
}
