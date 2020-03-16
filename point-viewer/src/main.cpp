#include <SDL.h>
#include <picCore.h>

#include "basicShader.h"
#include "basicGeom.h"
#include "camera.h"

#define DEBUG() std::cout<<"BREAK POINT: LINE "<<__LINE__<<" IN "<<__FILE__<<std::endl

static std::vector<float> particleData = {
    -1.0f, -1.0f, 0.0f,
     1.0f, -1.0f, 0.0f,
     0.0f,  1.0f, 0.0f
};

namespace pic{

auto sphereCollision = 
[](double pPos_x, double pPos_y, double pPos_z,
   double pVel_x, double pVel_y, double pVel_z)
-> std::array<double,3>
{
    //=========================================Vector Funcs=========================================

    //==============================================================================================

    double sCenter_x = 5;
    double sCenter_y = 5;
    double sCenter_z = 5;

    double sphereRadius = 1;

    double t;

    std::array<double,3> newPos = add({pPos_x, pPos_y, pPos_z}, {pVel_x, pVel_y, pVel_z});

    auto intersected = intersectRaySphere({pPos_x, pPos_y, pPos_z}, {pVel_x, pVel_y, pVel_z},
        {sCenter_x, sCenter_y, sCenter_z}, sphereRadius, t, newPos);

    return newPos;
};

}

static uint windowWidth = 640;
static uint windowHeight = 480;
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

int main()
{
    auto state = setupSDL();

    if (GLEW_OK != glewInit())
    {
        // GLEW failed!
        exit(1);
    }

    glClearColor(0.5f, 0.5f, 0.5f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    SDL_GL_SwapWindow(state.window);

    //Setup Camera
    auto camera = cameraData( windowWidth, windowHeight);
    rotateCamera(camera, glm::radians(-45.f));
    pitchCamera(camera,  glm::radians(-45.f));

    updateCamera(camera);

    GLuint my_basicShader = compileBasicShaderProgram();

    //Create center of world grid plain
    auto gridPlaine = gridData();
    std::vector<float> gridPlaineVertexData = generateGridVertexData(1, 5, 5);
    populateGrid(gridPlaine, gridPlaineVertexData);

    //Create particle system
    auto particleSystem = pic::ParticleAttributes();

    pic::randomisedParticleBB(particleSystem, 1, 5, 5, 5, 5, 5, 5);

    //Create Bounding box


    //Setup point drawer
    pointData points;

    populatePoints(points, 
        particleSystem.positions_x,
        particleSystem.positions_y,
        particleSystem.positions_z);

    //return 0;
    //Create cell grid

    pic::GridAttributes grid(10, 10, 10, 1);

    if (my_basicShader == 0)
    {
        std::cout<<"Shader failed to compile."<<std::endl;
        return 1;
    }

    bool quit = false;
    while(!quit)
    {
        SDL_Event event;
        if (SDL_PollEvent(&event) != 0)
        {
            switch (event.type)
            {
                case SDL_QUIT:
                    quit = true;
                    break;
                
                case SDL_KEYDOWN:
                    switch(event.key.keysym.sym)
                    {
                        case SDLK_LEFT:
                            moveCamera(camera, ORBIT_LEFT);
                            break;
                        case SDLK_RIGHT:
                            moveCamera(camera, ORBIT_RIGHT);
                            break;
                        case SDLK_UP:
                            moveCamera(camera, ORBIT_UP);
                            break;
                        case SDLK_DOWN:
                            moveCamera(camera, ORBIT_DOWN);
                            break;
                    }
                    
                    printf( "Key press detected\n" );
                    break;
        
                case SDL_KEYUP:
                    printf( "Key release detected\n" );
                    break;

            }
            if (event.type == SDL_QUIT)
            {
                quit = true;
            }
        }

        glClearColor(0.0f, 0.0f, 0.0f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        updateCamera(camera);

        //randomiseParticleStep(particleSystem);
        //updateParticleSystemVAO(particleSystem);
        glm::mat4 cameraVP = camera.projectionMat * camera.viewMat;

            //DEBUG();
        pic::transferAttributes(particleSystem, grid);
            //DEBUG();
        pic::transferAttributes(grid, particleSystem);
            //DEBUG();
        pic::timeStep(particleSystem, pic::sphereCollision, 0.001);
            //DEBUG();

        std::cout<<points.vertexData.size()<<std::endl;
        points.vertexData = flattenPointCoorAttr(
            particleSystem.positions_x,
            particleSystem.positions_y,
            particleSystem.positions_z);

        std::cout<<points.vertexData.size()<<std::endl;
            //DEBUG();
        updatePointsVAO(points);
            //DEBUG();
        drawPoints(points, my_basicShader, cameraVP);
            //DEBUG();
        
        drawGrid(gridPlaine, my_basicShader, cameraVP);
        
        SDL_GL_SwapWindow(state.window);
    }

    teardown(state);
}

