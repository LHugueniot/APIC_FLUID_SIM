#include <SDL.h>

#include "particleSystem.h"
#include "basicShader.h"
#include "basicGeom.h"
#include "camera.h"
#include "grid.h"

static std::vector<float> particleData = {
    -1.0f, -1.0f, 0.0f,
     1.0f, -1.0f, 0.0f,
     0.0f,  1.0f, 0.0f
};

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

    auto camera = cameraData( windowWidth, windowHeight);
    rotateCamera(camera, glm::radians(-45.f));
    pitchCamera(camera,  glm::radians(-45.f));

    updateCamera(camera);

    GLuint my_basicShader = compileBasicShaderProgram();

    auto particleSystem = particleSystemData();
    //std::vector<float> particleVertexData = generateRandomParticleData(0, 1, 0, 1, 1, 1, 10000);
    std::vector<float> particleVertexData = createCellGridVertexData(0, 0, 0, 10, 10, 10, 10, 10, 10);
    populateParticleSystem(particleSystem, particleVertexData);

    auto grid = gridData();
    std::vector<float> gridVertexData = generateGridVertexData(1, 5, 5);
    populateGrid(grid, gridVertexData);

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
        glm::mat4 cameraVP = camera.projectionMat * 
                              camera.viewMat;

        drawParticleSystem(particleSystem, my_basicShader, cameraVP);

        drawGrid(grid, my_basicShader, cameraVP);
        SDL_GL_SwapWindow(state.window);
    }

    teardown(state);
}

