#include <windows.h>	// for timeGetTime()
#include <mmsystem.h>	// ditto
#include <vector>
#include <iostream>		// I/O
#include "glew.h"
#include "glut.h"            // Glut
#include "textfile.h"
#include "Camera.h"		// Camera
#include "MathAux.h"	
#include "OBJLoader.h"
#include "Model.h"
#include "GLToMovie.h"
#include "Solver.h"

/******************************************
********* defines *************************
******************************************/
#define WINDOW_DEFAULT_WIDTH 700
#define WINDOW_DEFAULT_HEIGHT 700

#define VIPER_OBJ_FILE "../res/obj/viper.obj"

#define VIPER_TEXTURE_FILE "../res/obj/viper.png"

#define FLUIDSVS_FILE "../shaders/fluidDensity.vert"

#define FLUIDSFS_FILE "../shaders/fluidDensity.frag"

#define RECORDMOVIE false

#define WINDOW_TITLE "RTP Stable fluids"


#define FLUID_SIDE_DIM 64

#define FLUID3D_EDGE_DIM 64



enum Demo {
	FluidDensity
};

Demo currentDemo = FluidDensity;

bool simulate3D = true;
bool simulate2D = false;
bool simulate2DIn3D = false;

RTP::FluidSolver fluidSolver;

RTP::FluidSolver3D tridimFluidSolver;

/* Simulation switches */
bool drawFluidContainer = false;
bool drawFluidDensity = true;
bool drawFluidVelocity = false;
bool drawMarkerParticleSet = false;

bool lerpMarkerParticleVelocities = true;

/* Source input */
float densityX, densityY, uX, uY, vX, vY, uDelta, vDelta;

void draw2DFluidSimulation ();
// Should provide a plane to align the sim
void draw2DFluidSimulationIn3D ();

void draw3DFluidSimulation ();

void init2DFluidSolver ();
void init3DFluidSolver ();


/////////////// TEST ////////////////////
int win_x = WINDOW_DEFAULT_WIDTH;
int win_y = WINDOW_DEFAULT_HEIGHT;
int omx, omy, mx, my;
bool mouse_down[3];

/* Window and scene-related variables */
bool wireframe = false;
int windowId;
DWORD lastTickCount;
bool computeAndShowFPS = false;
/* FPS value */
char fPSStr[15] = {0};

/******************************************
********* FUNCTION DECLARATIONS ***********
******************************************/
/* Frame rate computation */
void calculateFrameRate (DWORD currentTime);

/* Scene-related function declarations */
void setShaders (char* vertexShaderFn, char* fragmentShaderFn);
void setupScene();
void updateScene();
void renderScene();
void setViewport (int w, int h);
void exitScene();

/* Shader loading related variables */
GLuint v,f,p;

void restartDemo();

/* Callback functions to process keyboard input */
void keypress (unsigned char key, int x, int y);
void specialkeypress (int key, int xx, int yy);
void keyup (unsigned char key, int x, int y);
void specialkeyup (int key, int xx, int yy);

/* We buffer keys, so to be able to check multiple key presses */
bool keys[256];
bool specialkeys[248];
void processKeyboard();
/* Callback functions to process mouse input */
void Mouse (int b, int s, int x, int y);
void Motion (int x, int y);

/* Simple on screen information */
void drawFPS ();
void drawStats ();
void drawText ();

//Mouse input state variables
int lastMouseX = 0;
int lastMouseY = 0;
unsigned char mouseButtons[3] = {0};

/* Demo-related variables */
bool showHelp = false;
//Camera
//RTA::Camera camera(130.0f, 0.0f, 10.0f, -10.0f, -1.0f);
RTA::Camera camera(140.0f, 0.0f, 10.0f, -10.0f, -1.0f);


/* Some control flags to decide what to display and animate */
bool pause = false;


/* Cubemap texture */
GLenum  cubemap[6] = {  GL_TEXTURE_CUBE_MAP_POSITIVE_X,
						GL_TEXTURE_CUBE_MAP_NEGATIVE_X,
						GL_TEXTURE_CUBE_MAP_POSITIVE_Y,
						GL_TEXTURE_CUBE_MAP_NEGATIVE_Y,
						GL_TEXTURE_CUBE_MAP_POSITIVE_Z,
						GL_TEXTURE_CUBE_MAP_NEGATIVE_Z };
char* cubemapTextures[6] = {"../res/textures/posx.png",
							"../res/textures/negx.png",
							"../res/textures/posy.png",
							"../res/textures/negy.png",
							"../res/textures/posz.png",
							"../res/textures/negz.png"};
float cubeMapExtent = 115.0f;

bool loadCubeMap();

void drawCubeMap();

// LIGHT 0
GLfloat lpos[3] = {20.0f, 0.0f, 10.0f};


/////////////////////
// MOVIE RECORDING //
/////////////////////

bool recordMovie = RECORDMOVIE;

CGLToMovie g_MovieRecorder = CGLToMovie(L"D:\\Eliminare\\StableFluids.avi", WINDOW_DEFAULT_WIDTH, WINDOW_DEFAULT_HEIGHT, 24, mmioFOURCC('X','V','I','D'), 15);

