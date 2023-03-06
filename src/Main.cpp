#include "Main.h"

void calculateFrameRate (DWORD currentTime) {
	/* Static variables, which maintain value between function calls */
	/* Used to create a delta for the measurement */
	static DWORD lastUpdateTime = 0;

	static float framesPerSecond;

	// Increase the frame counter
    ++framesPerSecond;

	/* Has 1s passed or not? */
    if (currentTime - lastUpdateTime > 1000) {
		lastUpdateTime = currentTime;

		/* Build up the result into the string */
		sprintf_s(fPSStr, "FPS: %d", int(framesPerSecond));

		/* Reset count */
        framesPerSecond = 0;
    }
}

void drawCubeMap () {

	/* Use the fixed pipeline */
	glUseProgram(0);

	 glBegin(GL_QUADS);
        //////////////////////////////////////////////
        // Negative X
        glTexCoord3f(-1.0f, -1.0f, 1.0f);
        glVertex3f(-cubeMapExtent, -cubeMapExtent, cubeMapExtent);
        
        glTexCoord3f(-1.0f, -1.0f, -1.0f);
        glVertex3f(-cubeMapExtent, -cubeMapExtent, -cubeMapExtent);
        
        glTexCoord3f(-1.0f, 1.0f, -1.0f);
        glVertex3f(-cubeMapExtent, cubeMapExtent, -cubeMapExtent);
        
        glTexCoord3f(-1.0f, 1.0f, 1.0f);
        glVertex3f(-cubeMapExtent, cubeMapExtent, cubeMapExtent);


        ///////////////////////////////////////////////
        //  Postive X
        glTexCoord3f(1.0f, -1.0f, -1.0f);
        glVertex3f(cubeMapExtent, -cubeMapExtent, -cubeMapExtent);
        
        glTexCoord3f(1.0f, -1.0f, 1.0f);
        glVertex3f(cubeMapExtent, -cubeMapExtent, cubeMapExtent);
        
        glTexCoord3f(1.0f, 1.0f, 1.0f);
        glVertex3f(cubeMapExtent, cubeMapExtent, cubeMapExtent);
        
        glTexCoord3f(1.0f, 1.0f, -1.0f);
        glVertex3f(cubeMapExtent, cubeMapExtent, -cubeMapExtent);
 

        ////////////////////////////////////////////////
        // Negative Z 
        glTexCoord3f(-1.0f, -1.0f, -1.0f);
        glVertex3f(-cubeMapExtent, -cubeMapExtent, -cubeMapExtent);
        
        glTexCoord3f(1.0f, -1.0f, -1.0f);
        glVertex3f(cubeMapExtent, -cubeMapExtent, -cubeMapExtent);
        
        glTexCoord3f(1.0f, 1.0f, -1.0f);
        glVertex3f(cubeMapExtent, cubeMapExtent, -cubeMapExtent);
        
        glTexCoord3f(-1.0f, 1.0f, -1.0f);
        glVertex3f(-cubeMapExtent, cubeMapExtent, -cubeMapExtent);


        ////////////////////////////////////////////////
        // Positive Z 
        glTexCoord3f(1.0f, -1.0f, 1.0f);
        glVertex3f(cubeMapExtent, -cubeMapExtent, cubeMapExtent);
        
        glTexCoord3f(-1.0f, -1.0f, 1.0f);
        glVertex3f(-cubeMapExtent, -cubeMapExtent, cubeMapExtent);
        
        glTexCoord3f(-1.0f, 1.0f, 1.0f);
        glVertex3f(-cubeMapExtent, cubeMapExtent, cubeMapExtent);
        
        glTexCoord3f(1.0f, 1.0f, 1.0f);
        glVertex3f(cubeMapExtent, cubeMapExtent, cubeMapExtent);


        //////////////////////////////////////////////////
        // Positive Y
        glTexCoord3f(-1.0f, 1.0f, 1.0f);
        glVertex3f(-cubeMapExtent, cubeMapExtent, cubeMapExtent);
        
        glTexCoord3f(-1.0f, 1.0f, -1.0f);
        glVertex3f(-cubeMapExtent, cubeMapExtent, -cubeMapExtent);
        
        glTexCoord3f(1.0f, 1.0f, -1.0f);
        glVertex3f(cubeMapExtent, cubeMapExtent, -cubeMapExtent);
        
        glTexCoord3f(1.0f, 1.0f, 1.0f);
        glVertex3f(cubeMapExtent, cubeMapExtent, cubeMapExtent);
  
    
        ///////////////////////////////////////////////////
        // Negative Y
        glTexCoord3f(-1.0f, -1.0f, -1.0f);
        glVertex3f(-cubeMapExtent, -cubeMapExtent, -cubeMapExtent);
        
        glTexCoord3f(-1.0f, -1.0f, 1.0f);
        glVertex3f(-cubeMapExtent, -cubeMapExtent, cubeMapExtent);
        
        glTexCoord3f(1.0f, -1.0f, 1.0f);
        glVertex3f(cubeMapExtent, -cubeMapExtent, cubeMapExtent);
        
        glTexCoord3f(1.0f, -1.0f, -1.0f);
        glVertex3f(cubeMapExtent, -cubeMapExtent, -cubeMapExtent);
    glEnd();
}

void renderScene(){
	if (!simulate2D) {
		glClearColor(0.0,0.0f,0.0f,1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glMatrixMode(GL_MODELVIEW);
		glEnable(GL_SMOOTH);
		glLoadIdentity();
		glPushMatrix();

		// Use camera directly
		glTranslatef(camera.getXTranslation(), camera.getYTranslation(), -camera.getZoom());
		glRotatef(camera.getXRotation(),1,0,0);
		glRotatef(camera.getYRotation(),0,1,0);
	}
	

	if (simulate3D) {
		//glEnable(GL_CULL_FACE);
		///* Cubemap */
		//// Enable cube mapping, and set texture environment to mudulate
		//glEnable(GL_TEXTURE_CUBE_MAP);
		//glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	
		//// Cubemap is manually textured
		//glDisable(GL_TEXTURE_GEN_S);
		//glDisable(GL_TEXTURE_GEN_T);
		//glDisable(GL_TEXTURE_GEN_R);        
		//// Perform drawing..
		//drawCubeMap();
		//
		//
		//////Use texgen to apply cube map
		////glEnable(GL_TEXTURE_GEN_S);
		////glEnable(GL_TEXTURE_GEN_T);
		////glEnable(GL_TEXTURE_GEN_R);


		//glDisable(GL_TEXTURE_CUBE_MAP);
		//glDisable(GL_CULL_FACE);
		//glUseProgram(p);
		draw3DFluidSimulation();
		//glUseProgram(0);
	}

	if (simulate2DIn3D) {
		draw2DFluidSimulationIn3D();
	}

	if (simulate2D) {
		draw2DFluidSimulation();
	}

	if (computeAndShowFPS)
		drawText();
	
	// Capture the Rendering into CGLToMovie's movie file
	if (recordMovie)
		g_MovieRecorder.RecordFrame();

	glutSwapBuffers();
        
}

void draw2DFluidSimulation () {
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	gluOrtho2D (0.0, 1.0, 0.0, 1.0);
	glClearColor (0.0f, 0.0f, 0.0f, 1.0f);
	glClear (GL_COLOR_BUFFER_BIT);

	if (drawFluidDensity) {
		fluidSolver.drawDensity();
	}

	if (drawFluidContainer) {
		fluidSolver.drawGrid();
	}

	if (drawFluidVelocity) {
		fluidSolver.drawVelocity();
	}

	
	if (drawMarkerParticleSet) {
		fluidSolver.drawMarkerParticleSet();
	}
}

void draw2DFluidSimulationIn3D () {
	if (drawFluidDensity) {
		fluidSolver.draw2DDensityIn3D();
	}
	if (drawFluidVelocity) {
		fluidSolver.draw2DVelocityIn3D();
	}
}

void draw3DFluidSimulation () {

	if (drawFluidVelocity) {
		tridimFluidSolver.drawVelocity();
	}

	glDisable(GL_DEPTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	if (drawFluidDensity) {
		tridimFluidSolver.drawDensityWithBoxes();
	}

	glEnable(GL_DEPTH);
	glDisable(GL_BLEND);

	if (drawMarkerParticleSet) {
		tridimFluidSolver.drawMarkerParticleSet();
	}

	if (drawFluidContainer) {
		tridimFluidSolver.drawContainer();
	}

}

void updateScene() {
	
	float elapsedTime = ((float)(timeGetTime() - lastTickCount)/1000.f);
	lastTickCount=timeGetTime();
	if (elapsedTime > 10.0f) // just discard the first fake calculation
		return;

	/* Process key state */
	processKeyboard();

	/* Compute frames per second */
	if (computeAndShowFPS)
		calculateFrameRate(timeGetTime());

	///////////   TEST     ////////////////
	if (simulate2D || simulate2DIn3D)
		fluidSolver.getFromUI(mouse_down, win_x, win_y, mx, my, omx, omy);

	if (simulate2D || simulate2DIn3D)
		fluidSolver.simulateStep();
	if (simulate3D)
		tridimFluidSolver.simulateStep();
    
	// Draw the next frame
    glutPostRedisplay();

}

void setupScene () {
	
	std::cout<<"Initializing scene..."<<std::endl; 

	srand(timeGetTime());

	if (simulate2D || simulate2DIn3D) {
		init2DFluidSolver();
	}

	if (simulate3D) {
		init3DFluidSolver();
	}

	/* Keyboard state reset */
	for (unsigned int j = 0; j < 256; j++)
		keys[j] = false;
	for (unsigned int j = 0; j < 248; j++)
		specialkeys[j] = false;

	// Basic GL setup..
	//glCullFace(GL_BACK);
	//glFrontFace(GL_CCW);
	//glEnable(GL_CULL_FACE);
	//glEnable(GL_DEPTH_TEST);

	if (simulate3D ) {
		//if (!loadCubeMap())
		//std::cout << "It was not possible to load cube map textures" << std::endl;
	}

	//Set up Lighting Stuff
	//glLightfv(GL_LIGHT0, GL_POSITION, lpos);
	//glEnable(GL_LIGHT0);

}

void init2DFluidSolver () {
	std::cout<<"Fluid solver init"<<std::endl;
	fluidSolver = RTP::FluidSolver(FLUID_SIDE_DIM, false, false, false);
	fluidSolver.setSimulationParameters();
}

void init3DFluidSolver () {
	std::cout<<"3d Fluid solver init"<<std::endl;
	tridimFluidSolver = RTP::FluidSolver3D(FLUID3D_EDGE_DIM, RTA::Vector3(), true, true, true);
	tridimFluidSolver.setSimulationParameters();
}

bool loadCubeMap () {

	// Set up texture maps        
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_NEAREST_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_NEAREST);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE); 
	

	// Load Cube Map images
    for (unsigned int i = 0; i < 6; i++) {        
		if (!RTR::loadPNGTextureCubeMapFace(cubemapTextures[i], cubemap[i]))
			return false;
	}
        
    glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_REFLECTION_MAP);
    glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_REFLECTION_MAP);
    glTexGeni(GL_R, GL_TEXTURE_GEN_MODE, GL_REFLECTION_MAP);

	return true;
}

void exitScene(){

    std::cout<<"Exiting scene..."<<std::endl;

    // Close window
    glutDestroyWindow(windowId);

    // Free any allocated memory

    // Exit program
    exit(0);
}

void setShaders (char* vertexShaderFn, char* fragmentShaderFn) {

	char *vs = NULL,*fs = NULL;

	v = glCreateShader(GL_VERTEX_SHADER);
	f = glCreateShader(GL_FRAGMENT_SHADER);
	
	vs = textFileRead(vertexShaderFn);
	fs = textFileRead(fragmentShaderFn);
	
	const char * ff = fs;
	const char * vv = vs;

	glShaderSource(v, 1, &vv, NULL);
	glShaderSource(f, 1, &ff, NULL);
	
	free(vs);free(fs);

	glCompileShader(v);
	glCompileShader(f);
	
	p = glCreateProgram();
	glAttachShader(p,f);
	glAttachShader(p,v);

	glLinkProgram(p);
	//glUseProgram(p);

	// shader logs
	int  vlength,    flength;
	char vlog[2048], flog[2048];
	glGetShaderInfoLog(v, 2048, &vlength, vlog);
	glGetShaderInfoLog(f, 2048, &flength, flog);
	std::cout << vlog << std::endl << std::endl << flog << std::endl << std::endl;

	switch (currentDemo) {
		case FluidDensity:

			break;
	}
}



void setViewport(int width, int height) {

    // Work out window ratio, avoid divide-by-zero when window min
    if(height==0)
		height=1;
	float ratio = float(width)/float(height);

	// Reset projection matrix
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	
	// Fill screen with viewport
	glViewport(0, 0, width, height);

	// Set a 45 degree perspective
	gluPerspective(45, ratio, .1, 1000);

}


/**********************************************************
***********************************************************
***********  INPUT CALLBACKS ****************************
***********************************************************/

/* Keyboard input callback functions */
void keypress (unsigned char key, int x, int y) {
	/* Buffer key press */
	keys[key] = true;

	if (key == 'h' || key == 'H')
		showHelp = !showHelp;

	if ((keys['c'] || keys['C'])) {
		if (simulate2D || simulate2DIn3D)
			fluidSolver.clearMatrices();
		if (simulate3D)
			tridimFluidSolver.clearMatrices();
	}

	if ((keys['g'] || keys['G'])) {
		drawFluidContainer = !drawFluidContainer;
	}

	if ((keys['d'] || keys['D'])) {
		drawFluidDensity = !drawFluidDensity;
	}

	if ((keys['v'] || keys['V'])) {
		drawFluidVelocity = !drawFluidVelocity;
	}

	
	if ((keys['m'] || keys['M'])) {
		drawMarkerParticleSet = !drawMarkerParticleSet;
	}

	if ((keys['l'] || keys['L'])) {
		lerpMarkerParticleVelocities = !lerpMarkerParticleVelocities;
		fluidSolver.switchLERPOnMarker(lerpMarkerParticleVelocities);
		tridimFluidSolver.switchLERPOnMarker(lerpMarkerParticleVelocities);
	}

	if (keys['1']) {
		simulate2DIn3D = !simulate2DIn3D;
		if (simulate2DIn3D && !fluidSolver.isReady())
			init2DFluidSolver();
	}

	if (keys['2']) {
		simulate2D = !simulate2D;
		if (simulate2D && !fluidSolver.isReady())
			init2DFluidSolver();
	}

	if (keys['3']) {
		simulate3D = !simulate3D;
		if (simulate3D && !tridimFluidSolver.isReady())
			init3DFluidSolver();
	}

}

void keyup (unsigned char key, int x, int y) {
	/* Key is not pressed any more */
	keys[key] = false;
}

void specialkeypress(int key, int xx, int yy) {
	/* Buffer special key press */
	specialkeys[key] = true;
	if (key == GLUT_KEY_F1)
		computeAndShowFPS = !computeAndShowFPS;
}

void specialkeyup (int key, int xx, int yy) {
	specialkeys[key] = false;
}

void processKeyboard () {
	
	// Test if user pressed ESCAPE (ascii 27)
	// If so, exit the program
    if(keys[27]){
		exitScene();
	}

	// 'W' key toggles wireframe mode on & off
	if(keys['w'] || keys['W']){
		wireframe=!wireframe;
		if(wireframe){
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		}
		else glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}

}

/* Mouse input callback functions */
void Mouse (int b, int s, int x, int y) {
	
	/////////    TEST     ///////////
	omx = mx = x;
	omx = my = y;

	mouse_down[b] = s == GLUT_DOWN;
	
	lastMouseX = x;
	lastMouseY = y;
	switch(b)
	{
	case GLUT_LEFT_BUTTON:
		mouseButtons[0] = ((GLUT_DOWN==s)?1:0);
		break;
	case GLUT_MIDDLE_BUTTON:
		mouseButtons[1] = ((GLUT_DOWN==s)?1:0);
		break;
	case GLUT_RIGHT_BUTTON:
		mouseButtons[2] = ((GLUT_DOWN==s)?1:0);
		break;
	default:
		break;		
	}
	glutPostRedisplay();
}

void Motion (int x, int y) {

	//////////   TEST     /////////////
	mx = x;
	my = y;
	float dx = (float)(x - lastMouseX);
	float dy = (float)(y - lastMouseY);
	lastMouseX = x;
	lastMouseY = y;

	/* Central button is currently down, or both left and right buttons are down (for those that do not have the third button */
	if (mouseButtons[1] || (mouseButtons[0] && mouseButtons[2])) {
		camera.translate(dx, dy);
		densityX = dx/(float)WINDOW_DEFAULT_WIDTH;
		densityY = ((float)WINDOW_DEFAULT_HEIGHT - dy)/(float)WINDOW_DEFAULT_HEIGHT;
	}
	/* Right button is currently down */
	else if (mouseButtons[2]) {
		camera.zoomUpdate(dx);
		vX = dx/(float)WINDOW_DEFAULT_WIDTH;
		vY = ((float)WINDOW_DEFAULT_HEIGHT - dy)/(float)WINDOW_DEFAULT_HEIGHT;
		vDelta = dy;
	}
	/* Left button is currently down */
	else if (mouseButtons[0]) {
		camera.rotate(dy, dx);
		uX = dx/(float)WINDOW_DEFAULT_WIDTH;
		uY = ((float)WINDOW_DEFAULT_HEIGHT - dy)/(float)WINDOW_DEFAULT_HEIGHT;
		uDelta = dx;
	}
	glutPostRedisplay();
}


void drawText () {
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0.0, WINDOW_DEFAULT_WIDTH, 0.0, WINDOW_DEFAULT_HEIGHT);

	
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	if (computeAndShowFPS) {
		drawFPS();
		drawStats();
	}

	/* Restore previous matrix stack */
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glColor3f(1.0f, 1.0f, 1.0f);
}

void drawFPS () {
	/* Set the color to red */
	glColor3f(1.0, 0.0, 0.0);

	int n = strlen(fPSStr);
	/* Set font and display string using glutBitmapCharacter */
	void * font = GLUT_BITMAP_9_BY_15;
	glRasterPos2i(2, WINDOW_DEFAULT_HEIGHT - 16);

	for (int i = 0; i < n; i++) {
		glutBitmapCharacter(font, fPSStr[i]);
	}


}

void drawStats () {
	/* Set the color to yellow */
	glColor3f(1.0, 1.0, 0.0);

	///* Set the window location where the text is to be displayed
	// * Raster position is set in screen coordinates, starting from the
	// * lower-left corner
	// */
	char stringBuff[150];

	bool vorticityConfinement = simulate3D? tridimFluidSolver.isVorticityConfinementOn() : fluidSolver.isVorticityConfinementOn();

	int n = sprintf_s(stringBuff, "simulate = %s, (vorticityConfinement + buoyancy) = %s, lerpMarkerParticleVelocities = %s\n", simulate3D? "3D" : "2D", vorticityConfinement? "on" : "off", lerpMarkerParticleVelocities? "on" : "off");
	/* Set font and display string using glutBitmapCharacter */
	void * font = GLUT_BITMAP_HELVETICA_12;
	glRasterPos2i(0, WINDOW_DEFAULT_HEIGHT - 44);

	for (int i = 0; i < n; i++) {
		glutBitmapCharacter(font, stringBuff[i]);
	}

	///* Set the color to green */
	//glColor3f(0.0, 1.0, 0.0);

	float diff = simulate3D ? tridimFluidSolver.getDiff() : fluidSolver.getDiff();
	float viscosity = simulate3D ? tridimFluidSolver.getViscosity() : fluidSolver.getViscosity();
	float dt = simulate3D ? tridimFluidSolver.getTimestep() : fluidSolver.getTimestep();
	float aTempr = simulate3D ? tridimFluidSolver.getAmbientTemperature() : fluidSolver.getAmbientTemperature();
	float alpha = simulate3D ? tridimFluidSolver.getBuoyancyAlphaFactor() : fluidSolver.getBuoyancyAlphaFactor();
	float beta = simulate3D ? tridimFluidSolver.getBuoyancyBetaFactor() : fluidSolver.getBuoyancyBetaFactor();
	float vorticityEpsilon = simulate3D ? tridimFluidSolver.getVorticityEps() : fluidSolver.getVorticityEps();
	n = sprintf_s(stringBuff, "diff = %4.3f, viscosity = %4.3f, dt = %3.2f\n", diff, viscosity, dt);

	/* Set font and display string using glutBitmapCharacter */
	glRasterPos2i(0, WINDOW_DEFAULT_HEIGHT - 64);

	for (int i = 0; i < n; i++) {
		glutBitmapCharacter(font, stringBuff[i]);
	}

	n = sprintf_s(stringBuff, "ambient temperature = %4.3f, buoyancy alpha = %4.3f, buoyancy beta = %4.3f, vorticity eps = %2.2f\n", aTempr, alpha, beta, vorticityEpsilon);

	/* Set font and display string using glutBitmapCharacter */
	glRasterPos2i(0, WINDOW_DEFAULT_HEIGHT - 84);

	for (int i = 0; i < n; i++) {
		glutBitmapCharacter(font, stringBuff[i]);
	}
}

/////////////////
// DEMO UPDATE //
/////////////////

void restartDemo () {




	switch (currentDemo) {
		case FluidDensity:
			//setShaders(FLUIDSVS_FILE, FLUIDSFS_FILE);
			break;
	}
}

/*
**********************************************
****************** ENTRY POINT ***************
*/

int main(int argc, char *argv[]){
        
    // Initialise OpenGL
    glutInit(&argc, argv); 

    // Set window position, size & create window
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH);
    glutInitWindowPosition(0,0);
	glutInitWindowSize(WINDOW_DEFAULT_WIDTH, WINDOW_DEFAULT_HEIGHT);
	windowId = glutCreateWindow(WINDOW_TITLE);

    // Set GLUT callback functions
    glutReshapeFunc(setViewport);
    glutDisplayFunc(renderScene);
    glutIdleFunc(updateScene);
	glutKeyboardFunc(keypress);
	glutSpecialFunc(specialkeypress);
	glutKeyboardUpFunc(keyup);
	glutSpecialUpFunc(specialkeyup);
	// Mouse handling (presses and motion WHEN buttons pressed)
	glutMouseFunc(Mouse);
	glutMotionFunc(Motion);

	/* Shader initialization */
	glewInit();
	
	if (glewIsSupported("GL_VERSION_4_2"))
		printf("Ready for OpenGL 4.2\n");
	else {
		printf("OpenGL 4.2 not supported\n");
		if (glewIsSupported("GL_VERSION_2_0"))
			printf("Ready for OpenGL 2.0\n");
		else {
			printf("OpenGL 2.0 not supported\n");
			exit(1);
		}
	}

	// Setup OpenGL state & scene resources (models, textures etc)
    setupScene();

	restartDemo();

    

    // Show window & start update loop
    glutMainLoop();    

    return 0;
    
}
