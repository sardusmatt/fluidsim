/*****************************************************************************************
* Utils.h - Provides utility methods and defines used throughout the code
* Author: Matteo Tanca [matteotanca@incomingconnection.net]
*****************************************************************************************/

#ifndef UTILS_H
#define UTILS_H

#include "lodepng.h"
#include "glew.h"
#include "glut.h"	// GLUT

 
using namespace std;

namespace RTA {

	#define TEXTURE_LOAD_ERROR 0

	 /*
	  * Loads a png file into an opengl texture object, using cstdio , libpng,
	  * and opengl.
	  * filename : the png file to be loaded
	  * GLuint : return an opengl texture id.  Will be 0 if there is a major
	  * error
	  */
	GLuint loadPNGTexture (const char* filename, GLenum minFilter = GL_NEAREST,
		GLenum magFilter = GL_NEAREST, GLenum wrapMode = GL_CLAMP_TO_EDGE);

	 
 }

namespace RTR {
	/* Loads png files into a cube map
	  * filename : the png file to be loaded
	  * faceID : di of the face of the cube
	  */
	bool loadPNGTextureCubeMapFace (const char* filename, GLenum faceID);
}
#endif //UTILS