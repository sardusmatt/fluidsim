/*****************************************************************************************
* Utils - implementation
* Author: Matteo Tanca [matteotanca@incomingconnection.net]
*****************************************************************************************/

#include "Utils.h"

namespace RTA {

	GLuint loadPNGTexture (const char* filename, GLenum minFilter,
		GLenum magFilter, GLenum wrapMode) {
            unsigned char* image;
            unsigned w, h;
            unsigned error;
            GLuint texName;
            error = LodePNG_decode32_file(&image, &w, &h, filename);
            if(error) {
				printf("Error reading in png image: %d\n", error);
				return TEXTURE_LOAD_ERROR;
            } 
			else {
				glPixelStorei(GL_UNPACK_ALIGNMENT,1);
				glGenTextures(1,&texName);
				glBindTexture(GL_TEXTURE_2D,texName);
				
				// The default GL_TEXTURE_WRAP_S and ""_WRAP_T property is
				// GL_REPEAT.
				// We need to turn this to GL_CLAMP_TO_EDGE, otherwise it
				// creates ugly seems
				// in our sky box.  GL_CLAMP_TO_EDGE does not repeat when bound
				// to an object.
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrapMode);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrapMode);
				glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER, minFilter);
				glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER, magFilter);
				gluBuild2DMipmaps(GL_TEXTURE_2D, 1, w, h, GL_RGB, GL_UNSIGNED_BYTE, image);
	
				glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,w,								h,0, GL_RGBA,GL_UNSIGNED_BYTE,image);
            }
            free(image);
            return texName;
	}

}

namespace RTR {

	bool loadPNGTextureCubeMapFace (const char* filename, GLenum faceID) {
		unsigned char* image;
		unsigned w, h;
		unsigned error;

		error = LodePNG_decode32_file(&image, &w, &h, filename);
		if(error) {
			printf("Error reading in png cubemap: %d\n", error);
			return false;
		} 
		else {
			glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_GENERATE_MIPMAP, GL_TRUE);
			glTexImage2D(faceID, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, image);
			free(image);
		}

		return true;
	}
}