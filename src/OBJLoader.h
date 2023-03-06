/*****************************************************************************************
* OBJLoader - Loader for OBJ format, supports textures and normals
* Author: Matteo Tanca [matteotanca@incomingconnection.net]
*****************************************************************************************/

#ifndef OBJLOADER_H
#define OBJLOADER_H

#include <string>
#include <cstring>
#include <vector>
#include <iostream>		// I/O
#include "MathAux.h"


using namespace RTA;

namespace RTP {

	/* True if loading fine */
	bool loadOBJ (const char * path, std::vector<Vector3>& vertices, 
		std::vector<Vector2>& uVCoords, std::vector<Vector3>& normals);


}

#endif