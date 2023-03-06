#include "OBJLoader.h"

namespace RTP {

	/* True if loading fine */
	bool loadOBJ (const char * path, std::vector<Vector3>& vertices,
		std::vector<Vector2>& uVCoords, std::vector<Vector3>& normals) {

		printf("Loading OBJ file %s...\n", path);
		std::vector<unsigned int> vertexIndices, uvIndices, normalIndices;
        std::vector<Vector3> temp_vertices; 
        std::vector<Vector2> temp_uvs;
        std::vector<Vector3> temp_normals;

		/* Open the file in read only mode */
		FILE * file = fopen(path, "r");
        if (file == NULL) {
			std::cout << "Cannot open file in OBJLoader" << std::endl;
            return false;
        }

		bool eOf = false;

		/* Read the OBJ file content */
		while (!eOf) {
			char lineHeader[128];
            // read the first word of the line
            int res = fscanf(file, "%s", lineHeader);
            if (res == EOF)
                eOf = true;
			else {
				// Vertices
				if ( strcmp( lineHeader, "v" ) == 0 ){
					Vector3 vertex;
					fscanf(file, "%f %f %f\n", &vertex.x, &vertex.y, &vertex.z );
					temp_vertices.push_back(vertex);
				// Texture coord
				}else if ( strcmp( lineHeader, "vt" ) == 0 ){
					Vector2 uv;
					fscanf(file, "%f %f\n", &uv.x, &uv.y );
					// Would make sense if we were using dds
					//uv.y = -uv.y;
					temp_uvs.push_back(uv);
				// Normals
				}else if ( strcmp( lineHeader, "vn" ) == 0 ){
					Vector3 normal;
					fscanf(file, "%f %f %f\n", &normal.x, &normal.y, &normal.z );
					temp_normals.push_back(normal);
				}
				// Face
				else if ( strcmp( lineHeader, "f" ) == 0 ){
					std::string vertex1, vertex2, vertex3;
					unsigned int vertexIndex[3], uvIndex[3], normalIndex[3];
					int matches = fscanf(file, "%d/%d/%d %d/%d/%d %d/%d/%d\n", &vertexIndex[0], &uvIndex[0], &normalIndex[0], &vertexIndex[1], &uvIndex[1], &normalIndex[1], &vertexIndex[2], &uvIndex[2], &normalIndex[2] );
					if (matches != 9){
						std::cout << "File can't be read by our simple parser. -_*" << std::endl;
						return false;
					}
					vertexIndices.push_back(vertexIndex[0]);
					vertexIndices.push_back(vertexIndex[1]);
					vertexIndices.push_back(vertexIndex[2]);
					uvIndices.push_back(uvIndex[0]);
					uvIndices.push_back(uvIndex[1]);
					uvIndices.push_back(uvIndex[2]);
					normalIndices.push_back(normalIndex[0]);
					normalIndices.push_back(normalIndex[1]);
					normalIndices.push_back(normalIndex[2]);
				}
				else {
					// Get rid of comments in the OBJ file (assuming it is a comment) 
					char commentBuffer[1000];
					fgets(commentBuffer, 1000, file);
				}	
			}
		}

		/* And now use the info to fill the buffers that were passed as arguments */
		// For each vertex of each triangle (note that OBJ indices start from 1, so
		// the minus one
        for (unsigned int i = 0; i < vertexIndices.size(); i++) {

			// Get the indices of its attributes
			unsigned int vertexIndex = vertexIndices[i];
			unsigned int uvIndex = uvIndices[i];
			unsigned int normalIndex = normalIndices[i];

			// Get the attributes through the index
			Vector3 vertex = temp_vertices[vertexIndex-1];
			Vector2 uv = temp_uvs[uvIndex-1];
			Vector3 normal = temp_normals[normalIndex-1];

			// Put the attributes in buffers
			vertices.push_back(vertex);
			uVCoords.push_back(uv);
			normals .push_back(normal);
        
        }

		fclose(file);
        return true;	// lkoaded


	}


}


