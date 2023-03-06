/*****************************************************************************************
* Model - subclass of RigidBody representing OBJ model files
* Author: Matteo Tanca [matteotanca@incomingconnection.net]
*****************************************************************************************/

#ifndef MODEL_H
#define MODEL_H

#include <vector>
#include "OBJLoader.h"
#include "Utils.h"

using namespace RTA;

namespace RTR {

	class Vertex {

		public:
			Vector3 pos;	// position in local space coordinate

			// Attributes
			Vector2 uv;
			Vector3 normal;

			Vertex () {
				pos = Vector3();
				uv = Vector2();
				normal = Vector3();
			}

			Vertex (const float x, const float y, const float z) {
				pos.x = x;
				pos.y = y;
				pos.z = z;
			}

			Vertex (const Vector3& pos) {
				this->pos.x = pos.x;
				this->pos.y = pos.y;
				this->pos.z = pos.z;
			}

			const float X () const {
				return pos.x;
			}

			const float Y () const {
				return pos.y;
			}

			const float Z () const {
				return pos.z;
			}

			const float NormalX () const {
				return normal.x;
			}

			const float NormalY () const {
				return normal.y;
			}

			const float NormalZ () const {
				return normal.z;
			}

			const float U () const {
				return uv.x;
			}

			const float V () const {
				return uv.y;
			}
	};

	class Model {

		protected:

			// Vertices (pos, normals, uv map)
			std::vector<Vertex> vertices;

			float scale; // used in initialization phase to tune size

			// Flag for the correct loading of file and texture
			bool loaded;

		public:

			// Default constructor
			Model () {
				scale = 1.0f;
			}

			~Model () {}

			Model (const char* objfilename, const float scaleFactor = 1.0f) {					scale = scaleFactor;
				loaded = loadObj(objfilename);
				if (loaded) {
					scaleVertices(scale);
				}
			}

			virtual void draw (GLuint textureID) const {

				if (loaded) {

					if (textureID != NULL) {
						glEnable(GL_TEXTURE_2D);
						glBindTexture(GL_TEXTURE_2D, textureID);
					}
					
					glColor3f(1.0f, 1.0f, 1.0f);
					glBegin(GL_TRIANGLES);		

					for (unsigned int i = 0; i < vertices.size(); i++) {
						glNormal3f(vertices[i].NormalX(), vertices[i].NormalY(), vertices[i].NormalZ());
						if (textureID != NULL)
							glTexCoord2f(vertices[i].U(), vertices[i].V());
						Vector3 pos = vertices[i].pos;
						
						glVertex3f(pos.x, pos.y, pos.z);
					}
					glEnd();

					if (textureID != NULL)
						glDisable(GL_TEXTURE_2D);
					
				}
			}


			/* Scale is relative to current size...not to the original one, save for the initial invocation */
			void setScale (const float scaleFactor) {
				if (scaleFactor != 1.0f) {
					scale = scaleFactor;
					scaleVertices(scale);
				}
			}

			const float getScale () const {
				return scale;
			}


		private:

			bool loadObj (const char* objfile) {
				std::vector<Vector3> coords;
				std::vector<Vector2> uvs;
				std::vector<Vector3> normals;
				bool res = RTP::loadOBJ(objfile, coords, uvs, normals);
				if (res) {
					for (unsigned int i = 0; i < coords.size(); i++) {
						Vertex v(coords[i]);
						v.normal = normals[i];
						v.uv = uvs[i];
						vertices.push_back(v);
					}
				}
				return res;
			}

			void scaleVertices (const float scale) {
				for (unsigned int i = 0; i < vertices.size(); i++) {
					vertices[i].pos *= scale;
				}
			}



	};

}

#endif