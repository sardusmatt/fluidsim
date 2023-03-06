/*****************************************************************************************
* Camera - a simple camera class enabling Maya-like pan, zoom, and strafe
* Author: Matteo Tanca [matteotanca@incomingconnection.net]
*****************************************************************************************/

#ifndef CAMERA_H
#define CAMERA_H


namespace RTA {

	/* Basic constraints to limit camera movements w.r.t. plane position */ 
	struct CameraConstraints {
		float minZoom;
		float maxZoom;
		float minXRot;
		float maxXRot;
		float minYRot;
		float maxYRot;
		float minXTran;
		float maxXTran;
		float minYTran;
		float maxYTran;
	};

	class Camera {
		private:
			// Camera zoom
			float zoom;
			// Rotations around the y and x axis
			float xRot;
			float yRot;
			// Translations along the x and y axis
			float xTran;
			float yTran;
			// Camera constraints
			CameraConstraints constr;

		public:

			/* Step factors */
			static const float ZOOMSTEP;
			static const float ROTATIONSTEP;
			static const float TRANSLATIONSTEP;

			/* default constructor creates a camera with all parameters set to 0 */
			 Camera ();

			 /* Since no memory allocation occurs inside the camera, default destructor
			    is fine
			  */
			 //~Camera();

			 /* parameterised constructor  */
			 Camera (const float cameraZoom, const float xRotation, const float yRotation, const float xTranslation, const float yTranslation);
			 
			 /* Subtracts the zoomDelta from the current value of the zoom instance var */
			 void zoomUpdate (const float zoomDelta);

			 /* Rotates the camera view based on the argument rotations for x and y axis */
			 void rotate (const float xRotationDelta, const float yRotationDelta);

			 /* Translates the camera view based on the arguments */
			 void translate (const float xTranslationDelta, const float yTranslationDelta);
			 
			 /* Sets basic camera constraints */
			 void setCameraConstraints (CameraConstraints const& cameraConstraints);

			 /* Inline methods */
			 float getZoom () const {return zoom;}
			 float getXRotation () const {return xRot;}
			 float getYRotation () const {return yRot;}
			 float getXTranslation () const {return xTran;}
			 float getYTranslation () const {return yTran;}
			 


	};

}


#endif //CAMERA