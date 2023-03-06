/*****************************************************************************************
* Camera - implementation
* Author: Matteo Tanca [matteotanca@incomingconnection.net]
*****************************************************************************************/

#include "Camera.h"

namespace RTA {

	/* Float const initialisation */
	const float Camera::ZOOMSTEP = 0.05f;
	const float Camera::ROTATIONSTEP = 0.5f;
	const float Camera::TRANSLATIONSTEP = 0.05f;

	/* Default constructor creates a camera with all parameters set to 0 */
	Camera::Camera () : zoom(0.0f), xRot(0.0f), xTran(0.0f), yRot(0.0f), yTran(0.0f) { }

	/* Since no memory allocation occurs inside the camera, default destructor
	   is fine
	 */
	//Camera::~Camera();
	
	/* Parameterised constructor  */
	Camera::Camera (const float cameraZoom, const float xRotation, const float yRotation, const float xTranslation, const float yTranslation) :
			zoom(cameraZoom), xRot(xRotation), xTran(xTranslation), yRot(yRotation), yTran(yTranslation) { }
		
	/* Subtracts the zoomDelta from the current value of the zoom instance var */
	void Camera::zoomUpdate (const float zoomDelta) {
		float v = zoom - ZOOMSTEP * zoomDelta;
		/*if (v < constr.minZoom)
			zoom = constr.minZoom;
		else if (v > constr.maxZoom)
			zoom = constr.maxZoom;
		else*/
			zoom = v;
	}
	
	/* Rotates the camera view based on the argument rotations for x and y axis */
	void Camera::rotate (const float xRotationDelta, const float yRotationDelta) {
		float v = xRot + ROTATIONSTEP * xRotationDelta;
		/* normalisation for angles >|< 360°/-360° */
		if (v > 360.0f)
			v -= 360.0f;
		if (v < -360.0f)
			v += 360.0f;

		/*if (v < constr.minXRot)
			xRot = constr.minXRot;
		else if (v > constr.maxXRot)
			xRot = constr.maxXRot;
		else*/
			xRot = v;

		v = yRot + ROTATIONSTEP * yRotationDelta;
		/* normalisation for angles >|< 360°/-360° */
		if (v > 360.0f)
			v -= 360.0f;
		if (v < -360.0f)
			v += 360.0f;

		/*if (v < constr.minYRot)
			yRot = constr.minYRot;
		else if (v > constr.maxYRot)
			yRot = constr.maxYRot;
		else*/
			yRot = v;
	}

	/* Translates the camera view based on the arguments */
	void Camera::translate (const float xTranslationDelta, const float yTranslationDelta) {
		float v = xTran + TRANSLATIONSTEP * xTranslationDelta;
		/*if (v < constr.minXTran)
			xTran = constr.minXTran;
		else if (v > constr.maxXTran)
			xTran = constr.maxXTran;
		else*/
			xTran = v;

		v = yTran + TRANSLATIONSTEP * yTranslationDelta;
		/*if (v < constr.minYTran)
			yTran = constr.minYTran;
		else if (v > constr.maxYTran)
			yTran = constr.maxYTran;
		else*/
			yTran = v;
	}

	/* Sets basic camera constraints */
	void Camera::setCameraConstraints (const CameraConstraints& cameraConstraints) {
		constr.minZoom = cameraConstraints.minZoom;
		constr.maxZoom = cameraConstraints.maxZoom;
		constr.minXRot = cameraConstraints.minXRot;
		constr.maxXRot = cameraConstraints.maxXRot;
		constr.minYRot = cameraConstraints.minYRot;
		constr.maxYRot = cameraConstraints.maxYRot;
		constr.minXTran = cameraConstraints.minXTran;
		constr.maxXTran = cameraConstraints.maxXTran;
		constr.minYTran = cameraConstraints.minYTran;
		constr.maxYTran = cameraConstraints.maxYTran;
	}
}