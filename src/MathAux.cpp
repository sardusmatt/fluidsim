/*****************************************************************************************
* MathAux - implementation
* Author: Matteo Tanca [matteotanca@incomingconnection.net]
*****************************************************************************************/

#include "MathAux.h"

namespace RTA {

/********************************** QUATERNION METHOD IMPLEMENTATIONS ****************/

/* The method takes an angle (expressed in degrees) and an axis of rotation (expressed in terms
 * of its components on the three axes) and applies the equation to go from the arguments to a
 * quaternion representing the rotation
 * this.x = rX * sin (theta / 2)
 * this.y = rY * sin (theta / 2)
 * this.z = rZ * sin (theta / 2)
 * this.w = cos (theta / 2)
 */
void Quaternion::CreateFromAxisAngle (float rX, float rY, float rZ, float thetaDegrees) { 


	// From degrees to radians
	float theta = float((thetaDegrees / 180.0f) * PI);

	// sin (theta / 2)
	float sinHalfTheta = (float)sin( theta / 2.0f );

	// Calculate the x, y and z of the quaternion
	x = float(rX * sinHalfTheta);
	y = float(rY * sinHalfTheta);
	z = float(rZ * sinHalfTheta);
		
	// cos( theta / 2 )
	w = (float)cos( theta / 2.0f );

	
}

void Quaternion::addScaledVector (const Vector3& vector, float scale) { 
	Quaternion r (vector.x * scale,
				  vector.y * scale,
				  vector.z * scale,
				  0.0f);
	r *= *this;
	
	w += r.w * 0.5f;
	x += r.x * 0.5f;
	y += r.y * 0.5f;
	z += r.z * 0.5f;
}

void Quaternion::operator*= (const Quaternion& o) {
	Quaternion r = *this;
	// The correct way?
	/*w = r.w * o.w - r.x * o.x - r.y * o.y - r.z * o.z;
	x = r.w * o.x + r.x * o.w - r.y * o.z - r.z * o.y;
	y = r.w * o.y - r.x * o.z + r.y * o.w - r.z * o.x;
	z = r.w * o.z + r.x * o.y - r.y * o.x + r.z * o.w;*/
	// Millington's code
	w = r.w * o.w - r.x * o.x - r.y * o.y - r.z * o.z;
	x = r.w * o.x + r.x * o.w + r.y * o.z - r.z * o.y;
	y = r.w * o.y + r.y * o.w + r.z * o.x - r.x * o.z;
	z = r.w * o.z + r.z * o.w + r.x * o.y - r.y * o.x;
}




/********************************** BEZIER METHOD IMPLEMENTATIONS ****************/
	BezierCurve::BezierCurve () { }

	/* Since no memory allocation occurs inside the camera, default destructor
	   is fine
	 */
	//BezierCurve::~BezierCurve();
	
	/* Parameterised constructor  */
	BezierCurve::BezierCurve (Vector3 const cp[4]) {
		for (int i = 0; i < 4; i++) {
			controlPoints[i].x = cp[i].x;
			controlPoints[i].y = cp[i].y;
			controlPoints[i].z = cp[i].z;
		}
	}

	/* Updates control points */
	void BezierCurve::updateControlPoints (Vector3 const cp[4]) {
		for (int i = 0; i < 4; i++) {
			controlPoints[i].x = cp[i].x;
			controlPoints[i].y = cp[i].y;
			controlPoints[i].z = cp[i].z;
		}
	}

	/* Return the current position on the curve at the exact time t,
	* based on the Bernstein Form
	*/
	Vector3 BezierCurve::pointOnCurveAtTime (float t) {
		 //r(t) = controlPoints[0] * ( 1 - t )^3 + cP[1] * 3 * t * ( 1 - t )^2 + cP[2] * 3 * t^2 * ( 1 - t ) + cP[3] * t^3
		float tempVal1, tempVal2, tempVal3;	// to save frequently computed values and not overcomplicate the formula

		// result
		Vector3 r;

		tempVal1 = 1.0f - t;
		tempVal2 = tempVal1 * tempVal1 * tempVal1;	// (1 - t) ^ 3
		tempVal3 = t * t * t; // t ^ 3

		/* Compute the result */
		r.x = controlPoints[0].x * tempVal2 + controlPoints[1].x * 3 * t * tempVal1 * tempVal1 + controlPoints[2].x * 3 * t * t * tempVal1 + controlPoints[3].x * tempVal3;
		r.y = controlPoints[0].y * tempVal2 + controlPoints[1].y * 3 * t * tempVal1 * tempVal1 + controlPoints[2].y * 3 * t * t * tempVal1 + controlPoints[3].y * tempVal3;
		r.z = controlPoints[0].z * tempVal2 + controlPoints[1].z * 3 * t * tempVal1 * tempVal1 + controlPoints[2].z * 3 * t * t * tempVal1 + controlPoints[3].z * tempVal3;
		return r;
	}

	/* Draws the whole Bézier curve, based on the number of steps */
	void BezierCurve::draw (float lineWidth) {

		//glColor3ub(0, 0, 0);

		Vector3 p;

		// Here we tell OpenGL to render lines with a greater thickness (default is 1.0)

		glLineWidth(lineWidth); // Increase the size of a line for visibility

		glBegin(GL_LINE_STRIP);	// Start drawing lines

		/* Go through the curve starting at 0, ending at 1 + another step.
		 * Since we are using line strips, we need to go past the end point by 1 step to
		 * end at the end of the curve
		 */
		for (float t = 0; t <= (1 + (1.0f / MAX_STEPS)); t += 1.0f / MAX_STEPS) {
			// Here we pass in our 4 points that make up the curve to PointOnCurve().
			// We also pass in "t", which is the current time from 0 to 1.  If we pass
			// in 0 for t, we should get the starting point of the curve, if we pass in
			// 1 for t, we should get the end point of the curve.  So anything in between
			// 0 and 1 will be another point along the curve

			p = pointOnCurveAtTime(t);

			// Draw the current point at distance "t" of the curve.
			glVertex3f(p.x, p.y, p.z);
		}

		glEnd();

	}

	/* Draws the control points P1 red, P2 light green, P3 blue, P4 green */
	void BezierCurve::drawControlPoints (float radius) {
		
		/* Declares a quadric to draw spheres */
		GLUquadricObj *pSphere = gluNewQuadric();
		
		for (int i = 0; i < 4; i++) {
			switch (i) {
				case 0:
					// RED
					glColor3ub(255, 0, 0);
					break;
				case 1:
					// LIGHT GREEN
					glColor3ub(0, 255, 0);
					break;
				case 2:
					// BLUE
					glColor3ub(0, 0, 255);
					break;
				case 3:
					// GREEN
					glColor3ub(0, 128, 0);
					break;
			}
			DrawSphere(controlPoints[i].x, controlPoints[i].y, controlPoints[i].z, radius, pSphere);
		}
	
		/* Frees the quadric */
		gluDeleteQuadric(pSphere);
	}
	
	/* Draws a sphere with the specified parameters relying on the quadric argument */
	void BezierCurve::DrawSphere(float x, float y, float z, float radius, GLUquadric *pSphere) {
	
		/* Push a new matrix. The translations should not affect the other drawing operations */
		glPushMatrix();

		glTranslatef(x, y, z);	// Move the sphere to the desired (x, y, z)

		/* Draws the a sphere with a given radius and a width and height detail of 15 (quite round and detailed) */
		gluSphere(pSphere, radius, 15, 15);

		/* Discard the translation */
		glPopMatrix();
	}

			 			 
	
}