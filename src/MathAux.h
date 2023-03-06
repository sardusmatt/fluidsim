/*****************************************************************************************
* MathAux.h - Defines a few support math classes
	- 3D Vector
	- Beziér Curve
	- Quaternion
* Author: Matteo Tanca [matteotanca@incomingconnection.net]
*****************************************************************************************/

#ifndef MATH_AUX_H
#define MATH_AUX_H

#include <math.h>
#include <stdlib.h>
#include "glut.h"	// GLUT

namespace RTA {

	/* PI */
	#define  PI 3.14159265

	#define MIN(a,b)     (((a) < (b)) ? (a) : (b))
	#define MAX(a,b)     (((a) < (b)) ? (b) : (a))

	/* Linear interpolation macro */
	#define LERP(a,b,c)     (((b) - (a)) * (c) + (a))

	const float radiusEpsilon = 1e-4f;   // NOTE: used to avoid numerical inaccuracies when computing Welzl's smallest bounding sphere

	const float EPSILON = 1e-4f;

	/* Random float generator. Pre-requisite: call srand once before use */
	static float RandomFloat (const float a, const float b) {
		float random = ((float) rand()) / (float) RAND_MAX;
		float diff = b - a;
		float r = random * diff;
		return a + r;
	}

	static float Clamp (const float f, const float min, const float max) {
		if (f < min) return min;
		if (f > max) return max;
		return f;
	}

	/* Random int generator. Pre-requisite: call srand once before use */
	static float RandomInt (const int a, const int b) {
		return (rand() % (b - a)) + a;
	}

	///// Inverts a 4x4 float matrix /////////////
	static bool invert4x4 (const float m[16], float invOut[16]) {
		float inv[16], det;
		int i;
		inv[0] = m[5]  * m[10] * m[15] - 
			m[5]  * m[11] * m[14] - 
			m[9]  * m[6]  * m[15] + 
			m[9]  * m[7]  * m[14] +
			m[13] * m[6]  * m[11] - 
			m[13] * m[7]  * m[10];

		inv[4] = -m[4]  * m[10] * m[15] + 
			m[4]  * m[11] * m[14] + 
			m[8]  * m[6]  * m[15] - 
			m[8]  * m[7]  * m[14] - 
			m[12] * m[6]  * m[11] + 
			m[12] * m[7]  * m[10];

		inv[8] = m[4]  * m[9] * m[15] - 
			m[4]  * m[11] * m[13] - 
			m[8]  * m[5] * m[15] + 
			m[8]  * m[7] * m[13] + 
			m[12] * m[5] * m[11] - 
			m[12] * m[7] * m[9];

		inv[12] = -m[4]  * m[9] * m[14] + 
			m[4]  * m[10] * m[13] +
			m[8]  * m[5] * m[14] - 
			m[8]  * m[6] * m[13] - 
			m[12] * m[5] * m[10] + 
			m[12] * m[6] * m[9];

		inv[1] = -m[1]  * m[10] * m[15] + 
			m[1]  * m[11] * m[14] + 
			m[9]  * m[2] * m[15] - 
			m[9]  * m[3] * m[14] - 
			m[13] * m[2] * m[11] + 
			m[13] * m[3] * m[10];

		inv[5] = m[0]  * m[10] * m[15] - 
			m[0]  * m[11] * m[14] - 
			m[8]  * m[2] * m[15] + 
			m[8]  * m[3] * m[14] + 
			m[12] * m[2] * m[11] - 
			m[12] * m[3] * m[10];

		inv[9] = -m[0]  * m[9] * m[15] + 
			m[0]  * m[11] * m[13] + 
			m[8]  * m[1] * m[15] - 
			m[8]  * m[3] * m[13] - 
			m[12] * m[1] * m[11] + 
			m[12] * m[3] * m[9];

		inv[13] = m[0]  * m[9] * m[14] - 
			m[0]  * m[10] * m[13] - 
			m[8]  * m[1] * m[14] + 
			m[8]  * m[2] * m[13] + 
			m[12] * m[1] * m[10] - 
			m[12] * m[2] * m[9];

		inv[2] = m[1]  * m[6] * m[15] - 
			m[1]  * m[7] * m[14] - 
			m[5]  * m[2] * m[15] + 
			m[5]  * m[3] * m[14] + 
			m[13] * m[2] * m[7] - 
			m[13] * m[3] * m[6];

		inv[6] = -m[0]  * m[6] * m[15] + 
			m[0]  * m[7] * m[14] + 
			m[4]  * m[2] * m[15] - 
			m[4]  * m[3] * m[14] - 
			m[12] * m[2] * m[7] + 
			m[12] * m[3] * m[6];

		inv[10] = m[0]  * m[5] * m[15] - 
			m[0]  * m[7] * m[13] - 
			m[4]  * m[1] * m[15] + 
			m[4]  * m[3] * m[13] + 
			m[12] * m[1] * m[7] - 
			m[12] * m[3] * m[5];

		inv[14] = -m[0]  * m[5] * m[14] + 
			m[0]  * m[6] * m[13] + 
			m[4]  * m[1] * m[14] - 
			m[4]  * m[2] * m[13] - 
			m[12] * m[1] * m[6] + 
			m[12] * m[2] * m[5];

		inv[3] = -m[1] * m[6] * m[11] + 
			m[1] * m[7] * m[10] + 
			m[5] * m[2] * m[11] - 
			m[5] * m[3] * m[10] - 
			m[9] * m[2] * m[7] + 
			m[9] * m[3] * m[6];

		inv[7] = m[0] * m[6] * m[11] - 
			m[0] * m[7] * m[10] - 
			m[4] * m[2] * m[11] + 
			m[4] * m[3] * m[10] + 
			m[8] * m[2] * m[7] - 
			m[8] * m[3] * m[6];

		inv[11] = -m[0] * m[5] * m[11] + 
			m[0] * m[7] * m[9] + 
			m[4] * m[1] * m[11] - 
			m[4] * m[3] * m[9] - 
			m[8] * m[1] * m[7] + 
			m[8] * m[3] * m[5];

		inv[15] = m[0] * m[5] * m[10] - 
			m[0] * m[6] * m[9] - 
			m[4] * m[1] * m[10] + 
			m[4] * m[2] * m[9] + 
			m[8] * m[1] * m[6] - 
			m[8] * m[2] * m[5];

		det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

		if (det == 0)
			return false;

		det = 1.0 / det;

		for (i = 0; i < 16; i++)
			invOut[i] = inv[i] * det;

		return true;
	}

	/* Basic 3DVector representation */ 
	class Vector3 {
		public:
			float x;
			float y;
			float z;

			Vector3 () : x(0.0f), y(0.0f), z(0.0f) {}

			Vector3 (const float cpx, const float cpy, const float cpz) : x(cpx), y(cpy), z(cpz) {}

			void normalize () {

				float norm = magnitude();

				/* Sanity check to avoid dividing by zero */
				if (norm > 0) {
					x /= norm;
					y /= norm;
					z /= norm;
				}
			}

			void invert () {
				x = -x;
				y = -y;
				z = -z;
			}

			void clear () {
				x = 0.0f;
				y = 0.0f;
				z = 0.0f;
			}

			// Takes into account a small tolerance threshold
			bool isNotZero () const {
				return (abs(x) > EPSILON ||
						abs(y) > EPSILON ||
						abs(z) > EPSILON    );
			}

			void truncate (const float maxValue) {
				if (abs(x) > maxValue)
					x = (x > 0.0f) ? maxValue : - maxValue;
				if (abs(y) > maxValue)
					y = (y > 0.0f) ? maxValue : - maxValue;
				if (abs(z) > maxValue)
					z = (z > 0.0f) ? maxValue : - maxValue;
			}

			float magnitude () const {
				return sqrtf(x*x + y*y + z*z);
			}

			/* Can be used for comparisons between magnitudes. Relative order is preserved
			 * but square root extra computational cost is avoided
			 */
			float squaredMagnitude () const {
				return x*x + y*y + z*z;
			}

			Vector3 normalized () const {
				Vector3 norma = Vector3(x, y, z);
				norma.normalize();
				return norma;
			}

			/***************************************************
			 ************* OPERATOR OVERLOADS ******************
			 **************************************************/
			
			// +
			Vector3 operator+ (const Vector3& vector) const {
				return Vector3(vector.x + x, vector.y + y, vector.z + z);
			}

			void operator+= (const Vector3& vector) {
				x += vector.x;
				y += vector.y;
				z += vector.z;
			}

			bool operator== (const Vector3& vector) const {
				return (x == vector.x
					&&	y == vector.y
					&&	z == vector.z);
			}

			void operator= (const Vector3& vector) {
				x = vector.x;
				y = vector.y;
				z = vector.z;
			}

			// -
			Vector3 operator- (const Vector3& vector) const {
				return Vector3(x - vector.x, y - vector.y, z - vector.z);
			}

			void operator-= (const Vector3& vector) {
				x -= vector.x;
				y -= vector.y;
				z -= vector.z;
			}

			// *
			Vector3 operator* (const float scalar) const {
				return Vector3(x * scalar, y * scalar, z * scalar);
			}

			/* Dot */
			float operator* (const Vector3& vector) const {
				return (x * vector.x + y * vector.y + z * vector.z);
			}

			void operator*= (const float scalar) {
				x *= scalar;
				y *= scalar;
				z *= scalar;
			}

			// /
			Vector3 operator/ (const float scalar) const {
				return Vector3(x/scalar, y/scalar, z/scalar);
			}

			// / =
			void operator/= (const float scalar) {
				x /= scalar;
				y /= scalar;
				z /= scalar;
			}

			/* Vector (X) product */
			Vector3 vectorProduct(const Vector3 &vector) const {
				return Vector3(y*vector.z-z*vector.y, 
					           z*vector.x-x*vector.z, 
						       x*vector.y-y*vector.x);
			}
	};

	/* Returns a random vector with x, y, z in the range of the
	 * corresponding components of a and b
	 */
	static Vector3 RandomVector (const Vector3& a, const Vector3& b) {
		Vector3 result;
		result.x = RandomFloat(a.x, b.x);
		result.y = RandomFloat(a.y, b.y);
		result.z = RandomFloat(a.z, b.z);
		return result;
	}


	/* Based on Ericson's pseudocode. Computes closest points c1, c2, of segments
	 * s1(s) = p1 + s*(q1 - p1) and s2(t) = p2 + t*(q2 - p2), returning s and t.
	 * Return value is the squared distance between s1(s) and s2(t)
	 */
	static float ClosestPtSegmentSegment (const Vector3& p1, const Vector3& q1, const Vector3& p2, const Vector3& q2, float &s, float &t, Vector3& c1, Vector3& c2) {

		// First calculate the direction vectors for s1, s2
		Vector3 d1 = q1 - p1;
		Vector3 d2 = q2 - p2;

		Vector3 r = p1 - p2;

		// Squared magnitude of s1, s2
		float a = d1.squaredMagnitude();
		float e = d2.squaredMagnitude();

		float f = d2 * r;

		// Check possible segment to point degeneration
		if (a <= EPSILON && e <= EPSILON) {
			// Both segemnts are points
			s = t = 0.0f;
			c1 = p1;
			c2 = p2;

			return (c1 - c2).squaredMagnitude();
		}

		if (a <= EPSILON) {
			// s1 is a point segment
			s = 0.0f;
			t = f / e;	// s = 0 -> t = (b * s + f) / e
			t = Clamp(t, 0.0f, 1.0f);
		}
		else {
			float c = d1 * r;
			if (e <= EPSILON) {
				// s2 degenerates into point
				t = 0.0f;
				s = Clamp(-c / a, 0.0f, 1.0f);
			}
			else {
				float b = d1 * d2;
				float den = a * e - b * b;


				// Segments not parallel. Computes closest point on l1 to l2 and clamps to
				// s1. Otherwise pick an arbitrary s (0)
				if (den != 0.0f) {
					s = Clamp((b * f - c * e) / den, 0.0f, 1.0f);
				}
				else
					s = 0.0f;

				// Computes the closest point on l2 to s1(s)
				t = (b * s + f) / e;

				// Clamps t and recomputes s if need be
				if (t < 0.0f) {
					t = 0.0f;
					s = Clamp(-c / a, 0.0f, 1.0f);
				}
				else if ( t> 1.0f) {
					t = 1.0f;
					s = Clamp((b - c) / a, 0.0f, 1.0f);
				}
			}
		}

		c1 = p1 + d1 * s;
		c2 = p2 + d2 * t;
		return (c1 - c2).squaredMagnitude();
	}


/*********************** VECTOR2 (UV coordinates) ********************/
	
	class Vector2 {
		public:
			float x;
			float y;

			Vector2 () : x(0.0f), y(0.0f) {}

			Vector2 (const float cpx, const float cpy) : x(cpx), y(cpy){}

	};

/*********************** QUATERION************************************/

	class Quaternion {
		
		public:
			float x;	// first complex component
			float y;	// second complex component
			float z;	// third complex component

			float w;	// real comp of the quaternion

			Quaternion () : x(0.0f), y(0.0f), z(0.0f), w(1.0f) {}

			Quaternion (float qx, float qy, float qz, float qw) : x(qx), y(qy), z(qz), w(qw) {}


			/* converts an axis angle rotation to a quaternion */
			void CreateFromAxisAngle(float rX, float rY, float rZ, float thetaDegrees);

			void rotateByVector(const Vector3& vector) {
				Quaternion q(vector.x, vector.y, vector.z, 0.f);
				(*this) *= q;
			}

			/* adds a scaled vector to the quaternion. Used to update current
			 * orientation (the quaternion) by an angular velocity (the vector) scaled
			 * by an amount of time
			 */
			void addScaledVector (const Vector3& vector, float scale);

			/* Combines this and o by multiplying them and storing the result in this
			 * As with transformation matrices, the result is equivalent to apply first
			 * the second factor (that is, the rotation o) and then the first one (this)
			 */
			void operator*= (const Quaternion& o);

			/* Normalises the quaternion to unit length, making it a valid
			 * orientation quaternion
			 */
			void normalise () {
				float d = w*w+x*x+y*y+z*z;

				// Check for zero length quaternion, and use the no-rotation
				// quaternion in that case.
				if (d == 0) { 
					w = 1; 
					return;
				}

				d = (1.0f)/sqrtf(d);
				w *= d;
				x *= d;
				y *= d;
				z *= d;
			}

	};


/*********************** MATRIX3 *************************************/

	class Matrix3 {
		public:
			float elements[9];	// stored row by row

			// Matrix filled with zeroes
			Matrix3 () {
				for (unsigned int i = 0; i < 9; i++)
					elements[i] = 0.0f;
			}

			Matrix3 (const float* elByRow) {
				for (unsigned int i = 0; i < 9; i++)
					elements[i] = elByRow[i];
			}

			/* this (rotation matrix) is set to correspond to the argument */
			void setOrientation (const Quaternion& q) {
				elements[0] = 1.0f - 2.0f * ( q.y * q.y + q.z * q.z ); 
				elements[1] = 2.0f * (q.x * q.y + q.z * q.w);
				elements[2] = 2.0f * (q.x * q.z - q.y * q.w);

				// Row 2
				elements[3] = 2.0f * (q.x * q.y - q.z * q.w);  
				elements[4] = 1.0f - 2.0f * (q.x * q.x + q.z * q.z); 
				elements[5] = 2.0f * (q.y * q.z + q.x * q.w);  

				// Row 3
				elements[6] = 2.0f * (q.x * q.z + q.y * q.w);
				elements[7] = 2.0f * (q.y * q.z - q.x * q.w);
				elements[8] = 1.0f - 2.0f * (q.x * q.x + q.y * q.y);  
				
			};

			Vector3 operator* (const Vector3& vector) const {
				return Vector3(vector.x * elements[0] +
							   vector.y * elements[1] +
							   vector.z * elements[2],

							   vector.x * elements[3] +
							   vector.y * elements[4] +
							   vector.z * elements[5],

							   vector.x * elements[6] +
							   vector.y * elements[7] +
							   vector.z * elements[8]);
			}

			void setIdentity () {
				elements[1] = elements[2] = elements[3] = elements[5] = elements[7] = elements[6] = 0.0f;
            elements[0] = elements[4] = elements[8] = 1.0f;
			}

			Matrix3 transpose () const {
				float byCol[9] = {elements[0], elements[3], elements[6], elements[1], elements[4], elements[7], elements[2], elements[5], elements[8]};
				return Matrix3(byCol);
			}

			 Matrix3 operator*(const Matrix3 &o) const {
				 float elems[9] = {elements[0]*o.elements[0] + elements[1]*o.elements[3] + elements[2]*o.elements[6],
					elements[0]*o.elements[1] + elements[1]*o.elements[4] + elements[2]*o.elements[7],
					elements[0]*o.elements[2] + elements[1]*o.elements[5] + elements[2]*o.elements[8],

					elements[3]*o.elements[0] + elements[4]*o.elements[3] + elements[5]*o.elements[6],
					elements[3]*o.elements[1] + elements[4]*o.elements[4] + elements[5]*o.elements[7],
					elements[3]*o.elements[2] + elements[4]*o.elements[5] + elements[5]*o.elements[8],

					elements[6]*o.elements[0] + elements[7]*o.elements[3] + elements[8]*o.elements[6],
					elements[6]*o.elements[1] + elements[7]*o.elements[4] + elements[8]*o.elements[7],
					elements[6]*o.elements[2] + elements[7]*o.elements[5] + elements[8]*o.elements[8]};
				 return Matrix3(elems);
			}

			/* Computation of a 3x3 determinant */
			static float Matrix3::determinant(float m11, float m12, float m13, 
											  float m21, float m22, float m23, 
											  float m31, float m32, float m33)	{
				return m11 * (m22 * m33 - m32 * m23) -
					   m21 * (m12 * m33 - m32 * m13) +
					   m31 * (m12 * m23 - m22 * m13);
			}

			Matrix3 getInverse () const {
				Matrix3 r;

				float t4 = elements[0] * elements[4];
				float t6 = elements[0] * elements[5];
				float t8 = elements[1] * elements[3];
				float t10 = elements[2] * elements[3];
				float t12 = elements[1] * elements[6];
				float t14 = elements[2] * elements[6];

				// Calculate the determinant
				float t16 = (t4 * elements[8] - t6 * elements[7] - t8 * elements[8] +
					t10 * elements[7] + t12 * elements[5] - t14 * elements[4]);

				// Make sure the determinant is non-zero.
				if (t16 == 0.0f)
					return r;
				float t17 = 1/t16;

				r.elements[0] = (elements[4]*elements[8]-elements[5]*elements[7])*t17;
				r.elements[1] = -(elements[1]*elements[8]-elements[2]*elements[7])*t17;
				r.elements[2] = (elements[1]*elements[5]-elements[2]*elements[4])*t17;
				r.elements[3] = -(elements[3]*elements[8]-elements[5]*elements[6])*t17;
				r.elements[4] = (elements[0]*elements[8]-t14)*t17;
				r.elements[5] = -(t6-t10)*t17;
				r.elements[6] = (elements[3]*elements[7]-elements[4]*elements[6])*t17;
				r.elements[7] = -(elements[0]*elements[7]-t12)*t17;
				r.elements[8] = (t4-t8)*t17;

				return r;
			}

	};


/*********************** MATRIX4 ************************************/

	class Matrix4 {
		public:
			float elements[12];	// last row/four elements are assumed (0,0,0,1) to produce a homogeneous 4x4 matrix

			/* Identity matrix */
			Matrix4 () {
				elements[1] = elements[2] = elements[3] = elements[4] = elements[6] = elements[7] = elements[8] = elements[9] = elements[11] = 0.0f;
            elements[0] = elements[5] = elements[10] = 1.0f;
			}

			float& operator[] (unsigned int index) {
				return elements[index];
			}

			Vector3 operator* (const Vector3& vector) const {
				return Vector3(vector.x * elements[0] +
							   vector.y * elements[1] +
							   vector.z * elements[2] +
							   elements[3],	// implicit

							   vector.x * elements[4] +
							   vector.y * elements[5] +
							   vector.z * elements[6] +
							   elements[7],	// implicit

							   vector.x * elements[8] +
							   vector.y * elements[9] +
							   vector.z * elements[10] +
							   elements[11]);	// implicit
			}

			/* this (rotation matrix) is set to correspond to the arguments (the quaternion gives the orientation, whilst the vector the position  */
			void setOrientationAndPosition (const Quaternion& q, const Vector3& p) {
				elements[0] = 1.0f - 2.0f * (q.y * q.y + q.z * q.z); 
				elements[1] = 2.0f * (q.x * q.y + q.z * q.w);
				elements[2] = 2.0f * (q.x * q.z - q.y * q.w);
				elements[3] = p.x;

				// Row 2
				elements[4] = 2.0f * (q.x * q.y - q.z * q.w);  
				elements[5] = 1.0f - 2.0f * (q.x * q.x + q.z * q.z); 
				elements[6] = 2.0f * (q.y * q.z + q.x * q.w);  
				elements[7] = p.y;

				// Row 3
				elements[8] = 2.0f * (q.x * q.z + q.y * q.w);
				elements[9] = 2.0f * (q.y * q.z - q.x * q.w);
				elements[10] = 1.0f - 2.0f * (q.x * q.x + q.y * q.y);  
				elements[11] = p.z;
				
			};

	};

/**************************** BEZIER**********************************/

	#define MAX_STEPS 25.0f	// This is the amount of steps we want to draw the curve in.

	class BezierCurve {
		private:
			// Control points. 4 of them, so cubic curve
			Vector3 controlPoints[4];

			/* Draws a sphere with the specified parameters relying on the quadric argument */
			void DrawSphere (float x, float y, float z, float radius, GLUquadric *pSphere);

		public:

			/* default constructor creates a BezierCurve with associated control points
			 * all set to the origin of the axes */
			 BezierCurve ();

			 /* Since no memory allocation occurs inside the class, default destructor
			    is fine
			  */
			 //~BezierCurve();

			 /* parameterised constructor  */
			 BezierCurve (Vector3 const cp[4]);

			 /* Updates control points */
			 void updateControlPoints (Vector3 const cp[4]);

			 /* Return the current position on the curve at the exact time t,
			  * based on the Bernstein Form
			  */
			 Vector3 pointOnCurveAtTime (float t);

			 /* Draws the whole Bézier curve, based on the number of steps */
			 void draw (float lineWidth = 1.5f);

			 /* Draws the control points P1 red, P2 light green, P3 blue, P4 green */
			 void drawControlPoints (float radius = 1.5f);
			 			 
	};

}


#endif //MATH_AUX