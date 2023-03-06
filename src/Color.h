/*****************************************************************************************
* Color - simple container struct for storing color RGB values
* Author: Matteo Tanca [matteotanca@incomingconnection.net]
*****************************************************************************************/

#ifndef COLOR_H
#define COLOR_H

namespace RTP {
	struct Color {
		float R;
		float G;
		float B;

		Color (const float red, const float green, const float blue) {
			R = red;
			G = green;
			B = blue;
		}

		Color () : R(1.0f), G(1.0f), B(1.0f) {}

		//* Used to darken (+ values) or light up the color
		void operator*= (const float factor) {
			float c = 1.0f - factor;
			R *= c;
			G *= c;
			B *= c;
		}
	};
}

#endif