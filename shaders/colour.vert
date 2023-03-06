// simple colour shader
// The gl_Color is the attribute value
// The gl_FrontColor and gl_BackColor are then determined in the vertex shader
// and used to compute the varying gl_Color for the fragment shader by interpolation

void main() {	
	gl_FrontColor = gl_Color;
	gl_Position = ftransform();
}
