// Toon shader second ver (the vertex form just sends the normal to the fragment



varying vec3 normal;

void main() {
	normal = gl_Normal;
	gl_Position = ftransform();
}




