uniform samplerCube CubeMap;

varying vec3 R; // refracted vector

void main () { 
	gl_FragColor = textureCube (CubeMap , R);
}