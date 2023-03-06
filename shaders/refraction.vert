
uniform float eta; // eta ratio

varying vec3 R; // refract vector

void main () {
	// Create incident and normal vectors
	vec4 V = gl_ModelViewMatrix * gl_Vertex;
	vec4 E = gl_ProjectionMatrixInverse	* vec4 (0 ,0 , -1 ,0);
	vec3 I = normalize (V.xyz - E.xyz);
	vec3 N = normalize ( gl_Normal );
	R = refract (I, N, eta);
	//R.y = -R.y;
	//R.x = -R.x;
	gl_Position = ftransform ();
}
