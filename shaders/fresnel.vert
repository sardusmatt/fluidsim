uniform float eta; // Ratio of indices of refraction
uniform float fresnelPower;
varying vec3 Reflect;
varying vec3 Refract;
varying float Ratio;
void main() {
	float F = ((1.0-eta) * (1.0-eta)) / ((1.0+eta) * (1.0+eta));
	// Create incident and normal vectors
	vec4 V = gl_ModelViewMatrix * gl_Vertex;
	vec4 E = gl_ProjectionMatrixInverse	* vec4 (0 ,0 , -1 ,0);

	vec3 I = normalize (V.xyz - E.xyz);

	vec3 N = normalize (gl_NormalMatrix * gl_Normal);

	Ratio = F + (1.0 - F) * pow((1.0 - dot(-I, N)), fresnelPower);

	Refract = refract(I, N, eta);
	Refract.y = -Refract.y;

	Refract.x = -Refract.x;
	Refract = vec3(gl_TextureMatrix[0] * vec4(Refract, 1.0));
	Reflect = reflect(I, N);
	Reflect = vec3(gl_TextureMatrix[0] * vec4(Reflect, 1.0));
	gl_Position = ftransform();
}