//
//uniform float eta; // eta ratio
//
//varying vec3 R; // refract vector
//
//void main () {
	//// Create incident and normal vectors
	//vec4 V = gl_ModelViewMatrix * gl_Vertex ;
	//vec4 E = gl_ProjectionMatrixInverse
	//* vec4 (0 ,0 , -1 ,0);
	//vec3 I = normalize (V.xyz - E.xyz);
	//vec3 N = normalize ( gl_Normal );
	//R = refract (I, N, eta);
	//gl_Position = ftransform ();
//}
//

varying vec3 ReflectDir;
varying float LightIntensity;
uniform vec3 LightPos;
void main()
{
gl_Position = ftransform();
vec3 normal = normalize(gl_NormalMatrix * gl_Normal);
vec4 pos = gl_ModelViewMatrix * gl_Vertex;
vec3 eyeDir = pos.xyz;
ReflectDir = reflect(eyeDir, normal);
LightIntensity = max(dot(normalize(gl_LightSource[0].position.xyz - eyeDir), normal),0.0);
}
