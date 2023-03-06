//uniform samplerCube CubeMap;
//
//varying vec3 R; // refracted vector
//
//void main () { 
	//gl_FragColor = textureCube (CubeMap , R);
//}

vec3 BaseColor = vec3(1.0,1.0,1.0);
float MixRatio = 0;
uniform samplerCube EnvMap;
varying vec3 ReflectDir;
varying float LightIntensity;
void main()
{
// Look up environment map value in cube map
vec3 envColor = vec3(textureCube(EnvMap, ReflectDir));
// Add lighting to base color and mix
vec3 base = LightIntensity * BaseColor;
envColor = mix(envColor, base, MixRatio);
gl_FragColor = vec4(envColor, 1.0);
}
