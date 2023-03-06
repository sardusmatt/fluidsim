// Toon shader third ver (gl lights are used...gl lights are not required to be enabled)

varying vec3 normal;


void main()
{
	normal = gl_NormalMatrix * gl_Normal;
	gl_Position = ftransform();

}
