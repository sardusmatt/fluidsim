// Toon shader first ver (light direction is computed in the application code and then
// passed in as a uniform

// The vertex shader has access to the normals, as specified in the OpenGL application, 
// through the attribute variable gl_Normal. This is the normal as defined in the OpenGL
// application with the glNormal function, hence in model local space. ASSUMING THAT
// the light direction is defined in world space and no rotations and scale transforms
// are used, then the normal is the same in normal space. Thus, the computation does not
// need to perform any transformation on these directions before using them

// We expect lightDir and normals to be normalised (that is, magnitude = 1), so the
// intensity can be computed as the dot product lightDir . normal

// intensity is declared as varying both in the vertex and fragment shader, so that they
// can work on the same variable (the vertex shader writes it, the fragment shader reads
// the value

uniform vec3 lightDir;
varying float intensity;

void main() {
	intensity = dot(lightDir,gl_Normal);
	gl_Position = ftransform();
}




