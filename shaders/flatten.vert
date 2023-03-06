// Flattens then model by zeroing the z coordinate
// gl_Vertex is an attribute, so it cannot be modified in place

void main(void) {
	vec4 v = vec4(gl_Vertex);
	v.y = 0.0;

	gl_Position = gl_ModelViewProjectionMatrix * v;
}
