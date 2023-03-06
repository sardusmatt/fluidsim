varying vec3 N;
varying vec3 View;
varying vec3 ScreenPos;
varying vec3 ecPosition3;


void main() {

	ecPosition3 = vec3(gl_ModelViewMatrix * gl_Vertex);
	View = normalize(-ecPosition3);
	N = normalize( gl_NormalMatrix * gl_Normal);
	gl_TexCoord[0] = gl_TextureMatrix[0] * vec4(gl_Normal,1);
	gl_Position = ftransform();
	ScreenPos = vec3(gl_Position);
}

