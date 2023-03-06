//varying vec3 lightDir;                                          
//varying vec3 viewDir;
//varying float intensity;
//
	//void main()
	//{
		//
		//lightDir = normalize(vec3(gl_LightSource[0].position));
	//
		//intensity = dot(lightDir,gl_Normal);
//
		//gl_Position = ftransform();
	//} 
//
//void main() {
	//
	//gl_Position = ftransform(); //gl_ModelViewProjectionMatrix * gl_Vertex;
//}


//varying vec3 texCoord;

void main() {
	//texCoord = gl_Normal.xyz;
	gl_Position = ftransform();		
}