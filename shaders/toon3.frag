// The light direction is now used in the fragment to compute the correct color value based
// on the normal value for the pixel

varying vec3 normal;

uniform float lowerThreshold, middleThreshold, higherThreshold;

void main()
{
	float intensity;
	vec4 color;
	vec3 n = normalize(normal);
	intensity = dot(vec3(gl_LightSource[0].position),n);

	if (intensity > higherThreshold)
		color = vec4(0.5,0.5,1.0,1.0);
	else if (intensity > middleThreshold)
		color = vec4(0.25,0.25,0.5,1.0);
	else if (intensity > lowerThreshold)
		color = vec4(0.125,0.125,0.25,1.0);
	else
		color = vec4(0.2,0.1,0.1,1.0);
	gl_FragColor = color;

}