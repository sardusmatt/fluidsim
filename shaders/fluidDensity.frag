//uniform sampler3D g_densityTex;
//uniform vec3 g_lightPos;
//uniform vec3 g_eyePos;
//
//varying vec3 lightDir;                                         
//varying vec3 viewDir;
//uniform float g_absorption;

uniform sampler3D g_densityTex;
//uniform vec3 g_lightPos;
//uniform vec3 g_lightIntensity;
//uniform vec3 g_eyePos;
//uniform float g_absorption;
//
void main()
{
	// HACKS
	const float intensity = 1.0f;
	const vec3 g_lightIntensity = vec3(intensity, intensity, intensity);
	vec3 g_lightPos = gl_LightSource[0].position.xyz;
	vec3 g_eyePos = vec3(0.0f, 0.0f, 0.0f);
	const float g_absorption = 0.5f;

	/////

	// diagonal of the cube
	const float maxDist = sqrt(3.0);

	const int numSamples = 128;
	const float scale = maxDist/float(numSamples);

	const int numLightSamples = 32;
	const float lscale = maxDist / float(numLightSamples);

	// assume all coordinates are in texture space
	vec3 pos = gl_TexCoord[0].xyz;
	vec3 eyeDir = normalize(pos-g_eyePos)*scale;

	// transmittance
	float T = 1.0;
	// in-scattered radiance
	vec3 Lo = vec3(0.0);

	for (int i=0; i < numSamples; ++i)
	{
		// sample density
		float density = texture3D(g_densityTex, pos).x;

		// skip empty space
		if (density > 0.0)
		{
			// attenuate ray-throughput
			T *= 1.0-density*scale*g_absorption;
			if (T <= 0.01)
				break;

			// point light dir in texture space
			vec3 lightDir = normalize(g_lightPos-pos)*lscale;

			// sample light
			float Tl = 1.0;	// transmittance along light ray
			vec3 lpos = pos + lightDir;

			for (int s=0; s < numLightSamples; ++s)
			{
				float ld = texture3D(g_densityTex, lpos).x;
				Tl *= 1.0-g_absorption*lscale*ld;

				if (Tl <= 0.01)
					break;

				lpos += lightDir;
			}

			vec3 Li = g_lightIntensity*Tl;

			Lo += Li*T*density*scale;
		}

		pos += eyeDir;
	}

	gl_FragColor.xyz = Lo;
	gl_FragColor.w = 1.0-T;
}

