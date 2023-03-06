varying vec3 v_N;

uniform float k; // minnaert roughness

// Minnaert limb darkening diffuse term
vec3 minnaert( vec3 L, vec3 Nf, float k) {
	float ndotl = max( 0.0, dot(L, Nf));
	return gl_LightSource[0].diffuse.rgb * pow( ndotl, k);
}

void main() {
	vec3 N = normalize(v_N);
	vec3 L = normalize(vec3(gl_LightSouce[0].position));

	vec3 ambient = gl_FrontMaterial.ambient.rgb;
	vec3 diffuse = gl_FrontMaterial.diffuse.rgb * minnaert( L, N, k);
	vec3 specular = gl_FrontMaterial.specular.rgb;

	gl_FragColor = vec4( ambient + diffuse + specular, 1.0);
}

