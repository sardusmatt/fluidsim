varying vec3 v_V;
varying vec3 v_N;

uniform float k; // minnaert roughness  1.5
uniform float roughness; // Ward isotropic specular roughness 0.2
uniform float ior; // Schlick's fresnel approximation index of refraction 1.5

// Minnaert limb darkening diffuse term
vec3 minnaert( vec3 L, vec3 Nf, float k) {
	float ndotl = max( 0.0, dot(L, Nf));
	return gl_LightSource[0].diffuse.rgb * pow( ndotl, k);
}

// Ward isotropic specular term
vec3 wardiso( vec3 Nf, vec3 Ln, vec3 Hn, float roughness, float ndotv ) {
	float ndoth = dot( Nf, Hn);
	float ndotl = dot( Nf, Ln);
	float tandelta = tan( acos(ndoth));
	return gl_LightSource[0].specular.rgb
		* exp( -( pow( tandelta, 2.0) / pow( roughness, 2.0)))
		* (1.0 / sqrt( ndotl * ndotv ))
		* (1.0 / (4.0 * pow( roughness, 2.0)));
	}
	
float schlick( vec3 Nf, vec3 Vf, float ior, float ndotv ) {
	float kr = (ior-1.0)/(ior+1.0);
	kr *= kr;
	return kr + (1.0-kr)*pow( 1.0 - ndotv, 5.0);
}
	
void main() {
	vec3 N = normalize(v_N);
	vec3 V = normalize(v_V);
	vec3 L = normalize(vec3(gl_LightSource[0].position));
	vec3 Vf = -V;
	float ndotv = dot(N, Vf);
	vec3 H = normalize(L+Vf);

	vec3 ambient = gl_FrontMaterial.ambient.rgb;
	vec3 diffuse = gl_FrontMaterial.diffuse.rgb * minnaert( L, N, k);
	float fresnel = schlick( N, V, ior, ndotv);
	vec3 specular = gl_FrontMaterial.specular.rgb * wardiso( N, L, H, roughness, ndotv) * fresnel;

	gl_FragColor = vec4( ambient + diffuse + specular, 1.0);
}

