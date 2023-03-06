varying vec3 N;
varying vec3 View;
uniform sampler2D Texture;
uniform samplerCube Environment;
float refraction_index = 1.9;

void myRefract(in vec3 incom, in vec3 normal, in float index_external, in float index_internal,
				out vec3 reflection, out vec3 refraction,
				out float reflectance, out float transmittance) 
{
				
	float eta = index_external/index_internal;
	float cos_theta1 = dot(incom, normal);
	float cos_theta2 = sqrt(1.0 - ((eta * eta) * ( 1.0 - (cos_theta1 * cos_theta1))));
	reflection = incom - 2.0 * cos_theta1 * normal;
	refraction = (eta * incom) + (cos_theta2 - eta * cos_theta1) * normal;
	refraction.x = -refraction.x;
	refraction.y = -refraction.y;

	float fresnel_rs = (index_external * cos_theta1 - index_internal * cos_theta2 ) /
		(index_external * cos_theta1 + index_internal * cos_theta2);


	float fresnel_rp = (index_internal * cos_theta1 - index_external * cos_theta2 ) /
			(index_internal * cos_theta1 + index_external * cos_theta2);

	reflectance =  (fresnel_rs * fresnel_rs + fresnel_rp * fresnel_rp) / 2.0;
	transmittance =  ((1.0-fresnel_rs) * (1.0-fresnel_rs) + (1.0-fresnel_rp) * (1.0-fresnel_rp)) / 2.0;
}

void main() {


	
	vec3 nN = 0.5 * (N + 2.0*(texture2D(Texture, gl_TexCoord[0].xy).rgb - 0.5));
	
	vec3 refraction_ray, reflection_ray;
	float fresnel_R, fresnel_T;
	
	myRefract(View, nN, 1.0, refraction_index,
		reflection_ray, refraction_ray, fresnel_R, fresnel_T);

	refraction_ray = -(gl_ModelViewMatrixTranspose * vec4(refraction_ray,0.0)).xyz;
	reflection_ray = -(gl_ModelViewMatrixTranspose * vec4(reflection_ray,0.0)).xyz;

	vec4 reflect_color = textureCube(Environment, reflection_ray);
	vec4 refract_color = textureCube(Environment, refraction_ray);
	
	fresnel_T = fresnel_T * 0.5;	
	fresnel_R = fresnel_R * 0.5;	

	gl_FragColor =
		  reflect_color * fresnel_R +
		  refract_color * fresnel_T;

}
