uniform samplerRect primePrev;
uniform samplerRect wind;
uniform samplerRect random;
uniform samplerRect lambda;
uniform float time_step;

void main(void)
{
	vec2 texCoord = gl_TexCoord[0].xy;
	vec3 prime = vec3(textureRect(primePrev, texCoord));

	//Currently this code just keeps passing the prime values,
	//keeping them the same.  
	//We need to add the equations to calculate the new values.
	gl_FragColor = vec4(prime, 1.0);


}