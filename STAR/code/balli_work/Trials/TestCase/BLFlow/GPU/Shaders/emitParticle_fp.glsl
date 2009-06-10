varying vec4 start_pos;
varying vec3 init_prime;

void main(void)
{
        
	//gl_FragColor = start_pos;
	gl_FragData[0] = start_pos;
	gl_FragData[1] = vec4(init_prime,1.0);
}