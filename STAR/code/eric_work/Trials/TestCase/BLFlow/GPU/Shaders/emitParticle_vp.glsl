varying vec4 start_pos;
varying vec3 init_prime;

void main(void)
{
	start_pos = gl_Color;
	init_prime = gl_Normal;
	gl_Position = ftransform();
}