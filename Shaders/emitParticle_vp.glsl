varying vec4 start_pos;

void main(void)
{
	start_pos = gl_Color;
	gl_Position = ftransform();
}