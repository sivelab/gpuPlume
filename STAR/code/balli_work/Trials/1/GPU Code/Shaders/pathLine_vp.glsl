varying vec4 pos;

void main(void)
{
	//gl_TexCoord[0] = gl_MultiTexCoord0;
	pos = gl_Color;
	gl_Position = ftransform();
}