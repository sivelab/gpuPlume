void main(void)
{
	vec4 pos = vec4(gl_Vertex);
	pos.w = 1.0;
	gl_Color = vec4(1.0, 1.0, 1.0, 1.0);
	gl_Position = gl_ModelViewProjectionMatrix * pos;
}