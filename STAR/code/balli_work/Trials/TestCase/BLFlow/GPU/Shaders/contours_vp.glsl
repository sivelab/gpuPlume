//varying vec4 position;

void main(void)
{
  //vec4 pos = vec4(gl_Vertex);
  //position = gl_ModelViewMatrix * pos;
  gl_TexCoord[0] = gl_MultiTexCoord0;
  gl_Position = ftransform();
	
}
