varying vec4 position;

void main()
{
  position = vec4(gl_Vertex);
  
  gl_Position = gl_ModelViewProjectionMatrix*gl_Vertex;
  //gl_Position = ftransform();

}
