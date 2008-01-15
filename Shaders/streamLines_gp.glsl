//#version 120
//#extension GL_EXT_geometry_shader4 : enable

void main()
{
  
  for(int i=0; i < gl_VerticesIn; i++)
  {
    gl_FrontColor = gl_FrontColorIn[i];
    gl_Position = gl_PositionIn[i];
    gl_FrontColor = vec4(0.0,0.0,1.0,1.0);
    EmitVertex();
    
  
    gl_Position.x = gl_Position.x + 1;
    ///gl_Position = gl_ModelViewProjectionMatrix*vec4(0.0,0.0,3.0,1.0);
    

    EmitVertex();
    
    gl_Position.y = gl_Position.y -1;
    //gl_Position = gl_ModelViewProjectionMatrix*vec4(

    EmitVertex();

    EndPrimitive();
  }
  //EndPrimitive();

  //vec4 v1 = gl_PositionIn[0];
  
  //gl_Position = gl_PositionIn[0];
  //gl_FrontColor = vec4(0.0,0,0.1,0,1.0);
  
  //EmitVertex();


}
