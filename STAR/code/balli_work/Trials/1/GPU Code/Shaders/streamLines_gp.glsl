//#version 120
//#extension GL_EXT_geometry_shader4 : enable

uniform sampler3D tau;
uniform sampler1D case_to_numpoly;
uniform sampler2D edge_connect_list;

void main()
{
   
  //First find the density function values for 
  //each corner of the cell
  vec3 texCoord;
  
  for(int i=0; i < gl_VerticesIn; i++){
    vec4 point = gl_PositionIn[i];
  }

  texCoord.x = point.x;
  texCoord.y = point.y;
  texCoord.z = point.z;
  vec4 value = vec4(texture3D(tau,texCoord));
  float v0 = value.x;

  texCoord.z += 1.0;
  value = vec4(texture3D(tau,texCoord));
  float v1 = value.x;

  texCoord.x += 1.0;
  value = vec4(texture3D(tau,texCoord));
  float v2 = value.x;

  texCoord.z -= 1.0;
  value = vec4(texture3D(tau,texCoord));
  float v3 = value.x;

  texCoord.x -= 1.0;
  texCoord.y += 1.0;
  value = vec4(texture3D(tau,texCoord));
  float v4 = value.x;

  texCoord.z += 1.0;
  value = vec4(texture3D(tau,texCoord));
  float v5 = value.x;

  texCoord.x += 1.0;
  value = vec4(texture3D(tau,texCoord));
  float v6 = value.x;

  texCoord.z -= 1.0;
  value = vec4(texture3D(tau,texCoord));
  float v7 = value.x;

  //Now find the case number
  int caseNum = 0;
  if(v0 >= 0.0){
    caseNum += 1;
  }
  if(v1 >= 0.0){
    caseNum += 2;
  }
  if(v2 >= 0.0){
    caseNum += 4;
  }
  if(v3 >= 0.0){
    caseNum += 8;
  }
  if(v4 >= 0.0){
    caseNum += 16;
  }
  if(v5 >= 0.0){
    caseNum += 32;
  }
  if(v6 >= 0.0){
    caseNum += 64;
  }
  if(v7 >= 0.0){
    caseNum += 128;
  }
  //Do a lookup into case_to_numpoly using the caseNum
  


  /*
  for(int i=0; i < gl_VerticesIn; i++)
  {
    //gl_FrontColor = gl_FrontColorIn[i];
    gl_FrontColor = vec4(0.0,0.0,1.0,1.0);
    vec4 point = gl_PositionIn[i];
    //gl_Position = gl_ModelViewProjectionMatrix*(point);
    gl_Position = point;
        
    EmitVertex();
       
    point.y = point.y + 0.5;
    gl_Position = point;
    // gl_Position.x = point.x + 1;
    ///gl_Position = gl_ModelViewProjectionMatrix*vec4(0.0,0.0,3.0,1.0);
    EmitVertex();
    

    point.z = point.z +0.5;
    gl_Position = point;
    EmitVertex();
    EndPrimitive();
    }*/
 
 
}
