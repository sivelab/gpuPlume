#version 120
#extension GL_EXT_geometry_shader4 : enable
uniform sampler3D tau;
uniform sampler1D case_to_numpoly;
uniform sampler2D edge_connect_list;

uniform float dx;
uniform float dy;
uniform float dz;

uniform float mesh;

float v0;
float v1;
float v2;
float v3;
float v4;
float v5;
float v6;
float v7;
vec3 cubePoints0;
vec3 cubePoints1;
vec3 cubePoints2;
vec3 cubePoints3;
vec3 cubePoints4;
vec3 cubePoints5;
vec3 cubePoints6;
vec3 cubePoints7;

vec3 getPointFromEdge(in int edge,in int c)
{
  float t;
  vec3 p = vec3(edge,50.0,20.0);
  if(c == 0)
    p = vec3(edge,50.0,20.0);
  else if(c == 1)
    p = vec3(100.0,edge,20.0);
  else if(c == 2)
    p = vec3(100.0,50.0,edge);

  //Position of point = direction*percentage(t) + first point location
  if(edge == 0){
    t = abs(v0)/(abs(v0)+abs(v1));   
    p = (cubePoints1-cubePoints0)*t + cubePoints0;   
  }
  else if(edge == 1){
    t = abs(v1)/(abs(v1)+abs(v2));   
    p = (cubePoints2-cubePoints1)*t + cubePoints1;   
    //p1 = 1;
    //p2 = 2;
  }
  else if(edge == 2){
    t = abs(v2)/(abs(v2)+abs(v3));   
    p = (cubePoints3-cubePoints2)*t + cubePoints2;   
    //p1 = 2;
    //p2 = 3;
  }
  else if(edge == 3){
    t = abs(v0)/(abs(v0)+abs(v3));   
    p = (cubePoints3-cubePoints0)*t + cubePoints0;    
    //p1 = 0;
    //p2 = 3;
  }
  else if(edge == 4){
    t = abs(v4)/(abs(v4)+abs(v5));   
    p = (cubePoints5-cubePoints4)*t + cubePoints4;   
    //p1 = 4;
    //p2 = 5;
  }
  else if(edge == 5){
    t = abs(v5)/(abs(v5)+abs(v6));   
    p = (cubePoints6-cubePoints5)*t + cubePoints5;   
    //p1 = 5;
    //p2 = 6;
  }
  else if(edge == 6){
    t = abs(v6)/(abs(v6)+abs(v7));   
    p = (cubePoints7-cubePoints6)*t + cubePoints6;   
    //p1 = 6;
    //p2 = 7;
  }
  else if(edge == 7){
    t = abs(v4)/(abs(v4)+abs(v7));   
    p = (cubePoints7-cubePoints4)*t + cubePoints4;   
    //p1 = 4;
    //p2 = 7;
  }
  else if(edge == 8){
    t = abs(v0)/(abs(v0)+abs(v4));   
    p = (cubePoints4-cubePoints0)*t + cubePoints0;   
    //p1 = 0;
    //p2 = 4;
  }
  else if(edge == 9){
    t = abs(v1)/(abs(v1)+abs(v5));   
    p = (cubePoints5-cubePoints1)*t + cubePoints1;   
    //p1 = 1;
    //p2 = 4;
  }
  else if(edge == 10){
    t = abs(v2)/(abs(v2)+abs(v6));   
    p = (cubePoints6-cubePoints2)*t + cubePoints2;   
    //p1 = 2;
    //p2 = 6;
  }
  else if(edge == 11){
    t = abs(v3)/(abs(v3)+abs(v7));   
    p = (cubePoints7-cubePoints3)*t + cubePoints3;   
    //p1 = 3;
    //p2 = 7;
  }

  return p;
}


void main(void)
{
   
  //First find the density function values for 
  //each corner of the cell
  vec3 texCoord;
  vec4 point;

  for(int i=0; i < 1; i++){
    point = gl_PositionIn[i];
  }


  texCoord.x = point.x/100.0;
  texCoord.y = point.y/50.0;
  texCoord.z = point.z/20.0;
  vec4 value = vec4(texture3D(tau,texCoord));
  v0 = value.x;
  cubePoints0 = vec3(point.x,point.y,point.z);
  
  texCoord.z += dz;
  value = vec4(texture3D(tau,texCoord));
  v1 = value.x;
  cubePoints1 = vec3(point.x,point.y,point.z+mesh);

  texCoord.x += dx;
  value = vec4(texture3D(tau,texCoord));
  v2 = value.x;
  cubePoints2 = vec3(point.x+mesh,point.y,point.z+mesh);

  texCoord.z -= dz;
  value = vec4(texture3D(tau,texCoord));
  v3 = value.x;
  cubePoints3 = vec3(point.x+mesh,point.y,point.z);

  texCoord.x -= dx;
  texCoord.y += dy;
  value = vec4(texture3D(tau,texCoord));
  v4 = value.x;
  cubePoints4 = vec3(point.x,point.y+mesh,point.z);

  texCoord.z += dz;
  value = vec4(texture3D(tau,texCoord));
  v5 = value.x;
  cubePoints5 = vec3(point.x,point.y+mesh,point.z+mesh);

  texCoord.x += dx;
  value = vec4(texture3D(tau,texCoord));
  v6 = value.x;
  cubePoints6 = vec3(point.x+mesh,point.y+mesh,point.z+mesh);

  texCoord.z -= dz;
  value = vec4(texture3D(tau,texCoord));
  v7 = value.x;
  cubePoints7 = vec3(point.x+mesh,point.y+mesh,point.z);

  // Now find the case number
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
  vec2 index;
  index.t = float(caseNum)/255.0;

  vec4 c = vec4(texture1D(case_to_numpoly,caseNum/255.0));
  int num_poly = int(c.x);
  
  //Now perform num_poly lookups into edge_connect_list
  
  float t;
  int i;
  vec3 edge_list;

  if(num_poly != 0 && (caseNum > 0) && (caseNum < 255)){

    for(i=0; i < 5; i++){
     
      if( i < num_poly ){

	index.s = float(i)/4.0;
	edge_list = vec3(texture2D(edge_connect_list,index));

	//for each set of edges, create a triangle
	//Create first point
         
	//Position of point = direction*percentage(t) + first point location
	vec3 pointPlace = getPointFromEdge(int(edge_list.x),0);  
   	gl_Position = vec4(pointPlace,1.0);
	EmitVertex();

      
	//Create second point
	pointPlace = getPointFromEdge(int(edge_list.y),1); 
	gl_Position = vec4(pointPlace,1.0);
	EmitVertex();

	//Create third point
	pointPlace = getPointFromEdge(int(edge_list.z),2);
	gl_Position = vec4(pointPlace,1.0);
	EmitVertex();
	EndPrimitive();

      }

    }
  
  }
 
}

