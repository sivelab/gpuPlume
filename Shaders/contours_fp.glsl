uniform int numContours;
uniform int tauValue;
//uniform samplerRect tau;
uniform samplerRect contourTex;
uniform sampler3D tau3d;
uniform int height;

void main(void)
{
  vec4 color1 = vec4(0.0,0.0,1.0,1.0);
  vec4 color2 = vec4(1.0,1.0,0.0,1.0);
  vec4 color3 = vec4(1.0,0.0,0.0,1.0);

  vec3 texCoord = gl_TexCoord[0].xyz;
  vec4 value = vec4(texture3D(tau3d, texCoord));

 
  float t;
  vec4 color = vec4(1.0,1.0,1.0,1.0);

  if(tauValue == 0){
    t = value.x;
  }
  else if(tauValue == 1)
    t = value.y;
  else if(tauValue == 2)
    t = value.z;
  else t = value.w;

  
  vec4 c1;
  vec4 c2;
  
  vec2 index;
  int scale = 0;

  int i=0;
  index.s = height*numContours;
  index.t = 0;

  c1 = vec4(textureRect(contourTex,index));
  
  //Tau11
  if(tauValue == 0){

    //Colors the smallest contour value
    if(t < c1.x)
      color = vec4(1.0,1.0,1.0,1.0);

    //if(t  < (c1.x + 0.0001) && t > (c1.x - 0.0001))
    //color = vec4(-1.0,0.0,0.0,0.0);
    

    //Colors all the tau values in between 1st contour line and last contour line
    for(i=0;i<(numContours-1); i++){
      index.s = height*numContours + i;
      index.t = 0;
      c1 = vec4(textureRect(contourTex,index));
		     
      index.s = height*numContours + i + 1;
      c2 = vec4(textureRect(contourTex,index));

      if((c1.x < t) && ( t < c2.x)){
	color = vec4(0.0,0.0,1.0,1.0);
	scale = i+1;
      }
   
    }
    //Color the largest contour value
    index.s = height*numContours + numContours-1;
    c1 = vec4(textureRect(contourTex,index));
    if(t > c1.x){
      color = vec4(0.0,0.0,0.0,1.0);
      scale = numContours;
    }
  }
  //Tau22
  else if(tauValue == 1){
    //Colors the smallest contour value

    if(t < c1.y)
      color = vec4(1.0,0.0,0.0,1.0);

    //Colors all the tau values in between 1st contour line and last contour line
    for(i=0;i<(numContours-1); i++){
      index.s = height*numContours + i;
      index.t = 0;
      c1 = vec4(textureRect(contourTex,index));
		     
      index.s = height*numContours + i + 1;
      c2 = vec4(textureRect(contourTex,index));

      if((c1.y < t) && ( t < c2.y)){
	color = vec4(0.0,0.0,1.0,1.0);
	scale = i+1;
      }
   
    }

    //Color the largest contour value
    index.s = height*numContours + numContours-1;
    c1 = vec4(textureRect(contourTex,index));
    if(t > c1.y){
      color = vec4(0.0,0.0,0.0,1.0);
      scale = numContours;
    }
  }
  //Tau33
  else if(tauValue == 2){
    //Colors the smallest contour value
    if(t < c1.z)
      color = vec4(1.0,0.0,0.0,1.0);

    //Colors all the tau values in between 1st contour line and last contour line
    for(i=0;i<(numContours-1); i++){
      index.s = height*numContours + i;
      index.t = 0;
      c1 = vec4(textureRect(contourTex,index));
		     
      index.s = height*numContours + i + 1;
      c2 = vec4(textureRect(contourTex,index));

      if((c1.z < t) && ( t < c2.z)){
	color = vec4(0.0,0.0,1.0,1.0);
	scale = i+1;
      }
   
    }
    //Color the largest contour value
    index.s = height*numContours + numContours-1;
    c1 = vec4(textureRect(contourTex,index));
    if(t > c1.z){
      color = vec4(0.0,0.0,0.0,1.0);
      scale = numContours;
    }
  }
  //Tau13
  else{
    //Colors the smallest contour value
    if(t < c1.w)
      color = vec4(1.0,0.0,0.0,1.0);

    //Colors all the tau values in between 1st contour line and last contour line
    for(i=0;i<(numContours-1); i++){
      index.s = height*numContours + i;
      index.t = 0;
      c1 = vec4(textureRect(contourTex,index));
		     
      index.s = height*numContours + i + 1;
      c2 = vec4(textureRect(contourTex,index));

      if((c1.w < t) && ( t < c2.w)){
	color = vec4(0.0,0.0,1.0,1.0);
	scale = i+1;
      }
   
    }
    //Color the largest contour value
    index.s = height*numContours + numContours-1;
    c1 = vec4(textureRect(contourTex,index));
    if(t > c1.w){
      color = vec4(0.0,0.0,0.0,1.0);
      scale = numContours;
    }
  }
    

  float s = float(scale)/float(numContours);
  float n = 1.0/2.0;

  //if(color.x == -1.0){
  //gl_FragColor = vec4(1.0,1.0,1.0,1.0);
  //}
  //else{
    if(s >= n){
      s = (s-n)/(1.0-n);

      color = color1*(1.0-s) + color2*s;
    }
    else{
      s = (n-s)/n;

      color = color3*s + color1*(1.0-s);
    }
    gl_FragColor = color;
    //}
 

}
