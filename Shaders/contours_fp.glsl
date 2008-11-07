uniform int numContours;
uniform vec4 tauValue;
//uniform samplerRect tau;
uniform sampler2DRect contourTex;
uniform sampler3D tau3d;
uniform int height;

void main(void)
{
  vec4 color1 = vec4(0.0,0.0,1.0,1.0);
  vec4 color2 = vec4(1.0,1.0,0.0,1.0);
  vec4 color3 = vec4(1.0,0.0,0.0,1.0);

  vec3 texCoord = gl_TexCoord[0].xyz;
  vec4 value = vec4(texture3D(tau3d, texCoord));
 
  vec4 color = vec4(1.0,1.0,1.0,1.0);
 
  //performs a dot product to get the current tau value
  float t = dot(value,tauValue);
  
  float c2;
 
  vec2 index;
  int scale = -1;

  int i=0;
  index.s = float(height*numContours);
  index.t = 0.0;
 
  float c1 =  dot(vec4(texture2DRect(contourTex,index)),tauValue);
  
  //Colors the smallest contour value
  if(t <= c1)
    scale = 0;
        
  //Colors all the tau values in between 1st contour line and last contour line
  for(i=0;i<(numContours-1); i++){
    index.s = float(height*numContours + i);
    index.t = 0.0;
    c1 =  dot(vec4(texture2DRect(contourTex,index)),tauValue);

    index.s = float(height*numContours + i + 1);
    c2 =  dot(vec4(texture2DRect(contourTex,index)),tauValue);

    if((c1 < t) && ( t <= c2)){
      scale = i+1;
    }
   
  }
  //Color the largest contour value
  index.s = float(height*numContours + numContours-1);
  c1 =  dot(vec4(texture2DRect(contourTex,index)),tauValue);

  if(t > c1){
    scale = numContours;
  }
  
  

  //The color is based on which contour region the value lies in.
  //The variable scale determines the region 
    
  float s = float(scale)/float(numContours);
  float n = 1.0/2.0;

  
  
  if(s >= n){
    s = (s-n)/(1.0-n);
    
    color = color1*(1.0-s) + color2*s;
  }
  else{
    s = (n-s)/n;
    
    color = color3*s + color1*(1.0-s);
  }
  gl_FragColor = color;
 
 

}
