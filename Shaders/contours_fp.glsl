uniform int numContours;
uniform int tauValue;
uniform float contourValues[];
uniform samplerRect tau;

void main(void)
{
  vec2 texCoord = gl_TexCoord[0].xy;
  vec4 value = vec4(textureRect(tau, texCoord));

  float t;
  vec4 color = vec4(0.0,0.0,0.0,1.0);

  if(tauValue == 0){
    t = value.x;
  }
  else if(tauValue == 1)
    t = value.y;
  else if(tauValue == 2)
    t = value.z;
  else t = value.w;

  

  
  int i=0;
  //Colors all the tau values in between 1st contour line and last contour line
  /*for(i=0;i<(numContours-1); i++){
    if( (contourValues[i] < t) && (t < contourValues[(i+1)]) ){
      color = vec4(0.0,0.0,1.0,1.0);
      break;
    }
    
    }*/
  /*while(i < (numContours-1)){
    if( (contourValues[i] < t) && (t < contourValues[i+1]) ){
      color = vec4(0.0,0.0,1.0,1.0);
    }
    
    i = i+1;
    }*/

  //Check 1st contour value
  i = (int)3;
  if( t < contourValues[0] ){
    color = vec4(0.0,1.0,0.0,1.0);
  }
  else if( t > contourValues[0] && t < contourValues[1]){
    color = vec4(0.0,0.0,1.0,1.0);
  }
  else if( t > contourValues[1] && t < contourValues[2]){
    color = vec4(1.0,1.0,0.0,1.0);
  }
  else if( t > contourValues[2] && t < contourValues[3]){
    color = vec4(1.0,0.0,1.0,1.0);
  }
  else if( t > contourValues[3] ){
    color = vec4(1.0,0.0,0.0,1.0);
  }

  /*
  if(numContours == 3)
    color.x = 1.0;
  if(numContours == 4)
    color.y = 1.0;
  if(numContours == 2)
  color.z = 1.0;*/
  
  gl_FragColor = color;

}
