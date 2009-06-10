
uniform sampler3D Tau;
uniform float max11;
uniform float max22;
uniform float max33;
uniform float max13;
uniform float min11;
uniform float min22;
uniform float min33;
uniform float min13;
uniform int controlTau;
uniform float slider;

void main(void)
{
  
  vec3 texCoord = gl_TexCoord[0].xyz;
 
  vec4 color = vec4(texture3D(Tau, texCoord));

    
  vec3 c1 = vec3(0.0,0.0,1.0);
  vec3 c2 = vec3(1.0,1.0,0.0);
  vec3 c3 = vec3(1.0,0.0,0.0);

  float t;

  float max = max11;
  if(max22>max)
    max = max22;
  if(max33>max)
    max = max33;
  if(max13>max)
    max = max13;
  float min = min11;
  if(min22<min)
    min = min22;
  if(min33<min)
    min = min33;
  if(min13<min)
    min = min13;

  vec3 col;
  float range = max - min;
  //float r = ((max-min)/2.0)/(max-min);
  //r = 0.05;

  if(controlTau == 1){
    //t = (color.x-min11)/(max11-min11);
    t = (color.x-min)/range; 
    
  }
  //Show only Tau22
  else if(controlTau == 2){
    //t = (color.y-min22)/(max22-min22);
    t = (color.y-min)/range;
  }
  else if(controlTau == 3){
    //t = (color.z-min33)/(max33-min33);
    t = (color.z-min)/range;
  }
  else if(controlTau == 4){
    //t = (color.w-min13)/(max13-min13);
    t = (color.w-min)/range;
  }

  
 
  if(t >= slider){
    t = (t-slider)/(1.0-slider);
    col = c1*(1.0-t) + c2*t;
  }
  else{
    t = (slider-t)/(slider);   
    col = c3*(t) + c1*(1.0-t);
    //gl_FragColor = c3*(1.0-t) + c1*(t);
  }
  
  gl_FragColor = vec4(col,0.8);
 
    
}
