uniform float xmin;
uniform float xmax;
uniform float tauMin;
uniform float tauMax;
uniform float slider;
varying vec4 position;

void main(void)
{
  vec4 c1 = vec4(0.0,0.0,1.0,1.0);
  vec4 c2 = vec4(1.0,1.0,0.0,1.0);
  vec4 c3 = vec4(1.0,0.0,0.0,1.0);

  float t = (position.x - xmin)/(xmax-xmin);
  
  //float r = ((xmax-xmin)/2.0)/(xmax-xmin);
  //gl_FragColor = c1*(1.0-t) + c2*t;
  //r = 0.05;

  if(t >= slider){
    t = (t-slider)/(1.0-slider);
    gl_FragColor = c1*(1.0-t) + c2*t;
  }
  else{
    t = (slider-t)/(slider);
    
    gl_FragColor = c3*(t) + c1*(1.0-t);
    //gl_FragColor = c3*(1.0-t) + c1*(t);
  }

}
