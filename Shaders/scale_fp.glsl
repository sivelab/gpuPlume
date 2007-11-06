uniform float xmin;
uniform float xmax;
uniform float tauMin;
uniform float tauMax;
varying vec4 position;

void main(void)
{
  vec4 c1 = vec4(0.0,0.0,1.0,1.0);
  vec4 c2 = vec4(1.0,1.0,0.0,1.0);
  vec4 c3 = vec4(0.0,1.0,0.0,1.0);

  float t = (position.x - xmin)/(xmax-xmin);
  
  float r = ((xmax-xmin)/2.0)/(xmax-xmin);

  gl_FragColor = c1*(1.0-t) + c2*t;
  /*if(t >= r){
    gl_FragColor = c1*(1.0-(t-r)) + c2*t;
  }
  else
    gl_FragColor = vec4(0.0,0.0,0.0,1.0);
    //gl_FragColor = c3*(1.0-t) + c1*t;
    */

}
