uniform int slice;
uniform sampler3D tau;


void main(void)
{

  vec4 color = vec4(1.0,0.0,0.0,1.0);
  /*
  if(slice == 0){
    color = vec4(0.0,1.0,0.0,1.0);
  }
  if(slice == 1){
    color = vec4(0.0,0.0,1.0,1.0);
  }
  if(slice == 2){
    color = vec4(1.0,0.0,1.0,1.0);
  }
  if(slice == 3){
    color = vec4(1.0,1.0,0.0,1.0);
    }	  
   if(slice == 19){
    color = vec4(0.0,1.0,1.0,1.0);
    }	   */
  vec3 texCoord = gl_TexCoord[0].xyz;
  vec4 value = vec4(texture3D(tau,texCoord));
		    
  //if(value.x > 1.58 && value.x < 1.59){
  //if(value.x <  1.0  ){
  //color = vec4(-1.0,0.0,0.0,1.0);
  //}
  
  //Density function

  //set contour value to 1
  //will need to pass in as uniform

  value.x = value.x - 1;

  color = value;

  gl_FragColor = color;

}
