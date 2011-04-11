
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
    
  vec3 layerColor; 
  float tke = (color.x + color.y + color.z) * 0.5;

  float minTKE = (min11 + min22 + min33) * 0.5;
  float maxTKE = (max11 + max22 + max33) * 0.5;
  float rangeTKE = maxTKE - minTKE;

  tke = (tke - minTKE) / rangeTKE;

  vec3 DKRED = vec3(0.5, 0.0, 0.0);
  vec3 RED = vec3(1.0, 0.0, 0.0);
  vec3 YELLOW = vec3(1.0, 1.0, 0.0);
  vec3 CYAN = vec3(0.0, 1.0, 1.0);
  vec3 BLUE = vec3(0.0, 0.0, 1.0);
  vec3 DKBLUE = vec3(0.0, 0.0, 0.5);

  // Mapping the following
  // Dark Blue <--> Blue <--> Cyan <--> Yellow <--> Red <--> Dark Red
  //     0          0.125     0.375     0.625       0.875      1.0  
  float t = 0.0;
  if (tke <= 0.125)
    {
      t = tke * 8.0;  // scale to 0-1
      layerColor = mix(DKBLUE, BLUE, t);
    }    
  else if (tke <= 0.375) 
    { 
      t = (tke - 0.125) * 4.0;
      layerColor = mix(BLUE, CYAN, t);
    }
  else if (tke <= 0.625) 
    { 
      t = (tke - 0.375) * 4.0;
      layerColor = mix(CYAN, YELLOW, t);
    }
  else if (tke <= 0.875)
    { 
      t = (tke - 0.625) * 4.0;
      layerColor = mix(YELLOW, RED, t);
    }
  else
    {	
      t = (tke - 0.875) * 8.0;
      layerColor = mix(RED, DKRED, t);
    }

  gl_FragColor = vec4(layerColor,0.75);
}
