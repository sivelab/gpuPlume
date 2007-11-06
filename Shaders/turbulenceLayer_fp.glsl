//uniform sampler2D pointsprite_texunit;
//uniform sampler2D pointspritenormal_texunit;
//uniform int point_visuals;
uniform samplerRect Tau;
uniform float max11;
uniform float max22;
uniform float max33;
uniform float max13;
uniform float min11;
uniform float min22;
uniform float min33;
uniform float min13;
uniform int controlTau;

varying vec4 pcolor;

void main(void)
{
  vec2 texCoord = gl_TexCoord[0].xy;
  vec4 color = vec4(textureRect(Tau, texCoord)); 
  
  vec3 c1 = vec3(0.0,0.0,1.0);
  vec3 c2 = vec3(1.0,1.0,0.0);
  float t;

  float max = max11;
  if(max22>max)
    max = max22;
  if(max33>max)
    max = max33;
  if(max13>max)
    max = max13;
  float min = min11;
  if(min22>min)
    min = min22;
  if(min33>min)
    min = min33;
  if(min13>min)
    min = min13;

  float range = max - min;


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
    t = (color.z-min)/range;
  }

  gl_FragColor = vec4((c1*(1-t) + c2*t),0.8);

  /*
  //Show only Tau11
  if(controlTau == 1){
    color.x = (color.x-min11)/(max11-min11);
    if(color.x > 0.05)
      gl_FragColor = vec4(color.x,0.0,0.0,1.0);
    else
      gl_FragColor = vec4(color.x,0.0,0.0,0.7);
  }
  //Show only Tau22
  else if(controlTau == 2){
    color.y = (color.y-min22)/(max22-min22);
    if(color.y > 0.05)
      gl_FragColor = vec4(0.0,color.y,0.0,1.0);
    else
      gl_FragColor = vec4(0.0,color.y,0.0,0.7);
  }
  //Show only Tau33
  else if(controlTau == 3){
    color.z = (color.z-min33)/(max33-min33);
    if(color.z > 0.05)
      gl_FragColor = vec4(0.0,0.0,color.z,1.0);
    else
      gl_FragColor = vec4(0.0,0.0,color.z,0.7);
  }
  else if(controlTau == 4){
    color.w = (color.w-min13)/(max13-min13);
    if(color.w > 0.05)
      gl_FragColor = vec4(color.w,color.w,0.0,1.0);
    else
      gl_FragColor = vec4(color.w,color.w,0.0,0.7);
  }
  //Show all together
  else{
    
    float max = max11;
    if(max22>max)
      max = max22;
    if(max33>max)
      max = max33;
    if(max13>max)
      max = max13;
    float min = min11;
    if(min22>min)
      min = min22;
    if(min33>min)
      min = min33;
    if(min13>min)
      min = min13;

    float range = max - min;

    color.x = (color.x-min)/range;
    color.y = (color.y-min)/range;
    color.z = (color.z-min)/range;
    color.w = 1.0-((color.w-min)/range);


    gl_FragColor = color;
    }*/
  /*if (point_visuals == 1)
    {
    //vec4 d = vec4(texture2D(pointspritenormal_texunit, gl_TexCoord[0].st));
    vec4 normal = vec4(texture2D(pointspritenormal_texunit, gl_TexCoord[0].st));
    if (normal.w == 0.0){
    discard;
    }	
    //else if (normal.x == 0.0 && normal.y == 0.0 && normal.z == 1.0)
    //{
    // vec4 color = vec4(0.0, 0.0, 1.0, 1.0);
    //gl_FragColor = color;
    //}
    else
    {
    vec4 l = normalize(vec4(1.0, -1.0, 1.0, 0.0));
    vec4 n = (normal * 2.0) - vec4(1.0);
    n.w = 0.0;
    vec4 h = normalize(l + vec4(0.0,0.0,1.0,0.0));

    float p = 64.0;
    float cp = 1.0;
    vec4 cr = pcolor;
    vec4 cl = vec4(1.0,1.0,1.0,1.0);
    vec4 ca = vec4(0.4,0.4,0.4,1.0);

    vec4 color = cr * (ca + cl * max(0.0,dot(n,l))) + cp * cl * pow(max(0.0,dot(h,n)), p);  
    gl_FragColor = color;
    }
    }
    else 
    {
    // 
    // If using simple point based particles, this is all that needs to be done
    // 
    gl_FragColor = pcolor;
    }*/
}
