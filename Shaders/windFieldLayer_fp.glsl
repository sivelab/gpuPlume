//uniform sampler2D pointsprite_texunit;
//uniform sampler2D pointspritenormal_texunit;
//uniform int point_visuals;
uniform samplerRect Wind;
uniform float max_x;
uniform float max_y;
uniform float max_z;
uniform float max_c;
uniform float min_x;
uniform float min_y;
uniform float min_z;
uniform float min_c;
uniform int controlWind;
uniform float slider;

varying vec4 pcolor;

void main(void)
{
  vec2 texCoord = gl_TexCoord[0].xy;
  vec4 color = vec4(textureRect(Wind, texCoord)); 
  
  vec3 c1 = vec3(0.0,0.0,1.0);
  vec3 c2 = vec3(1.0,1.0,0.0);
  vec3 c3 = vec3(1.0,0.0,0.0);

  float t;

  float max = max_x;
  if(max_y>max)
    max = max_y;
  if(max_z>max)
    max = max_z;
  if(max_c>max)
    max = max_c;
  float min = min_x;
  if(min_y>min)
    min = min_y;
  if(min_z>min)
    min = min_z;
  if(min_c>min)
    min = min_c;

  float range = max - min;
  //float r = ((max-min)/2.0)/(max-min);
  //r = 0.05;

  if(controlWind == 0){
    //t = (color.x-min11)/(max11-min11);
    t = (color.x-min)/range;
  }
  //Show only Tau22
  else if(controlWind == 1){
    //t = (color.y-min22)/(max22-min22);
    t = (color.y-min)/range;
  }
  else if(controlWind == 2){
    //t = (color.z-min33)/(max33-min33);
    t = (color.z-min)/range;
  }
  else if(controlWind == 3){
    //t = (color.w-min13)/(max13-min13);
    t = (color.w-min)/range;
  }

  
  vec3 col;
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
  //gl_FragColor = vec4((c1*(1-t) + c2*t),0.8);

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
