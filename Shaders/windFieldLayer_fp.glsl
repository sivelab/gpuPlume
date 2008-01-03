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
  vec4 wind_dir = vec4(textureRect(Wind, texCoord)); 
  
  mat3 opponent2rgb = mat3(1.0,  0.1140,  0.7436, 
			   1.0,  0.1140, -0.4111, 
			   1.0, -0.8860, 0.1663);
  vec3 opponent_color = normalize(wind_dir);


  // opponent_color = normalize(wind);

  vec3 color = opponent2rgb * opponent_color;

  // now, we need to apply a non-affine transform to make the red/green axis perpendicular to the blue/yellow axis
  float pi = 3.14159;
  float theta = atan(color.z, color.y);
  float theta_0 = 0.0;
  if (theta < (pi/3.0))
    theta_0 = (3.0 / 2.0) * theta;
  else if (theta <= pi && theta >= (pi/3.0))
    theta_0 = pi/2.0 + 0.75*(theta - pi/3.0);
      
  vec2 oRGB_new, oRGB;
  oRGB = vec2(color.y, color.z);
  mat2 rot = mat2(cos(theta), -sin(theta), sin(theta), cos(theta));
  oRGB_new = rot * oRGB;
      
  color.y = oRGB_new.x;  
  color.z = oRGB_new.y;

  gl_FragColor = vec4(color,1.0);

  
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
