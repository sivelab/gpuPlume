//uniform sampler2D pointsprite_texunit;
//uniform sampler2D pointspritenormal_texunit;
//uniform int point_visuals;
uniform samplerRect Wind;

varying vec4 pcolor;

void main(void)
{
  vec2 texCoord = gl_TexCoord[0].xy;
  vec4 wind_dir = vec4(textureRect(Wind, texCoord)); 
  
  // mat3 opponent2rgb = mat3(1.0,  0.1140,  0.7436, 
  // 1.0,  0.1140, -0.4111, 
  // 1.0, -0.8860, 0.1663);
  mat3 opponent2rgb = mat3(0.299,  0.587,  0.114, 
			   0.500,  0.500, -1.000, 
			   0.866, -0.886,  0.000);
  vec3 opponent_color = normalize(wind_dir).zxy;

  // opponent_color = normalize(wind);

  opponent_color.x = (1.0 + opponent_color.x)/2.0;
  opponent_color.y = opponent_color.y * sqrt(2.0);
  opponent_color.z = opponent_color.z * sqrt(2.0);

  vec3 color = opponent2rgb * opponent_color;

      float pi = 3.141592654;
      float pi_3 = 1.047197551;
      float theta_0 = 0.0;

      vec2 oRGB_new, oRGB;
      float theta = atan(color.z, color.y);
      if (theta >= 0.0)
      	{
	  // value of theta will/should be in [0 to pi] range
	  if (theta < pi_3)
	    theta_0 = 1.5 * theta;
	  else if (theta <= pi && theta >= pi_3)
	    theta_0 = pi/4.0 + 0.75*theta;
	  //	theta_0 = pi/2.0 + 0.75*(theta - pi/3.0);
	  
	  oRGB = vec2(color.y, color.z);
	  mat2 rot = mat2(cos(theta_0 - theta), -sin(theta_0 - theta), sin(theta_0 - theta), cos(theta_0 - theta));
	  oRGB_new = rot * oRGB;
	}
      else 
	{
	  // value of theta will/should be in [0 to -pi] range
	  theta = -1.0 * theta;
	  if (theta < pi_3)
	    theta_0 = 1.5 * theta;
	  else if (theta <= pi && theta >= pi_3)
	    theta_0 = pi/4.0 + 0.75*theta;
	  //	theta_0 = pi/2.0 + 0.75*(theta - pi/3.0);
	  
	  oRGB = vec2(color.y, color.z);
	  theta_0 = -1.0 * theta_0;
	  theta = -1.0 * theta;
	  mat2 rot = mat2(cos(theta_0 - theta), -sin(theta_0 - theta), sin(theta_0 - theta), cos(theta_0 - theta));
	  oRGB_new = rot * oRGB;
	}

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
