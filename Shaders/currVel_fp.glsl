uniform samplerRect currPrime;
uniform samplerRect position;
uniform samplerRect position_prev;
uniform samplerRect windVel;

uniform int numInRow;

uniform int nx;
uniform int ny;
uniform int nz;

void main(void)
{
   vec2 texCoord = gl_TexCoord[0].xy;
   vec3 primeVel = vec3(textureRect(currPrime, texCoord));
   vec3 prevPos = vec3(textureRect(position_prev, texCoord));
   vec3 pos = vec3(textureRect(position, texCoord));

   //   int i = floor(pos.y);
   //   int j = floor(pos.x);
   // Is this right???
   int j = floor(pos.y);
   int i = floor(pos.x);
   int k = floor(pos.z);

   //vec3 color = vec3(0.8, 0.8, 1.0);
	vec3 color_prime;
	//    vec3 color = vec3(0.0, 0.0, 0.0);
	//   vec3 color = vec3((float)j/(float)ny, (float)i/(float)ny, (float)k/(float)nz);
   vec3 color = vec3(1.0,1.0,1.0);
   if((i < ny) && (j < nx) && (k < nz) && (i >= 0) && (j >= 0) && (k >= 0))
   {
      vec2 index;
      index.s = float(j) + float(mod(float(k),float(numInRow))*nx);
      index.t = float(i) + float(floor(float(k)/float(numInRow))*ny);

      vec3 wind = vec3(textureRect(windVel,index));

      // This example assumes the difference in position sits in the
      // (Hering-style) opponent color space.  We then convert the
      // opponent color space to RGB to display on the screen.  This
      // should produce colors that oppose each other relative to the
      // direction. In other words, Red/Magenta opposes Green, Blue/Cyan
      // opposes Yellow, and Black opposes White.  We'll see how it
      // works.

      // We are taking the normalized difference (velocity direction)
      // of the particle positions as the L'C'C' intermediate
      // parameters of the oRGB conversion.

      // We map height (Z) to the luma (L') channel, direction in X to
      // the C'_{1} channel and direction in Y to the C'_{2} channel.
      vec3 opponent_color = normalize(pos - prevPos).zxy;
      // vec3 opponent_color = wind.zxy;

      color_prime = opponent_color;

      opponent_color.y = opponent_color.y * sqrt(2.0);
      opponent_color.z = opponent_color.z * sqrt(2.0);

      // The above computation will place the opponent_color values in
      // the range [-1, 1] since directions are within the unit
      // sphere.  The L'C'C' 
      // 
      // Map luma to the [0, 1] range
      opponent_color.x = (1.0 + opponent_color.x)/2.0;

      // The next step is to perform a linear transform 
      //      mat3 opponent2rgb = mat3(1.0,  0.1140,  0.7436, 
      //      			       1.0,  0.1140, -0.4111, 
      //                               1.0, -0.8860,  0.1663);
      //      color = opponent2rgb * opponent_color;

      // why is this done when we're already taking the points as being in the linear LCC space...
      mat3 opponent2rgb = mat3(0.299,  0.587,  0.114, 
      			       0.500,  0.500, -1.000, 
      			       0.866, -0.886,  0.000);
      color = opponent2rgb * opponent_color;

      // Now, we need to apply a non-affine transform to make the
      // red/green axis perpendicular to the blue/yellow axis

      float pi = 3.141592654;
      float pi_3 = 1.047197551;
      float theta_0 = 0.0;

      // In GLSL, the atan function with 2 arguments computes the
      // angle whose tangent is y/x.  The signs of x and y determine
      // what quadrant the angle is in and the range of values
      // returned are [-pi, pi].  Results are undefined if x and y are
      // both 0.  This behavior is similar to atan2 in the Unix math
      // library, but not the same...

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

	  // value of theta will/should be in [0 to -pi] range
	  // if (theta >= pi_3)
	  // theta_0 = 1.5 * theta;
	  // else if (theta >= pi && theta < pi_3)
	  // theta_0 = pi/4.0 + 0.75*theta;
	  //	theta_0 = pi/2.0 + 0.75*(theta - pi/3.0);
	  
	  // oRGB = vec2(color.y, color.z);
	  // mat2 rot = mat2(cos(theta - theta_0), -sin(theta - theta_0), sin(theta - theta_0), cos(theta - theta_0));
	  // oRGB_new = rot * oRGB;
	  
	}
      
      color.y = oRGB_new.x;  
      color.z = oRGB_new.y;
   }
   
   gl_FragColor = vec4(color,1.0);
}
