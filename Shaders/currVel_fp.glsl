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

   int i = floor(pos.y);
   int j = floor(pos.x);
   int k = floor(pos.z);

   vec3 color = vec3(0.8, 0.8, 1.0);

   if((i < ny) && (j < nx) && (k < nz) && (i >= 0) && (j >= 0) && (k >= 0))
   {
      vec2 index;
      index.s = j + mod(k,float(numInRow))*nx;
      index.t = i + floor(k/float(numInRow))*ny;
	
      vec3 wind = vec3(textureRect(windVel,index));

      // This example assumes the difference in position sits in the
      // (Hering-style) opponent color space.  We then convert the
      // opponent color space to RGB to display on the screen.  This
      // should produce colors that oppose each other relative to the
      // direction. In other words, Red/Magenta opposes Green, Blue/Cyan
      // opposes Yellow, and Black opposes White.  We'll see how it
      // works.

      mat3 opponent2rgb = mat3(1.0,  0.1140,  0.7436, 
                               1.0,  0.1140, -0.4111, 
                               1.0, -0.8860, 0.1663);
      vec3 opponent_color = normalize(pos - prevPos);

      // or, get the wind field
      // opponent_color = normalize(wind);

      opponent_color.x = abs(opponent_color.x);
      color = opponent2rgb * opponent_color;

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
   }
   
   gl_FragColor = vec4(color,1.0);
}
