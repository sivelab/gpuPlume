uniform samplerRect currPrime;
uniform samplerRect position;
uniform samplerRect position_prev;
uniform samplerRect windVel;

uniform float numInRow;

uniform int nx;
uniform int ny;
uniform int nz;

void main(void)
{
   vec2 texCoord = gl_TexCoord[0].xy;
   vec3 primeVel = vec3(textureRect(currPrime, texCoord));
   vec3 prevPos = vec3(textureRect(position_prev, texCoord));
   vec3 pos = vec3(textureRect(position, texCoord));

   //int i = floor(pos.y);
   //int j = floor(pos.x);
   //int k = floor(pos.z);

   //vec4 color = vec4(1.0,1.0,1.0,1.0);

   //if((i < ny) && (j < nx) && (k < nz) && (i >= 0) && (j >= 0) && (k >= 0)){
   //vec2 index;
   //index.s = j + mod(k,numInRow)*nx;
   //index.t = i + floor(k/numInRow)*ny;
	
   //vec3 wind = vec3(textureRect(windVel,index));
   //color = vec4(wind+primeVel,1.0);
   //color = abs(normalize(color));
   // vec3 color = abs(normalize(pos - prevPos));

   // This example assumes the difference in position sits in the
   // (Helmholtz style) opponent color space.  We then convert the
   // opponent color space to RGB to display on the screen.  This
   // should produce colors that oppose each other relative to the
   // direction. In other words, Red/Magenta opposes Green, Blue/Cyan
   // opposes Yellow, and Black opposes White.  We'll see how it
   // works.

   mat3 opponent2rgb = mat3(1.1677, -6.4315, -0.5044, 0.9014, 2.5970, 0.0159, 0.7214, 0.1257, 2.0517);
   vec3 opponent_color = normalize(pos - prevPos);
   vec3 color = opponent2rgb * opponent_color;
   	
   //}

   gl_FragColor = vec4(color,1.0);

}