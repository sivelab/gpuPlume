uniform samplerRect prevMean;
uniform samplerRect currVel;
uniform samplerRect position;
uniform samplerRect windVel;

uniform int numInRow;
uniform int nx;
uniform int ny;
uniform int nz;

void main(void)
{
   vec2 texCoord = gl_TexCoord[0].xy;
   vec3 velocity = vec3(textureRect(currVel, texCoord));
   vec4 mean = vec4(textureRect(prevMean, texCoord));
   vec3 pos = vec3(textureRect(position, texCoord));

   int i = floor(pos.y); 
   int j = floor(pos.x);
   int k = floor(pos.z);

   if((i < ny) && (j < nx) && (k < nz) && (i >= 0) && (j >= 0) && (k >= 0)){
	vec2 index;
	index.s = float(j) +float(mod(float(k),float(numInRow))*nx);
	index.t = float(i) + float(floor(float(k)/float(numInRow))*ny);
	
	vec3 wind = vec3(textureRect(windVel,index));

   	mean = mean + vec4(velocity,1.0) + vec4(wind,0.0);
   }

   gl_FragColor = mean;

}
