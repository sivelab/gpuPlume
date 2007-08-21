uniform samplerRect currPrime;
uniform samplerRect position;
uniform samplerRect windVel;

uniform float numInRow;

uniform int nx;
uniform int ny;
uniform int nz;

void main(void)
{
   vec2 texCoord = gl_TexCoord[0].xy;
   vec3 primeVel = vec3(textureRect(currPrime, texCoord));
   //vec4 windV = vec4(textureRect(windVel, texCoord));
   vec3 pos = vec3(textureRect(position, texCoord));

   int i = floor(pos.y);
   int j = floor(pos.x);
   int k = floor(pos.z);

   vec4 color = vec4(1.0,1.0,1.0,1.0);

   //if((i < ny) && (j < nx) && (k < nz) && (i >= 0) && (j >= 0) && (k >= 0)){
	vec2 index;
	index.s = j + mod(k,numInRow)*nx;
	index.t = i + floor(k/numInRow)*ny;
	
	vec3 wind = vec3(textureRect(windVel,index));
	color = vec4(wind+primeVel,1.0);
	color = abs(normalize(color));// + vec4(1.0,1.0,1.0,0.0);
	//color = color/2.0;
   	
   //}

   gl_FragColor = color;

}