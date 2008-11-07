uniform sampler2DRect prevMean;
uniform sampler2DRect currVel;
uniform sampler2DRect position;
uniform sampler2DRect windVel;

uniform int numInRow;
uniform int nx;
uniform int ny;
uniform int nz;

void main(void)
{
   vec2 texCoord = gl_TexCoord[0].xy;
   vec3 velocity = vec3(texture2DRect(currVel, texCoord));
   vec4 mean = vec4(texture2DRect(prevMean, texCoord));
   vec3 pos = vec3(texture2DRect(position, texCoord));

   int i = floor(pos.y); 
   int j = floor(pos.x);
   int k = floor(pos.z);

   if((i < ny) && (j < nx) && (k < nz) && (i >= 0) && (j >= 0) && (k >= 0)){
	vec2 index;
	index.s = float(j) +float(mod(float(k),float(numInRow))*nx);
	index.t = float(i) + float(floor(float(k)/float(numInRow))*ny);
	
	vec3 wind = vec3(texture2DRect(windVel,index));

   	mean = mean + vec4(velocity,1.0) + vec4(wind,0.0);
   }

   gl_FragColor = mean;

}
