uniform samplerRect pos_texunit;
uniform samplerRect wind_texunit;
uniform samplerRect random_texunit;
uniform float nx;
uniform float ny;
uniform float nz;
uniform float numInRow;
uniform float time_step;
uniform vec2 random_texCoordOffset;

void main(void)
{
   //This gets the position of the particle in 3D space.
   vec2 texCoord = gl_TexCoord[0].xy;
   vec4 pos = vec4(textureRect(pos_texunit, texCoord));
    
   // Read out a random value based on the particle's texture coordinate for use with 
   // the turbulence model.
   float turbulence_sigma = 0.5;

   // random values are being generated between -1.0 and 1.0
   vec3 turbulence = vec3(textureRect(random_texunit, texCoord));

   //The floor of the position in 3D space is needed to find the index into
   //the 2D Texture.
   float i = floor(pos.x);
   float j = floor(pos.z);
   float k = floor(pos.y);
	

   //This statement doesn't allow particles to move outside the domain.
   if((i < nx) && (j < nz) && (k < ny) && (i >= 0) && (j >= 0) && (k >=0)){

	//This is the initial lookup into the 2D texture that holds the wind field.
 	vec2 index;
   	index.s = j + mod(k,numInRow)*nz;
   	index.t = i + floor(k/numInRow)*nx;
	float s = index.s;
	float t = index.t;
	vec3 wind2;
   	vec3 wind = vec3(textureRect(wind_texunit, index));

	//Calculates distances to the edges of the surrounding cells.
	float dr = pos.z - mod(s,nz); 	//distance to the left edge.
	float dl = 1 - dr; 		//distance to the right edge.
	float da = pos.x - mod(t,nx); 	//distance to the edge below(in the same layer).
	float db = 1 - da; 		//distance to the edge above(in the same layer).
	float dk = pos.y - k; 		//distance to the layer below.
	float dd = 1 - dk; 		//distance to the layer above.

	//Perform lookups into the 2D Texture of the surrounding cells.
	//Then accumulate the weighted values.
 
	//To the Left
	if(mod(s,nz) != 0){
		index.s = s-1;
		index.t = t;
		wind2 = vec3(textureRect(wind_texunit, index));
		wind = wind + dl*wind2;
	}
	//To the Right
	if(mod(s,nz) != nz-1){
		index.s = s + 1;
		index.t = t;
		wind2 = vec3(textureRect(wind_texunit, index));
		wind = wind+ dr*wind2;
	}
	//Right Above
	if(mod(t,nx) != nx-1){
		index.s = s;
		index.t = t+1;
		wind2 = vec3(textureRect(wind_texunit, index));
		wind = wind+ wind2*da;
	}
	//Right Below
	if(mod(t,nx) != 0){
		index.s = s;
		index.t = t-1;
		wind2 = vec3(textureRect(wind_texunit, index));
		wind = wind+ wind2*db;
	}	
	//K level up
	if(k != ny-1){
		index.s = s;
		index.t = t + nx;
		wind2 = vec3(textureRect(wind_texunit, index));
		wind = wind+ wind2*dk;
	}
	//K level down
	if(k != 0){
		index.s = s;
		index.t = t-nx;
		wind2 = vec3(textureRect(wind_texunit, index));
		wind = wind+ wind2*dd;
	}		

	//Now move the particle by adding the direction.
   	pos = pos + vec4(wind,0.0)*time_step + vec4(turbulence,0.0)*time_step;

   }

   gl_FragColor = pos;
}
