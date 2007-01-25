
uniform samplerRect pos_texunit;
uniform samplerRect wind_texunit;
uniform float nx;
uniform float ny;
uniform float nz;
uniform float numInRow;
uniform float time_step;

void main(void)
{
   vec2 texCoord = gl_TexCoord[0].xy;
   vec4 pos = vec4(textureRect(pos_texunit, texCoord));
    
   float i = floor(pos.x);
   float j = floor(pos.z);
   float k = floor(pos.y);
	
   if((i < nx) && (j < nz) && (k < ny) && (i >= 0) && (j >= 0) && (k >=0)){

 	vec2 index;
   	index.s = j + mod(k,numInRow)*nz;
   	index.t = i + floor(k/numInRow)*nx;
	float s = index.s;
	float t = index.t;
	vec3 wind2;
   	vec3 wind = vec3(textureRect(wind_texunit, index));

	//Calculates distances to the edges of a cell
	float dr = pos.z - mod(s,nz);
	float dl = 1 - dr;
	float da = pos.x - mod(t,nx);
	float db = 1 - da;
	float dk = pos.y - k;
	float dd = 1- dk;

	//To the Left
	if(mod(s,nz) != 0){
		index.s = s-1;
		index.t = t;
		wind2 = vec3(textureRect(wind_texunit, index));
		wind = wind + dl*wind2;
		//wind = normalize(wind + wind2);
	}
	//To the Right
	if(mod(s,nz) != nz-1){
		index.s = s + 1;
		index.t = t;
		wind2 = vec3(textureRect(wind_texunit, index));
		wind = wind+ dr*wind2;
		//wind = normalize(wind + wind2);
	}
	//Right Above
	if(mod(t,nx) != nx-1){
		index.s = s;
		index.t = t+1;
		wind2 = vec3(textureRect(wind_texunit, index));
		wind = wind+ wind2*da;
		//wind = normalize(wind + wind2);
	}
	//Right Below
	if(mod(t,nx) != 0){
		index.s = s;
		index.t = t-1;
		wind2 = vec3(textureRect(wind_texunit, index));
		wind = wind+ wind2*db;
		//wind = normalize(wind + wind2);
	}	
	//K level up
	if(k != ny-1){
		index.s = s;
		index.t = t + nx;
		wind2 = vec3(textureRect(wind_texunit, index));
		wind = wind+ wind2*dk;
		//wind = normalize(wind + wind2);
	}
	//K level down
	if(k != 0){
		index.s = s;
		index.t = t-nx;
		wind2 = vec3(textureRect(wind_texunit, index));
		wind = wind+ wind2*dd;
		//wind = normalize(wind + wind2);
	}		

	//wind = normalize(wind);
   	pos = pos + vec4(wind,0.0)*time_step;  

   }

   gl_FragColor = pos;
}
