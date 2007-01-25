
uniform float simtime;
uniform samplerRect dir_texunit;
uniform samplerRect pos_texunit;
uniform samplerRect wind_texunit;
uniform float nx;
uniform float ny;
uniform float nz;
uniform float numInRow;

void main(void)
{
   vec2 texCoord = gl_TexCoord[0].xy;
   vec3 pos = textureRect(pos_texunit, texCoord);
   
   int i = floor(pos.x);
   int j = floor(pos.z);
   int k = floor(pos.y);

   if((i != nx) && (j != nz) && (k != ny) && (i >= 0) && (j >= 0) && (k >=0)){

 	ivec2 index;
   	index.s = (int)(j + (mod(k,numInRow))*nz);
   	index.t = (int)(i + floor(k/numInRow)*nx);
	int s = index.s;
	int t = index.t;
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

	wind = normalize(wind);
   	pos = pos + wind*.001;  

   }
/*
   vec3 dir = textureRect(dir_texunit, texCoord);
   pos = pos + dir*.01;
   
   float d = sqrt(pos.x*pos.x + pos.y*pos.y + pos.z*pos.z);
   pos = (pos/d)*r;

	
   vec2 center = vec2(0.5, 0.5);
   float radius = sqrt((texCoord.x - center.x) * (texCoord.x - center.x) +
			(texCoord.y - center.y) * (texCoord.y - center.y) );

   vec4 dir = vec4(cos(radius), sin(radius), 0.0, 0.0);
   pos = pos + dir * simtime;*/
   // pos = pos + vec4(1.0, 0.0, 0.0, 0.0);
  
   gl_FragColor = vec4(pos, 1.0);
}
