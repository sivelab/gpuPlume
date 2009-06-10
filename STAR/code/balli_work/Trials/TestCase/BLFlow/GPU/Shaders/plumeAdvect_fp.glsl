uniform samplerRect pos_texunit;
uniform samplerRect primePrev;
uniform samplerRect primeCurr;
uniform samplerRect wind_texunit;
uniform samplerRect random_texunit;


uniform int nx;
uniform int ny;
uniform int nz;
uniform int nxdx;
uniform int nydy;
uniform int nzdz;
uniform int numInRow;
uniform float time_step;
uniform float life_time;

void main(void)
{
   vec2 texCoord = gl_TexCoord[0].xy;
   vec3 prmPrev = vec3(textureRect(primePrev, texCoord));
   vec3 prmCurr = vec3(textureRect(primeCurr, texCoord));
   
   //This gets the position of the particle in 3D space.
   vec4 pos = vec4(textureRect(pos_texunit, texCoord));
     
      
   //The floor of the position in 3D space is needed to find the index into
   //the 2D Texture.
   int i = floor(pos.y);
   int j = floor(pos.x);
   int k = floor(pos.z);
  
   if(!(life_time <= 0)){
     if(pos.a != life_time+1.0){
       //Decrement the alpha value based on the time_step
       pos.a = pos.a - time_step;

     }
   }   
	
   //This statement doesn't allow particles to move outside the domain.
   //So, the particles inside the domain are the only ones operated on. 
   if((i < ny) && (j < nx) && (k < nz) && (i >= 0) && (j >= 0) && (k >= 0)){
   

    //This is the initial lookup into the 2D texture that holds the wind field.
 	vec2 index;
   	index.s = j + mod(float(k),float(numInRow))*float(nxdx);
   	index.t = i + floor(float(k)/float(numInRow))*float(nydy);
	int s = index.s;
	int t = index.t;
	vec3 wind2;
   	vec3 wind = vec3(textureRect(wind_texunit, index));

	//Calculates distances to the edges of the surrounding cells.
	//float dr = pos.x - mod(s,nx); 	//distance to the left edge.
	//float dl = 1 - dr; 		//distance to the right edge.
	//float da = pos.y - mod(t,ny); 	//distance to the edge below(in the same layer).
	//float db = 1 - da; 		//distance to the edge above(in the same layer).
	//float dk = pos.z - k; 		//distance to the layer below.
	//float dd = 1 - dk; 		//distance to the layer above.


      // FOLLOWING IS COMMENTED BECAUSE IT DOUBLES THE WIND SPEED AND NOT REQUIRED FOR UNIFORM FLOW CASE--BALLI(05/01/07)
	//Perform lookups into the 2D Texture of the surrounding cells.
	//Then accumulate the weighted values.
 
	//To the Left
	//if(mod(s,nx) != 0){
	//	index.s = s-1;
   	//	index.t = t;
	//	wind2 = vec3(textureRect(wind_texunit, index));
	//	wind = wind + dl*wind2;
	//}
	//To the Right
	//if(mod(s,nx) != nx-1){
	//	index.s = s + 1;
	//	index.t = t;
	//	wind2 = vec3(textureRect(wind_texunit, index));
	//	wind = wind+ dr*wind2;
	//}
	//Right Above
	//if(mod(t,ny) != ny-1){
	//	index.s = s;
	//	index.t = t+1;
	//	wind2 = vec3(textureRect(wind_texunit, index));
	//	wind = wind+ wind2*da;
	//}
	//Right Below
	//if(mod(t,ny) != 0){
	//	index.s = s;
	//	index.t = t-1;
	//	wind2 = vec3(textureRect(wind_texunit, index));
	//	wind = wind+ wind2*db;
	//}	
	//K level up
	//if(k != nz-1){
	//	index.s = s;
	//	index.t = t + ny;
	//	wind2 = vec3(textureRect(wind_texunit, index));
	//	wind = wind+ wind2*dk;
	//}
	//K level down
	//if(k != 0){
	//	index.s = s;
	//	index.t = t-ny;
	//	wind2 = vec3(textureRect(wind_texunit, index));
	//	wind = wind+ wind2*dd;
	//}		
    //COMMENTING COMPLETE--BALLI(05/01/07)
	//Now move the particle by adding the direction.
   	pos = pos + vec4(wind,0.0)*time_step + vec4(0.5*(prmPrev+prmCurr),0.0)*time_step;
	
	
   }
   if(pos.a <= 0 && (!(life_time <= 0))){
      gl_FragColor = vec4(100.0, 100.0, 100.0, life_time+1.0);
   }
   else{
      gl_FragColor = pos;
   }
   
}
