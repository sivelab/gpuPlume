uniform sampler2DRect pos_texunit;
uniform sampler2DRect primePrev;
uniform sampler2DRect wind_texunit;
uniform sampler2DRect lambda;
uniform sampler2DRect random;
uniform sampler2DRect tau_dz;
uniform sampler2DRect duvw_dz;

uniform int nx;
uniform int ny;
uniform int nz;
uniform int nxdx;
uniform int nydy;
uniform int nzdz;
uniform int numInRow;
uniform float time_step;
uniform float life_time;

//
// This variable contains a uniformally generated random 2D vector in
// the range [0,1].  The value is generated on the CPU each iteration
// and is used to offset the random texture (used in this shader).  By
// offsetting the random texture lookup, the random numbers will
// reflect a more uniform random value rather than the same value for
// each particle.
//
uniform int random_texWidth;
uniform int random_texHeight;
uniform vec2 random_texCoordOffset;

void main(void)
{
   vec2 texCoord = gl_TexCoord[0].xy;
   vec3 prmPrev = vec3(texture2DRect(primePrev, texCoord));
   vec3 prmCurr = prmPrev;
   
   //This gets the position of the particle in 3D space.
   vec4 pos = vec4(texture2DRect(pos_texunit, texCoord));
     
   // Add the random texture coordinate offset to the
   // texture coordinate.  The texture is set to perform wrapping
   // so texture coordinates outside the range will be valid.
   //
   // GL_TEXTURE_RECTANGLE_ARB does not support GL_REPEAT (at least on windows)
   // so we will do the normalization here given the width and length of the 
   // random texture.
   vec2 rTexCoord = texCoord + random_texCoordOffset;

   // bring the texture coordinate back within the (0,W)x(0,H) range
   if (rTexCoord.s > float(random_texWidth))
      rTexCoord.s = rTexCoord.s - float(random_texWidth);
   if (rTexCoord.t > float(random_texHeight))
      rTexCoord.t = rTexCoord.t - float(random_texHeight);

   // lookup the random value to be used for this particle in
   // this timestep
   vec3 randn = vec3(texture2DRect(random, rTexCoord));

   float xRandom = randn.x;
   float yRandom = randn.y;
   float zRandom = randn.z;     

   float upPrev=prmPrev.x;
   float vpPrev=prmPrev.y;
   float wpPrev=prmPrev.z;	


   //The floor of the position in 3D space is needed to find the index into
   //the 2D Texture.
   int i = int(floor(pos.y));
   int j = int(floor(pos.x));
   int k = int(floor(pos.z));
  
   if(!(life_time <= 0.0)){
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
   	index.s = float(j) + float(mod(float(k),float(numInRow)))*float(nxdx);
   	index.t = float(i) + float(floor(float(k)/float(numInRow)))*float(nydy);
	
   	vec3 wind = vec3(texture2DRect(wind_texunit, index));
	vec4 wind_tex = vec4(texture2DRect(wind_texunit, index));

	vec4 lam = vec4(texture2DRect(lambda, index));
	vec4 tau = vec4(texture2DRect(tau_dz, index));
	vec4 ddz = vec4(texture2DRect(duvw_dz, index));
	
	float dudz = ddz.x;

	float Tau11 = tau.x;
	float Tau22 = tau.y;
	float Tau33 = tau.z;
	float Tau13 = tau.w;

	float Lam11=lam.x;
	float Lam22=lam.y;
	float Lam33=lam.z;
	float Lam13=lam.w;
	
	float CoEps_D2=wind_tex.w; // grabing 4th vector,--Co*Eps/2 in the wind texture
      
	float du = (-CoEps_D2*(Lam11*upPrev+Lam13*wpPrev) + dudz*wpPrev + 0.5*Tau13)*time_step
		+ (Tau11*(Lam11*upPrev + Lam13*wpPrev) + Tau13*(Lam13*upPrev + Lam33*wpPrev))*
		(wpPrev/2.0)*time_step + pow((2.0*CoEps_D2*time_step),0.5)*xRandom;
	
	float dv = (-CoEps_D2*(Lam22*vpPrev) + Tau22*Lam22*vpPrev*(wpPrev/2.0))*time_step + 
			pow((2.0*CoEps_D2*time_step),0.5)* yRandom;

	float dw = (-CoEps_D2*(Lam13*upPrev + Lam33*wpPrev) + 0.5*Tau33)*time_step +
		(Tau13*(Lam11*upPrev + Lam13*wpPrev) + Tau33*(Lam13*upPrev + Lam33*wpPrev))*
		(wpPrev/2.0)*time_step + pow((2.0*CoEps_D2*time_step),0.5)* zRandom;

        //float du= -CoEps_D2*(Lam11*upPrev+Lam13*wpPrev)*time_step + pow((2*CoEps_D2*time_step),0.5)* xRandom;
        //float dv= -CoEps_D2*(Lam22*vpPrev)*time_step + pow((2*CoEps_D2*time_step),0.5)* yRandom;
        //float dw= -CoEps_D2*(Lam13*upPrev+Lam33*wpPrev)*time_step + pow((2*CoEps_D2*time_step),0.5)* zRandom;

	prmCurr = vec3(upPrev+du,vpPrev+dv,wpPrev+dw);

	//Now move the particle by adding the direction.
   	pos = pos + vec4(wind,0.0)*time_step + vec4(0.5*(prmPrev+prmCurr),0.0)*time_step;
	
	
	//Reflection off ground	
	vec4 n = vec4(0.0,0.0,1.0,0.0);
	float rdot = dot(pos,n);		

	if(pos.z < 0.0){
		pos.z = -pos.z;
		prmCurr.z = -prmCurr.z;
		//pos = reflect(pos,n);
		//prmCurr = reflect(prmCurr,vec3(0.0,0.0,1.0));
	}

   }
   if(pos.a <= 0.0 && (!(life_time <= 0.0))){
      gl_FragData[0] = vec4(100.0, 100.0, 100.0, life_time+1.0);
      gl_FragData[1] = vec4(prmCurr, 1.0);
   }
   else{
      gl_FragData[0] = pos;
      gl_FragData[1] = vec4(prmCurr, 1.0);
   }
   
}
