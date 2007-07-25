uniform samplerRect pos_texunit;
uniform samplerRect primePrev;
uniform samplerRect wind_texunit;
uniform samplerRect lambda;
uniform samplerRect random;
uniform samplerRect tau_dz;
uniform samplerRect duvw_dz;

uniform int nx;
uniform int ny;
uniform int nz;
uniform float numInRow;
uniform float time_step;
uniform float life_time;

//Building variables
uniform float xfo;
uniform float yfo;
uniform float zfo;
uniform float ht;
uniform float wti;
uniform float lti;

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
   vec3 prmPrev = vec3(textureRect(primePrev, texCoord));
   vec3 prmCurr = prmPrev;
   
   //This gets the position of the particle in 3D space.
   vec4 pos = vec4(textureRect(pos_texunit, texCoord));
   vec3 prevPos = vec3(pos.x,pos.y,pos.z);
     
   // Add the random texture coordinate offset to the
   // texture coordinate.  The texture is set to perform wrapping
   // so texture coordinates outside the range will be valid.
   //
   // GL_TEXTURE_RECTANGLE_ARB does not support GL_REPEAT (at least on windows)
   // so we will do the normalization here given the width and length of the 
   // random texture.
   vec2 rTexCoord = texCoord + random_texCoordOffset;

   // bring the texture coordinate back within the (0,W)x(0,H) range
   if (rTexCoord.s > random_texWidth)
      rTexCoord.s = rTexCoord.s - random_texWidth;
   if (rTexCoord.t > random_texHeight)
      rTexCoord.t = rTexCoord.t - random_texHeight;

   // lookup the random value to be used for this particle in
   // this timestep
   vec3 randn = vec3(textureRect(random, rTexCoord));

   float xRandom = randn.x;
   float yRandom = randn.y;
   float zRandom = randn.z;     

   float upPrev=prmPrev.x;
   float vpPrev=prmPrev.y;
   float wpPrev=prmPrev.z;	


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
   	index.s = j + mod(k,numInRow)*nx;
   	index.t = i + floor(k/numInRow)*ny;
	
   	vec3 wind = vec3(textureRect(wind_texunit, index));    
	vec4 wind_tex = vec4(textureRect(wind_texunit, index));

	vec4 lam = vec4(textureRect(lambda, index));
	vec4 tau = vec4(textureRect(tau_dz, index));
	vec4 ddz = vec4(textureRect(duvw_dz, index));
	
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
			pow((2*CoEps_D2*time_step),0.5)* yRandom;

	float dw = (-CoEps_D2*(Lam13*upPrev + Lam33*wpPrev) + 0.5*Tau33)*time_step +
		(Tau13*(Lam11*upPrev + Lam13*wpPrev) + Tau33*(Lam13*upPrev + Lam33*wpPrev))*
		(wpPrev/2.0)*time_step + pow((2*CoEps_D2*time_step),0.5)* zRandom;

        //float du= -CoEps_D2*(Lam11*upPrev+Lam13*wpPrev)*time_step + pow((2*CoEps_D2*time_step),0.5)* xRandom;
        //float dv= -CoEps_D2*(Lam22*vpPrev)*time_step + pow((2*CoEps_D2*time_step),0.5)* yRandom;
        //float dw= -CoEps_D2*(Lam13*upPrev+Lam33*wpPrev)*time_step + pow((2*CoEps_D2*time_step),0.5)* zRandom;

	prmCurr = vec3(upPrev+du,vpPrev+dv,wpPrev+dw);

	//Now move the particle by adding the direction.
   	pos = pos + vec4(wind,0.0)*time_step + vec4(0.5*(prmPrev+prmCurr),0.0)*time_step;
	
	
	//Reflection off ground		
	vec3 u;
	vec3 w;
	//point of intersection
	vec3 pI;	


        while((pos.z < 0) || ((pos.x >= xfo) && (pos.x <= xfo+lti) && (pos.y >= yfo-(wti/2.0)) && 
		(pos.y <= yfo+(wti/2.0)) && (pos.z <= zfo+ht))){

		if(pos.z < 0){
			pos.z = -pos.z;
			prmCurr.z = -prmCurr.z;
			//pos = reflect(pos,n);
			//prmCurr = reflect(prmCurr,vec3(0.0,0.0,1.0));
		}

		//Reflection off building
		//Check to see if particle is inside building
		else if((pos.x >= xfo) && (pos.x <= xfo+lti) && (pos.y >= yfo-(wti/2.0)) && (pos.y <= yfo+(wti/2.0)) && (pos.z <= zfo+ht)){
		
			u = vec3(pos.x,pos.y,pos.z) - prevPos;	
			//plane facing -x direction
			w = prevPos - vec3(xfo,0.0,0.0);
			float s1 = dot(vec3(1.0,0.0,0.0),w)/dot(vec3(-1.0,0.0,0.0),u);
			//plane facing +x direction
			w = prevPos - vec3(xfo+lti,0.0,0.0);
			float s2 = dot(vec3(-1.0,0.0,0.0),w)/dot(vec3(1.0,0.0,0.0),u);
			//plane facing +y direction
			w = prevPos - vec3(xfo,(yfo+wti/2.0),0.0);
			float s3 = dot(vec3(0.0,-1.0,0.0),w)/dot(vec3(0.0,1.0,0.0),u);
			//plane facing -y direction
			w = prevPos - vec3(xfo,(yfo-wti/2.0),0.0);
			float s4 = dot(vec3(0.0,1.0,0.0),w)/dot(vec3(0.0,-1.0,0.0),u);
			//plane facing +z direction
			w = prevPos -vec3(xfo,0.0,(zfo+ht));
			float s5 = dot(vec3(0.0,0.0,-1.0),w)/dot(vec3(0.0,0.0,1.0),u);

			//incident vector
			vec3 l;
			//normal vector
			vec3 normal;
	
			if((s1 >= 0.0) && (s1 <= 1.0)){
				pI = s1*u + prevPos;
				normal = vec3(-1.0,0.0,0.0);				
			}
			else if((s2 >= 0.0) && (s2 <= 1.0)){
				pI = s2*u + prevPos;
				normal = vec3(1.0,0.0,0.0);
			}
			else if((s3 >= 0.0) && (s3 <= 1.0)){
				pI = s3*u + prevPos;
				normal = vec3(0.0,1.0,0.0);
			}
			else if((s4 >= 0.0) && (s4 <= 1.0)){
				pI = s4*u + prevPos;
				normal = vec3(0.0,-1.0,0.0);
			}
			else if((s5 >= 0.0) && (s5 <= 1.0)){
				pI = s5*u + prevPos;
				normal = vec3(0.0,0.0,1.0);
			}
			l = prevPos - pI;
			l = normalize(l);
			pos = vec4(pI,0.0) - vec4(reflect(l,normal),0.0);
			prmCurr = pI - reflect(l,normal);

		}
	}
	
   }
   if(pos.a <= 0 && (!(life_time <= 0))){
      gl_FragData[0] = vec4(100.0, 100.0, 100.0, life_time+1.0);
      gl_FragData[1] = vec4(prmCurr, 1.0);
   }
   else{
      gl_FragData[0] = pos;
      gl_FragData[1] = vec4(prmCurr, 1.0);
   }
   
}