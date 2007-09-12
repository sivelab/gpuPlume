uniform samplerRect pos_texunit;
uniform samplerRect primePrev;
uniform samplerRect wind_texunit;
uniform samplerRect lambda;
uniform samplerRect random;
uniform samplerRect tau_dz;
uniform samplerRect duvw_dz;
uniform samplerRect buildings;
uniform samplerRect cellType;

uniform int nx;
uniform int ny;
uniform int nz;
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
    index.s = float(j) + float(mod(float(k),float(numInRow)))*float(nx);
    index.t = float(i) + float(floor(float(k)/float(numInRow)))*float(ny);

    vec4 wind = vec4(textureRect(wind_texunit, index));

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
	
    float CoEps_D2=wind.w; // grabing 4th vector,--Co*Eps/2 in the wind texture
      
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
    pos = pos + vec4(wind.x,wind.y,wind.z,0.0)*time_step + vec4(0.5*(prmPrev+prmCurr),0.0)*time_step;
	
    //Now do Reflection		
    vec3 u;
    //point of intersection
    vec3 pI;	
    //incident vector
    vec3 l;
    //reflection vector
    vec3 r;
    //normal vector
    vec3 normal;
    //distance from reflected surface
    float dis;
                  
    float denom;
    float numer;
	
    ivec2 cIndex;	

    i = int(floor(pos.y));
    j = int(floor(pos.x));
    k = int(floor(pos.z));
    int count = 0;	
    float eps = 0.0;
    float eps_S = 0.0001;
    float eps_d = 0.01;
    float smallestS = 500.0;

    if((i < ny) && (j < nx) && (k < nz) && (i >= 0) && (j >= 0)){
      vec4 cell_type = vec4(1.0,0.0,0.0,1.0);
	
      if(k < 0)
	k = 0;

      if(k >= 0){
	//cIndex.s = float(j) + float(mod(k,numInRow))*float(nx);
	//cIndex.t = float(i) + float(floor(float(k)/float(numInRow)))*float(ny);
	cIndex.s = j + mod(k,numInRow)*nx;
	cIndex.t = i + int(floor(float(k)/float(numInRow)))*ny;
	//Perform lookup into wind texture to see if new position is inside a building
	cell_type = vec4(textureRect(cellType, cIndex));  
      }

      while(((cell_type.x == 0.0 && cell_type.y == 0.0 && cell_type.z == 0.0) || (pos.z < eps)) && (count < 20)){
	count = count + 1;
	u = vec3(pos) - prevPos;

	float d;
	vec3 n;
	float s1 = -1.0;
	float s2 = -1.0;
	float s3 = -1.0;
	float s4 = -1.0;
	float s5 = -1.0;
	float s6 = -1.0;
	float s7 = -1.0;
	
	smallestS = 100.0;
	
	int id = int(cell_type.w);
	ivec2 bindex;
	bindex.s = 0;
	bindex.t = id;
		
	vec3 bcoords = vec3(textureRect(buildings, bindex));
	float xfo = bcoords.x;
	float yfo = bcoords.y;
	float zfo = bcoords.z;

	bindex.x = 1;
	vec3 bdim = vec3(textureRect(buildings,bindex));
	float ht = bdim.x;
	float wti = bdim.y;
	float lti = bdim.z;
	  
	//-x normal  
	n = vec3(-1.0,0.0,0.0);
	d = -dot(n,vec3(xfo,0.0,0.0));
	denom = dot(n,u);
	numer = dot(n,prevPos) + d;
	s1 = -numer/denom;
      
	//+x normal
	n = vec3(1.0,0.0,0.0);
	d = -dot(n,vec3(xfo+lti,0.0,0.0));
	denom = dot(n,u);
	numer = dot(n,prevPos) + d;
	s2 = -numer/denom;
        
	//+y normal
	n = vec3(0.0,1.0,0.0);
	d = -dot(n,vec3(xfo,yfo+(wti/2.0),0.0));
	denom = dot(n,u);
	numer = dot(n,prevPos) + d;
	s3 = -numer/denom;
           
	//-y normal
	n = vec3(0.0,-1.0,0.0);
	d = -dot(n,vec3(xfo,yfo-(wti/2.0),0.0));
	denom = dot(n,u);
	numer = dot(n,prevPos) + d;
	s4 = -numer/denom;
      	
	//+z normal
	n = vec3(0.0,0.0,1.0);
	d = -dot(n,vec3(xfo,0.0,zfo+ht));
	denom = dot(n,u);
	numer = dot(n,prevPos) + d;
	s5 = -numer/denom;
	  
	//-z normal
	/*n = vec3(0.0,0.0,-1.0);
	  d = -dot(n,vec3(xfo,0.0,zfo));
	  denom = dot(n,u);
	  numer = dot(n,prevPos) + d;
	  s6 = -numer/denom;*/
	//Ground plane
	n = vec3(0.0,0.0,1.0);
	numer = dot(n,prevPos);
	denom = dot(n,u);
	s7 = -numer/denom;

	if((s1 < smallestS) && (s1 >= -eps_S)){
	  smallestS = s1;
	  normal = vec3(-1.0,0.0,0.0);
	}	
	if((s2 < smallestS) && (s2 >= -eps_S)){
	  normal = vec3(1.0,0.0,0.0);
	  smallestS = s2;
	}
	if((s3 < smallestS) && (s3 >= -eps_S)){
	  normal = vec3(0.0,1.0,0.0);
	  smallestS = s3;
	}	
	if((s4 < smallestS) && (s4 >= -eps_S)){
	  normal = vec3(0.0,-1.0,0.0);
	  smallestS = s4;
	}	   
	if((s5 < smallestS) && (s5 >= -eps_S)){
	  normal = vec3(0.0,0.0,1.0);
	  smallestS = s5;
	}
	/*if((s6 < smallestS) && (s6 >= -eps_S)){
	  normal = vec3(0.0,0.0,-1.0);
	  smallestS = s6;
	  }*/

	//Detect Edge Collision
	float edgeS = abs(smallestS-s7);
	if((edgeS < eps_d)){
	  //smallestS = (s7 + smallestS)/2.0;
	  normal = normalize(normal+vec3(0.0,0.0,1.0));
	}
	else if((s7 < smallestS) && (s7 >= -eps_S)){
	  normal = vec3(0.0,0.0,1.0);
	  smallestS = s7;                  
	}	
       	  

	pI = smallestS*u + prevPos;
	if((smallestS >= -eps_S) && (smallestS <= eps_S)){
	  pI = prevPos;
	  r = normal;
	}	
	else{
	  l = normalize(pI-prevPos);
	  r = normalize(reflect(l,normal));
	  
	}
	dis = distance(pI,vec3(pos)); 
			
      	
	prevPos = pI;
	
	pos = vec4(pI+(dis*r),pos.a);
	prmCurr = reflect(prmCurr,normal);
      
	//After Reflection, see where the particle is in preparation for continuing while loop or not
	i = int(floor(pos.y));
	j = int(floor(pos.x));
	k = int(floor(pos.z));

	//NOTE: Consider what happens if building is too close to domain.
	//Do check to make sure i,j,k's are valid;
	cell_type = vec4(1.0,1.0,1.0,1.0);
	if(k < 0)
	  k = 0;
	if(k >= 0){
	  //cIndex.s = float(j) + float(mod(k,numInRow))*float(nx);
	  //cIndex.t = float(i) + float(floor(float(k)/float(numInRow)))*float(ny);
	  cIndex.s = j + mod(k,numInRow)*nx;
	  cIndex.t = i + int(floor(float(k)/float(numInRow)))*ny;
	  cell_type = vec4(textureRect(cellType, cIndex));
	}
      }
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
