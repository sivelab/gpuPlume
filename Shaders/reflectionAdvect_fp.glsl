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
uniform int nxdx;
uniform int nydy;
uniform int nzdz;
uniform int numInRow;
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
float check=1.0;  // to check how many times CFL condition got applied

void main(void){

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
  //===============For CFL Condition===========================
  float dx=1.0;  //defining grid size
  float dy=1.0;
  float dz=1.0;
  float epsilon=0.0001; //to avoid precision issues when checking how much time is remaining
   
  float timeStepRem = time_step; // time step remaining at this point
  float timeStepSim = time_step;  //time step to be used in the simulation
  float timeStepUsed = 0.0;   //used for storing how much time has been used from a time step during simulation
  bool loopThrough=true; 
   
  check=1.0;   

  if(!(life_time <= 0.0)){
    if(pos.a != life_time+1.0){
      //Decrement the alpha value based on the time_step
      pos.a = pos.a - time_step;
          
    }
  }   
  
  while(loopThrough){
    loopThrough=false;  
        
          
    //The floor of the position in 3D space is needed to find the index into
    //the 2D Texture.
    int i = int(floor(pos.y));
    int j = int(floor(pos.x));
    int k = int(floor(pos.z));
     	
    //This statement doesn't allow particles to move outside the domain.
    //So, the particles inside the domain are the only ones operated on. 
    if((i < ny) && (j < nx) && (k < nz) && (i >= 0) && (j >= 0) && (k >= 0)){
        
      //This is the initial lookup into the 2D texture that holds the wind field.
      vec2 index;
      index.s = float(j) + float(mod(float(k),numInRow))*float(nxdx);
      index.t = float(i) + float(floor(float(k)/numInRow))*float(nydy);
	      
      vec4 wind = vec4(textureRect(wind_texunit, index));          
      vec4 lam = vec4(textureRect(lambda, index));
      vec4 tau = vec4(textureRect(tau_dz, index));
      vec4 ddz = vec4(textureRect(duvw_dz, index));
	     
      float dudz = ddz.x;
      float ustar = ddz.w; //temporarily used to get ustr to calculate sigmas in shader
  	      
      float sigU = 2.5*ustar;
      float sigV = 2.0*ustar;
      float sigW = 1.3*ustar; 
	     
      float Tau11 = tau.x;
      float Tau22 = tau.y;
      float Tau33 = tau.z;
      float Tau13 = tau.w;
	     
      float Lam11=lam.x;
      float Lam22=lam.y;
      float Lam33=lam.z;
      float Lam13=lam.w;
	     
      float CoEps_D2=wind.w; // grabing 4th vector,--Co*Eps/2 in the wind texture
           
      float du = (-CoEps_D2*(Lam11*upPrev+Lam13*wpPrev) + dudz*wpPrev + 0.5*Tau13)*timeStepSim
	+ (Tau11*(Lam11*upPrev + Lam13*wpPrev) + Tau13*(Lam13*upPrev + Lam33*wpPrev))*
	(wpPrev/2.0)*timeStepSim + pow((2.0*CoEps_D2*timeStepSim),0.5)*xRandom;
	         
      float dv = (-CoEps_D2*(Lam22*vpPrev) + Tau22*Lam22*vpPrev*(wpPrev/2.0))*timeStepSim + 
	pow((2.0*CoEps_D2*timeStepSim),0.5)* yRandom;
	     
      float dw = (-CoEps_D2*(Lam13*upPrev + Lam33*wpPrev) + 0.5*Tau33)*timeStepSim +
	(Tau13*(Lam11*upPrev + Lam13*wpPrev) + Tau33*(Lam13*upPrev + Lam33*wpPrev))*
	(wpPrev/2.0)*timeStepSim + pow((2.0*CoEps_D2*timeStepSim),0.5)* zRandom;
         
      float disX = (wind.x+du)*timeStepSim; //calculating distance travelled by the particle
      float disY = (wind.y+dv)*timeStepSim;
      float disZ = (wind.z+dw)*timeStepSim;
         
      if(disX<0.0)disX=-disX;//calculating the absolute value
      if(disY<0.0)disY=-disY;
      if(disZ<0.0)disZ=-disZ;
         
      if( disX>0.7*dx || disY>0.7*dy || disZ>0.7*dz ){  //CFL condition
	timeStepSim = timeStepSim/2.0;
          
	//generating random number
	rTexCoord.s = rTexCoord.s + 1.0;
	if (rTexCoord.s > random_texWidth)rTexCoord.s = rTexCoord.s - random_texWidth;
	vec3 randnum = vec3(textureRect(random, rTexCoord));
            
	upPrev = sigU * randnum.x;
	vpPrev = sigV * randnum.y;
	wpPrev = sigW * randnum.z;
           
	loopThrough=true; //make loopThrough true so that it calculate distance travelled again by new timeStepSim
	check=check+1.0;
      }
         
         
      if(loopThrough==false){ //move particles only if CFL condition allows 
         
	prmCurr = vec3(upPrev+du,vpPrev+dv,wpPrev+dw);
            
	//Now move the particle by adding the direction.
	pos = pos + vec4(wind.x,wind.y,wind.z,0.0)*timeStepSim + vec4(0.5*(prmPrev+prmCurr),0.0)*timeStepSim;
	    

	i = int(floor(pos.y));
	j = int(floor(pos.x));
	k = int(floor(pos.z));
    	
	vec3 u;
	vec3 w;
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
	float smallestS = 0.0;
	float eps_S = 0.00001;
	if((i < ny) && (j < nx) && (k < nz) && (i >= 0) && (j >= 0)){
	  while((pos.z < 0.0) || ((pos.x > xfo) && (pos.x < xfo+lti) && (pos.y > yfo-(wti/2.0)) && 
				  (pos.y < yfo+(wti/2.0)) && (pos.z < zfo+ht))  ){
	     
	    u = vec3(pos) - prevPos;	
	    //Reflection off ground
	    if(pos.z < 0.0){
	      normal = vec3(0.0,0.0,1.0);
	      numer = dot(normal,prevPos);
	      denom = dot(normal,u);
	      float s = -numer/denom;

	      pI = s*u + prevPos;
	      if(s > -eps_S && s < eps_S){
		r = normal;
	      }
	      else{
		l = normalize(pI-prevPos);
		r = reflect(l,normal);
	      }
	      dis = distance(pI,vec3(pos));

	      //Set prevPos to point of intersection
	      prevPos = pI;
                  
	      pos = vec4(pI+(dis*r),pos.a);
	      prmCurr = reflect(prmCurr,normal);
	    }	      
               
	    //Reflection off building
	    //Check to see if particle is inside building
	    else{
                       		     
	      float d;
	      vec3 n;
	      //-x normal  
	      n = vec3(-1.0,0.0,0.0);
	      d = -dot(n,vec3(xfo,0.0,0.0));
	      denom = dot(n,u);
	      numer = dot(n,prevPos) + d;
	      float s1 = -numer/denom;
      
	      //+x normal
	      n = vec3(1.0,0.0,0.0);
	      d = -dot(n,vec3(xfo+lti,0.0,0.0));
	      denom = dot(n,u);
	      numer = dot(n,prevPos) + d;
	      float s2 = -numer/denom;
        
	      //+y normal
	      n = vec3(0.0,1.0,0.0);
	      d = -dot(n,vec3(xfo,yfo+(wti/2.0),0.0));
	      denom = dot(n,u);
	      numer = dot(n,prevPos) + d;
	      float s3 = -numer/denom;
           
	      //-y normal
	      n = vec3(0.0,-1.0,0.0);
	      d = -dot(n,vec3(xfo,yfo-(wti/2.0),0.0));
	      denom = dot(n,u);
	      numer = dot(n,prevPos) + d;
	      float s4 = numer/denom;
      	
	      //+z normal
	      n = vec3(0.0,0.0,1.0);
	      d = -dot(n,vec3(xfo,0.0,zfo+ht));
	      denom = dot(n,u);
	      numer = dot(n,prevPos) + d;
	      float s5 = numer/denom;
	    	
	      smallestS = 500.0;
	      if(s1 >= 0.0){
		smallestS = s1;
		normal = vec3(-1.0,0.0,0.0);
	      }
	      if(s2 < smallestS && s2 >=0.0){
		normal = vec3(1.0,0.0,0.0);
		smallestS = s2;
	      }
	      if(s3 < smallestS && s3 >=0.0){
		normal = vec3(0.0,1.0,0.0);
		smallestS = s3;
	      }	
	      if(s4 < smallestS && s4 >=0.0){
		normal = vec3(0.0,-1.0,0.0);
		smallestS = s4;
	      }	   
	      if(s5 < smallestS && s5 >=0.0){
		normal = vec3(0.0,0.0,1.0);
		smallestS = s5;
	      }	 	
	      pI = smallestS*u + prevPos;
	      if(smallestS > -eps_S && smallestS < eps_S){
		r = normal;
	      }	
	      else{
		l = normalize(pI-prevPos);
		r = reflect(l,normal);
	      }
	      dis = distance(pI,vec3(pos));
		   	         
	      //Set the previous position to the point of intersection.
	      prevPos = pI;        
	      pos = vec4(pI+(dis*r),pos.a);	   	         
	      prmCurr = reflect(prmCurr,normal);
		   	       
	    } 
		       
	  }//while for reflection
	}    
	timeStepUsed = timeStepUsed + timeStepSim;// stores total time used
	timeStepRem = time_step - timeStepUsed; //stores time remaining in the time_step
	timeStepSim = timeStepRem; // time stepfor next iteration
	if(timeStepRem>epsilon)loopThrough=true;
	          
      }//if for moving particles if time is remaining 
	     
    }// if for domain
      
  }//while for loopThrough 
   
  if(pos.a <= 0.0 && (!(life_time <= 0.0))){
    gl_FragData[0] = vec4(100.0, 100.0, 100.0, life_time+1.0);
    gl_FragData[1] = vec4(prmCurr, 1.0);
  }
  else{
    gl_FragData[0] = pos;
    gl_FragData[1] = vec4(prmCurr, 1.0);
  }
   
}// Main
