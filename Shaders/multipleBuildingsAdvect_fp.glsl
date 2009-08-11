// The version numbers are required and since Apple supports version
// 1.20 of the GLSL spec, this is the version we currently support.
// By not supplying a version, the compiler assumes version 1.10.

#version 120

#extension GL_ARB_texture_rectangle : enable

uniform sampler2DRect pos_texunit;
uniform sampler2DRect primePrev;
uniform sampler2DRect wind_texunit;
uniform sampler2DRect lambda;
uniform sampler2DRect random;
uniform sampler2DRect tau_dz;
uniform sampler2DRect duvw_dz;
uniform sampler2DRect buildings;
uniform sampler2DRect cellType;

uniform int nx;
uniform int ny;
uniform int nz;
uniform float dx;
uniform float dy;
uniform float dz;
uniform int nxdx;
uniform int nydy;
uniform int nzdz;
uniform int numInRow;
uniform float time_step;
uniform float life_time;

uniform int color_advect_terms;
uniform int velocity_to_color_by;

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

int check; //to check how many times CFL condition got applied
//function returns the celltype of the location i,j,k
float ReturnCellType(int i,int j,int k)
{
          vec4 cell_type=vec4(1.0,1.0,1.0,1.0);
          vec2 cIndex;
          if(k < 0)
	    k = 0;

	  if(k >= 0){
	    cIndex.s = float(i) + float(mod(float(k),float(numInRow)))*float(nxdx);
	    cIndex.t = float(j) + float(floor(float(k)/float(numInRow)))*float(nydy);
	    cell_type = vec4(texture2DRect(cellType, cIndex));  
	  }
          return cell_type.x;
}
float absolute(float i)
{
      if(i<0) return -i;
      else return i;
}
int maximum(int i,int j)
{
     if(i>j) return i;
     else return j;
}
void main(void)
{
  vec4 poi=vec4(-1.0,-1.0,-1.0,-1.0);
  vec3 drift_term;
  vec3 memory_term;
  vec3 random_term;
  vec4 color = vec4(1.0,1.0,1.0,1.0);

  vec2 texCoord = gl_TexCoord[0].xy;
  vec3 prmPrev = vec3(texture2DRect(primePrev, texCoord));
  vec3 prmCurr = prmPrev;
   
  //This gets the position of the particle in 3D space.
  vec4 pos = vec4(texture2DRect(pos_texunit, texCoord));
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
  if (rTexCoord.s > float(random_texWidth))
    rTexCoord.s = rTexCoord.s - float(random_texWidth);
  if (rTexCoord.t > float(random_texHeight))
    rTexCoord.t = rTexCoord.t - float(random_texHeight);

  // lookup the random value to be used for this particle in
  // this timestep
  vec3 randn = vec3(texture2DRect(random, rTexCoord));

  //float xRandom = randn.x;
  //float yRandom = randn.y;
  //float zRandom = randn.z;     

  float upPrev=prmPrev.x;
  float vpPrev=prmPrev.y;
  float wpPrev=prmPrev.z;	

  //===============For CFL Condition===========================
  //float dx=1.0;  //defining grid size
  //float dy=1.0;
  //float dz=1.0;
  float epsilon=0.0001; //to avoid precision issues when checking how much time is remaining
   
  float timeStepRem = time_step; // time step remaining at this point
  float timeStepSim = time_step;  //time step to be used in the simulation
  float timeStepUsed = 0.0;   //used for storing how much time has been used from a time step during simulation
  bool loopThrough=true; 
   
  check=1;   

  if(!(life_time <= 0.0)){
    if(pos.a != life_time+1.0){
      //Decrement the alpha value based on the time_step
      pos.a = pos.a - time_step;

    }
  }   

  //vec2 index;
  //vec4 wind;
  //vec4 lam;
  //vec4 tau;
  //vec4 ddz;

  bool single = true;
  while(loopThrough){
    loopThrough = false;

    //The floor of the position in 3D space is needed to find the index into
    //the 2D Texture.
    int j = int(floor(pos.y));
    int i = int(floor(pos.x));
    int k = int(floor(pos.z));
	
    //This statement doesn't allow particles to move outside the domain.
    //So, the particles inside the domain are the only ones operated on. 
    if((i < nx) && (j < ny) && (k < nz) && (i >= 0) && (j >= 0) && (k >= 0)){
   
      
      //This is the initial lookup into the 2D texture that holds the wind field.
      vec2 index;
      index.s = float(i) + float(mod(float(k),float(numInRow)))*float(nxdx);
      index.t = float(j) + float(floor(float(k)/float(numInRow)))*float(nydy);

      vec4 wind = vec4(texture2DRect(wind_texunit, index));
      vec4 lam = vec4(texture2DRect(lambda, index));
      vec4 tau = vec4(texture2DRect(tau_dz, index));
      vec4 ddz = vec4(texture2DRect(duvw_dz, index));

      float dudz = ddz.x;
      float ustar = ddz.w; //temporarily used to get ustar to calculate sigmas in shader

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
       

      memory_term.x = -CoEps_D2*(Lam11*upPrev+Lam13*wpPrev)*timeStepSim;
      memory_term.y = -CoEps_D2*(Lam22*vpPrev)*timeStepSim;
      memory_term.z = -CoEps_D2*(Lam13*upPrev + Lam33*wpPrev)*timeStepSim;

      drift_term.x = (dudz*wpPrev + 0.5*Tau13)*timeStepSim + (Tau11*(Lam11*upPrev + Lam13*wpPrev) + Tau13*(Lam13*upPrev + Lam33*wpPrev))*(wpPrev/2.0)*timeStepSim;
      drift_term.y = (Tau22*Lam22*vpPrev*(wpPrev/2.0))*timeStepSim;
      drift_term.z = (0.5*Tau33)*timeStepSim + (Tau13*(Lam11*upPrev + Lam13*wpPrev) + Tau33*(Lam13*upPrev + Lam33*wpPrev))*
	(wpPrev/2.0)*timeStepSim;

      random_term.x = pow((2.0*CoEps_D2*timeStepSim),0.5)*randn.x;
      random_term.y = pow((2.0*CoEps_D2*timeStepSim),0.5)*randn.y;
      random_term.z = pow((2.0*CoEps_D2*timeStepSim),0.5)*randn.z;

      /*float du = (-CoEps_D2*(Lam11*upPrev+Lam13*wpPrev) + dudz*wpPrev + 0.5*Tau13)*timeStepSim
	+ (Tau11*(Lam11*upPrev + Lam13*wpPrev) + Tau13*(Lam13*upPrev + Lam33*wpPrev))*
	(wpPrev/2.0)*timeStepSim + pow((2.0*CoEps_D2*timeStepSim),0.5)*randn.x;
	
      float dv = (-CoEps_D2*(Lam22*vpPrev) + Tau22*Lam22*vpPrev*(wpPrev/2.0))*timeStepSim + 
	pow((2.0*CoEps_D2*timeStepSim),0.5)*randn.y;

      float dw = (-CoEps_D2*(Lam13*upPrev + Lam33*wpPrev) + 0.5*Tau33)*timeStepSim +
	(Tau13*(Lam11*upPrev + Lam13*wpPrev) + Tau33*(Lam13*upPrev + Lam33*wpPrev))*
	(wpPrev/2.0)*timeStepSim + pow((2.0*CoEps_D2*timeStepSim),0.5)*randn.z;*/
      float du = memory_term.x + drift_term.x + random_term.x;
      float dv = memory_term.y + drift_term.y + random_term.y;
      float dw = memory_term.z + drift_term.z + random_term.z;
      
      float totVel= pow((upPrev*upPrev+vpPrev*vpPrev+wpPrev*wpPrev),0.5);
      float totVelNew=pow( ( (upPrev+du)*(upPrev+du)+(vpPrev+dv)*(vpPrev+dv)+(wpPrev+dw)*(wpPrev+dw)),0.5);
      float totVelComp = totVel;
         
      if(totVelComp<totVelNew)
	totVelComp = totVelNew;
      
      if((totVelComp>2.0*sigU) && single){
         
         //generating random number
         rTexCoord.s = rTexCoord.s + 1.0;
         if(rTexCoord.s > float(random_texWidth))
	   rTexCoord.s = rTexCoord.s - float(random_texWidth);

         vec3 randnum = vec3(texture2DRect(random, rTexCoord));
            
         upPrev = sigU * randnum.x;
         vpPrev = sigV * randnum.y;
         wpPrev = sigW * randnum.z; 
         loopThrough=true;
         single=false;     
      }
      

      float disX = (wind.x+du)*timeStepSim; //calculating distance travelled by the particle
      float disY = (wind.y+dv)*timeStepSim;
      float disZ = (wind.z+dw)*timeStepSim;
         
      if(disX<0.0)disX=-disX;//calculating the absolute value
      if(disY<0.0)disY=-disY;
      if(disZ<0.0)disZ=-disZ;
         
      //if( (disX>0.7*dx) || (disY>0.7*dy) || (disZ>0.7*dz) ){  //CFL condition
      if( disX>1.4*dx || disY>1.4*dy || disZ>1.4*dz 
         || pow((CoEps_D2*timeStepSim),0.5)*timeStepSim > 1.4*dx 
         || pow((CoEps_D2*timeStepSim),0.5)*timeStepSim > 1.4*dy
         || pow((CoEps_D2*timeStepSim),0.5)*timeStepSim > 1.4*dz){  //CFL condition
	   

	timeStepSim = timeStepSim/2.0;
	loopThrough=true; //make loopThrough true so that it calculate distance travelled again by new timeStepSim
	check=check+1;
      }
      if(loopThrough == false){
	
	prmCurr = vec3(upPrev+du,vpPrev+dv,wpPrev+dw);

	//Now move the particle by adding the direction.
	pos = pos + vec4(wind.x,wind.y,wind.z,0.0)*timeStepSim + vec4(0.5*(prmPrev+prmCurr),0.0)*timeStepSim;
        //if(ReturnCellType(i,j,k)==0) pos=pos+vec4(1.0,-2.0,0.0,0.0);
        vec3 n;
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
        float d;
	float denom;
	float numer;
	
	vec2 cIndex;	

	j = int(floor(pos.y));
	i = int(floor(pos.x));
	k = int(floor(pos.z));
	int cnt = 0;	
	float eps_S = 0.0001;
	float eps_d = 0.01;
	float smallestS = 100.0;
        
        float s1 = -1.0;
	float s2 = -1.0;
        float s3 = -1.0;
	float s4 = -1.0;
	float s5 = -1.0;
	float s6 = -1.0;
	float s7 = -1.0;
         
        if(k<0)
        {
               u = vec3(pos) - prevPos;
               normal=vec3(0.0,0.0,1.0);
               n = vec3(0.0,0.0,1.0);
	       numer = dot(n,prevPos);
	       denom = dot(n,u);
	       float s1 = -1.0;
               s1=-numer/denom;
               pI = s1*u + prevPos;
	       if((s1 >= -0.0001) && (s1 <= 0.0001)){
	            pI = prevPos;
	            r = normal;
	       }	
	       else{
	            l = normalize(pI-prevPos);
	            r = normalize(reflect(l,normal));
	         }
	    dis = distance(pI,vec3(pos));		
	    pos = vec4(pI+(dis*r),pos.a);
	    prmCurr = reflect(prmCurr,normal);
        }       
	if((i < nx) && (j < ny) && (k < nz) && (i >= 0) && (j >= 0)){ 
          /*vec4 cell_type = vec4(1.0,1.0,1.0,1.0);
	            if(k < 0)
	               k = 0;
	            if(k >= 0){
	                cIndex.s = float(i) + float(mod(float(k),float(numInRow)))*float(nxdx);
	                cIndex.t = float(j) + float(floor(float(k)/float(numInRow)))*float(nydy);
	                cell_type = vec4(texture2DRect(cellType, cIndex));}*/
          //if(ReturnCellType(i,j,k)==0){pos.x=prevPos.x;pos.y=prevPos.y;pos.z=prevPos.z;}
          if(ReturnCellType(i,j,k)==0)
          {   
            //Reflection using the celltypes of the building 
            float isign,jsign,ksign;
            int imm,jmm,kmm;
            while((ReturnCellType(i,j,k)==0) && (cnt<1))
            {
            //imm=int(floor(pos.x)); 
            //jmm=int(floor(pos.y)); 
            //kmm=int(floor(pos.z));  

            u = vec3(pos)-prevPos;
            vec3 prevPos1;
            prevPos1=prevPos;
            //calculating S values for each face of the cell type.
                //-x normal  
	  	n = vec3(-1.0,0.0,0.0);
	   	d = -dot(n,vec3(i,j+0.5,k+0.5));
	   	denom = dot(n,u);
	    	numer = dot(n,prevPos) + d;
	  	s1 = -numer/denom;
      
	    	//+x normal
	    	n = vec3(1.0,0.0,0.0);
	    	d = -dot(n,vec3(i+1.0,j+0.5,k+0.5));
	    	denom = dot(n,u);
	    	numer = dot(n,prevPos) + d;
	    	s2 = -numer/denom;
        
	    	//+y normal
	    	n = vec3(0.0,1.0,0.0);
	    	d = -dot(n,vec3(i+0.5,j+1,k+0.5));
	    	denom = dot(n,u);
	    	numer = dot(n,prevPos) + d;
	    	s3 = -numer/denom;
           
	    	//-y normal
	    	n = vec3(0.0,-1.0,0.0);
	    	d = -dot(n,vec3(i+0.5,j,k+0.5));
	    	denom = dot(n,u);
	    	numer = dot(n,prevPos) + d;
	    	s4 = -numer/denom;
      	
	    	//+z normal
	    	n = vec3(0.0,0.0,1.0);
	    	d = -dot(n,vec3(i+0.5,j+0.5,k+1));
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
	    
	        /* //Detect Edge Collision
	        float edgeS = abs(smallestS-s7);
	    	if((edgeS < eps_d)){
	      	  //smallestS = s6;
	      	  normal = normalize(normal+vec3(0.0,0.0,1.0));
	    	}
	    	else */if((s7 < smallestS) && (s7 >= -eps_S)){
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
      	    poi=vec4(pI,2.0);
	    prevPos = pI;
	    pos = vec4(pI+(dis*r),pos.a);
	    prmCurr = reflect(prmCurr,normal);
      	      
	    j = int(floor(pos.y));
	    i = int(floor(pos.x));
	    k = int(floor(pos.z));
            if(ReturnCellType(i,j,k)==0){pos.x=prevPos1.x; pos.y=prevPos1.y; pos.z=prevPos1.z;}

	    //NOTE: Consider what happens if building is too close to domain.
	    //Do check to make sure i,j,k's are valid;
	    /*vec4 cell_type = vec4(1.0,1.0,1.0,1.0);
	    if(k < 0)
	      k = 0;
	    if(k >= 0){
	      cIndex.s = float(i) + float(mod(float(k),float(numInRow)))*float(nxdx);
	      cIndex.t = float(j) + float(floor(float(k)/float(numInRow)))*float(nydy);
	      //cIndex.s = j + mod(k,numInRow)*nx;
	      //cIndex.t = i + int(floor(float(k)/float(numInRow)))*ny;
	      cell_type = vec4(texture2DRect(cellType, cIndex));
	    }*/
            cnt++;
           }
          }//end reflection-if loop
	  timeStepUsed = timeStepUsed + timeStepSim;// stores total time used
	  timeStepRem = time_step - timeStepUsed; //stores time remaining in the time_step
	  timeStepSim = timeStepRem; // time stepfor next iteration

	  if(timeStepRem>epsilon){
	    loopThrough=true;
	    //I think we need to do this????
	    //update prevPos and prevPrime
	    prevPos = vec3(pos);
	    prmPrev = prmCurr;
	    upPrev=prmPrev.x;
	    vpPrev=prmPrev.y;
	    wpPrev=prmPrev.z;
	    //generating random number
	    rTexCoord.s = rTexCoord.s + 1.0;
	    rTexCoord.t = rTexCoord.t + 1.0;
	    if (rTexCoord.s > float(random_texWidth))
	      rTexCoord.s = rTexCoord.s - float(random_texWidth);
	    if (rTexCoord.t > float(random_texHeight))
	      rTexCoord.t = rTexCoord.t - float(random_texHeight);

	    randn = vec3(texture2DRect(random, rTexCoord));
	    //check = 1;
	  }

	}//make sure particle is in domain still
	//if it isn't, don't perform any more advection
	else{
	  loopThrough = false;
	}	

      }//loopThrough if
     
      if(color_advect_terms == 1){
	//Find largest advect term and set color
	color = vec4(1.0,0.0,0.0,1.0);

	//float me = memory_term.x;//length(memory_term);
	//float dr = drift_term.x;//length(drift_term);
	//float ra = random_term.x;//length(random_term);

	float me = length(memory_term);
	float dr = length(drift_term);
	float ra = length(random_term);

	float largest = me;

	if(dr > largest){
	  largest = dr;
	  color = vec4(0.0,1.0,0.0,1.0);
	}
	if(ra > largest){
	  color = vec4(0.0,0.0,1.0,1.0);
	}
      }
      
    }//if on domain check

  }//while loopthrough condition
  
 
  
  if(color_advect_terms == 1){
    color=poi;   
    if(pos.a <= 0.0 && (!(life_time <= 0.0))){
      gl_FragData[0] = vec4(100.0, 100.0, 100.0, life_time+1.0);
      gl_FragData[1] = vec4(prmCurr, 1.0);
      gl_FragData[2] = color;
    }
    else{
      gl_FragData[0] = pos;
      gl_FragData[1] = vec4(prmCurr, 1.0);
      gl_FragData[2] = color;
    }

  }
  else{
    if(pos.a <= 0.0 && (!(life_time <= 0.0))){
      gl_FragData[0] = vec4(100.0, 100.0, 100.0, life_time+1.0);
      gl_FragData[1] = vec4(prmCurr, 1.0);
    }
    else{
      gl_FragData[0] = pos;
      gl_FragData[1] = vec4(prmCurr, 1.0);
    }
  }
   
}


//CFL Constant=1.4
//Other constant 2.0
