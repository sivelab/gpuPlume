// The version numbers are required and since Apple supports version
// 1.20 of the GLSL spec, this is the version we currently support.
// By not supplying a version, the compiler assumes version 1.10.

// #version 120

#extension GL_ARB_texture_rectangle : enable

uniform sampler2DRect pos_texunit;
uniform sampler2DRect primePrev;
uniform sampler2DRect wind_texunit;
uniform sampler2DRect lambda;
uniform sampler2DRect random;
uniform sampler2DRect tau_dz;
uniform sampler2DRect duvw_dz;
uniform sampler2DRect dxy_wall;
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

//Balli: Declaring variables here globally to check their values in prime texture output in "dumpcontents"- may be removed later
float dwall,coeps,dzm,dym,dyp,dxm,dzc,dsigwdn,du,dv,dw;
float alph1ij,alph2ij,alph3ij,bet1ij,bet2ij,bet3ij,gam1ij,gam2ij,gam3ij,alphn1ij,alphn2ij,alphn3ij,betn1ij,betn2ij;
float dfzm,drzp,dwxm,dwxp,dwym,dwyp,u1,v1,w1;
float dutotdni,dutotdsi,ani,bni,cni,dutotds;
float cosomeg,sinomeg,tanomeg,omeg1,cosomeg1,sinomeg1,dutotdn1,dutotdn;

//Hard wired variables !!MUST BE READ INTO SHADER FROM PARTICLE CONTROL  ***!!! ATTENTION!!!***
float h=25.;
float rcl=0.;
float z0 = 0.1;
float z0coeff = 1.01;
float taylor_flag = 0.;
//Hardwired variables end

float xnu = 0.;//Balli: initialized as zero; can change later if we add some more capabilities from QP
float kkar = 0.4;

float z0fac = z0*z0coeff;

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

// int maximum(int i,int j)
// {
//   if(i>j) return i;
//   else return j;
// }

float sign(const float A,const float B){
  float R;
  if(B>=0.f)
    R=absolute(A);
  else
    R=-absolute(A);
  return R;	
}

void main(void)
{
  vec4 poi=vec4(-1.0,-1.0,-1.0,-1.0);
  // vec3 drift_term;
  // vec3 memory_term;
  // vec3 random_term;
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
    
  int k=0;
  float sigu,sigv,sigw,ustar,tau11,tau22,tau33,tau13,tau13sq,lam11,lam22,lam33,lam13;
  int lthrough=0;
  bool single = true;
  while(loopThrough){
    lthrough = lthrough + 1;
    loopThrough = false;

    //The floor of the position in 3D space is needed to find the index into
    //the 2D Texture.
    int j = int(floor(pos.y));
    int i = int(floor(pos.x));
    int k = int(floor(pos.z));

    float dzc = pos.z - (k +0.5*dz);

    float dxc = pos.x - .5*dx-dx*i;
    float dyc = pos.y - .5*dy-dy*j;

    //This statement doesn't allow particles to move outside the domain.
    //So, the particles inside the domain are the only ones operated on. 

    if((i < nx) && (j < ny) && (k < nz) && (i >= 0) && (j >= 0) && (k >= 0)){
      //This is the initial lookup into the 2D texture that holds the wind field.

      vec2 index;
      index.s = float(i) + float(mod(float(k),float(numInRow)))*float(nxdx);
      index.t = float(j) + float(floor(float(k)/float(numInRow)))*float(nydy);

      vec4 wind = vec4(texture2DRect(wind_texunit, index));

      vec4 dxyz = vec4(texture2DRect(dxy_wall, index));
      float dzm = dxyz.x;
      float dzp = dxyz.y;
      float dym = dxyz.z;
      float dyp = dxyz.w;
      //Balli: Name of the following texture can be changed to better reflect its contents
      vec4 dxyz2 = vec4(texture2DRect(tau_dz, index));
      float dxp = dxyz2.x;
      float dxm = dxyz2.y;
            
      dfzm = dzm + dzc;
      drzp = dzp - dzc;
      dwxm = dxm + dxc;
      dwxp = dxp - dxc;
      dwym = dym + dyc;
      dwyp = dyp - dyc;

      //Balli: if the particle is too close to a wall
      if(dfzm<z0fac){
	pos.z = pos.z+z0fac;
	dzc   = pos.z - (k +0.5*dz);
	dfzm  = dzm + dzc;
      }
      if(dwxp<z0fac){
	pos.x=pos.x-z0fac;
	dxc=pos.x-.5*dx-dx*float(i);
	dwxp=dxp-dxc;
      }
      if(dwxm<z0fac){
	pos.x=pos.x+z0fac;
	dxc=pos.x-.5*dx-dx*float(i);
	dwxm=dxm + dxc;
      }
      if(dwyp<z0fac){
	pos.y=pos.y-z0fac;
	dyc=pos.y-.5*dy-dy*float(j);
	dwyp=dyp-dyc;
      }
      if(dwym<z0fac){
	pos.y=pos.y+z0fac;
	dyc=pos.y-.5*dy-dy*float(j);
	dwym=dym+dyc;
      }
      if(drzp<z0*z0fac){
	pos.z=pos.z-z0fac;
	dzc=pos.z-(k+0.5*dz);
	drzp=dzp-dzc;
      }
            
      //Balli: Name of the following texture can be changed to better reflect its contents
      vec4 sig_ustar = vec4(texture2DRect(lambda, index));
      sigu  = sig_ustar.x;
      sigv  = sig_ustar.y;
      sigw  = sig_ustar.z;
      ustar = sig_ustar.w;
            
      vec4 dutot = vec4(texture2DRect(duvw_dz, index));
      float dutotdxi = dutot.x;
      float dutotdyi = dutot.y;
      float dutotdzi = dutot.z;
      dsigwdn = dutot.w;
            
      tau11=sigu*sigu;
      tau22=sigv*sigv;
      tau33=sigw*sigw;
            
      //Balli: finding minimum distance to a wall
      float temp  = min(dfzm,drzp);
      temp        = min(temp,dwxm);
      temp        = min(temp,dwxp);
      temp        = min(temp,dwym);
      dwall       = min(temp,dwyp);

      float ustar3 = ustar*ustar*ustar;
            
      /*
	Balli: Following is a copy of the fortran code for "detang"
	subroutine for obtaining rotational parameters for the
	advection process.  It will be better if we can move this
	above as a function so that advection code looks clean and
	simple
      */

      float UV       = wind.x*wind.x + wind.y*wind.y;
      float UVW      = wind.z*wind.z + UV;
      // float tVel     = pow(UVW,0.5);
      float tVel     = sqrt(UVW);
      float tVelI    = 1./(tVel+1.E-10);
      // float totUV    = pow(UV,0.5);
      float totUV    = sqrt(UV);
      float totUVI   = 1./totUV +1.E-10;
            
      int iomega = 0;
      float cospsi,sinpsi,cosphiw,sinphiw,omenum,omeden;
      float omeg2,cosomeg2,sinomeg2,dutotdn2,omeg;
      float pi=acos(-1.);
      float betn3ij,gamn1ij,gamn2ij,gamn3ij;

      float ddxddy2,ddxddyddz2,ddxddy,ddxddyddz;
            
      if(totUV>1.e-05){
	cospsi  = wind.x*totUVI;
	sinpsi  = wind.y*totUV;
	sinphiw = wind.z*tVelI;
	cosphiw = totUV*tVelI;
	omenum  = -dutotdxi*sinpsi + dutotdyi*cospsi;
	omeden  = dutotdxi*cospsi*sinphiw + dutotdyi*sinpsi*sinphiw - dutotdzi*cosphiw;
	if(iomega==0){
	  if(absolute(omeden)<1.e-10){
	    cosomeg = 0.;
	    sinomeg = 1.;
	    dutotdn = dutotdxi*(sinpsi*sinomeg-cospsi*sinphiw*cosomeg) 
	      -dutotdyi*(cospsi*sinomeg+sinpsi*sinphiw*cosomeg) 
	      +dutotdzi*cosomeg*cosphiw;
	    if(dutotdn<0.)sinomeg = -1.;
	  }
	  else{
	    tanomeg  = omenum/omeden;
	    omeg1    = atan(tanomeg);
	    cosomeg1 = cos(omeg1);
	    sinomeg1 = sin(omeg1);
	    dutotdn1 = dutotdxi*(sinpsi*sinomeg1-cospsi*sinphiw*cosomeg1) 
	      -dutotdyi*(cospsi*sinomeg1+sinpsi*sinphiw*cosomeg1) 
	      +dutotdzi*cosomeg1*cosphiw;
	    omeg2    = omeg1+pi;
	    cosomeg2 = cos(omeg2);
	    sinomeg2 = sin(omeg2);
	    dutotdn2 = dutotdxi*(sinpsi*sinomeg2-cospsi*sinphiw*cosomeg2) 
	      -dutotdyi*(cospsi*sinomeg2+sinpsi*sinphiw*cosomeg2) 
	      +dutotdzi*cosomeg2*cosphiw;
	    if(dutotdn2>dutotdn1){
	      dutotdn = dutotdn2;
	      omeg    = omeg2;
	      cosomeg = cosomeg2;
	      sinomeg = sinomeg2;
	    }
	    else{
	      dutotdn = dutotdn1;
	      omeg    = omeg1;
	      cosomeg = cosomeg1;
	      sinomeg = sinomeg1;
	    }
	  }
	}
                
	else{
	  if(iomega==1){
	    omeg    = 0.;
	    cosomeg = 1.;
	    sinomeg = 0.;
	  }
	  else{
	    if(absolute(cospsi)>0.5){
	      omeg    = 0.5*pi;
	      cosomeg = 0.;
	      sinomeg = -sign(cospsi,1.0);
	      if(iomega==3){
		sinomeg = sign(cospsi,1.0);
		omeg    = -omeg;
	      }
	    }
	    /*else{
	      return;
	      }*/
	  }
	}
                
	alph1ij  = cospsi*cosphiw;
	alph2ij  = -sinpsi*cosomeg-cospsi*sinphiw*sinomeg;
	alph3ij  = sinpsi*sinomeg-cospsi*sinphiw*cosomeg;
	bet1ij   = sinpsi*cosphiw;
	bet2ij   = cospsi*cosomeg-sinpsi*sinphiw*sinomeg;
	bet3ij   = -cospsi*sinomeg-sinpsi*sinphiw*cosomeg;
	gam1ij   = sinphiw;
	gam2ij   = cosphiw*sinomeg;
	gam3ij   = cosphiw*cosomeg;
	alphn1ij = cospsi*cosphiw;
	alphn2ij = sinpsi*cosphiw;
	alphn3ij = sinphiw;
	betn1ij  = -sinpsi*cosomeg-cospsi*sinphiw*sinomeg;
	betn2ij  = cospsi*cosomeg-sinpsi*sinphiw*sinomeg;
	betn3ij  = cosphiw*sinomeg;
	gamn1ij  = sinpsi*sinomeg-cospsi*sinphiw*cosomeg;
	gamn2ij  = -cospsi*sinomeg-sinpsi*sinphiw*cosomeg;
	gamn3ij  = cosphiw*cosomeg;
	dutotdn  = dutotdxi*(sinpsi*sinomeg-cospsi*sinphiw*cosomeg) 
	  -dutotdyi*(cospsi*sinomeg+sinpsi*sinphiw*cosomeg) 
	  +dutotdzi*cosomeg*cosphiw;
	dutotds  = dutotdxi*cospsi*cosphiw+dutotdyi*
	  sinpsi*cosphiw+dutotdzi*sinphiw;
	dutotdni = dutotdn;
	dutotdsi = dutotds;
                
	ani      = (sinpsi*sinomeg-cospsi*sinphiw*cosomeg);
	bni      = -(cospsi*sinomeg+sinpsi*sinphiw*cosomeg);
	cni      = cosomeg*cosphiw;
      }
      else{
	if(absolute(wind.z)<1.e-05){
	  ddxddy2    = dutotdxi*dutotdxi+dutotdyi*dutotdyi;
	  ddxddyddz2 = ddxddy+dutotdzi*dutotdzi;
	  //	  ddxddy     = pow(ddxddy2,0.5);
	  //      ddxddyddz  = pow(ddxddyddz2,0.5);
	  ddxddy     = sqrt(ddxddy2);
	  ddxddyddz  = sqrt(ddxddyddz2);
	  float ddxddyddzI=1./ddxddyddz;
	  if(ddxddy>0.){
	    cospsi = dutotdxi/ddxddy;
	    sinpsi = dutotdyi/ddxddy;
	  }
	  else{
	    cospsi = 1.;
	    sinpsi = 0.;
	  }
	  if(ddxddyddz>0.){
	    cosphiw=dutotdzi*ddxddyddzI;
	    sinphiw=ddxddy*ddxddyddzI;
	  }
	  else{
	    cosphiw=1.;
	    sinphiw=0.;
	  }
	  alphn1ij = cospsi*cosphiw;
	  alphn2ij = sinpsi*cosphiw;
	  alphn3ij = -sinphiw;
	  betn1ij  = -sinpsi;
	  betn2ij  = cospsi;
	  betn3ij  = 0.;
	  gamn1ij  = sinphiw*cospsi;
	  gamn2ij  = sinphiw*sinpsi;
	  gamn3ij  = cosphiw;
	  alph1ij  = cospsi*cosphiw;
	  alph2ij  = -sinpsi;
	  alph3ij  = sinphiw*cospsi;
	  bet1ij   = sinpsi*cosphiw;
	  bet2ij   = cospsi;
	  bet3ij   = sinphiw*sinpsi;
	  gam1ij   = -sinphiw;
	  gam2ij   = 0.;
	  gam3ij   = cosphiw;
	  dutotdni = ddxddyddz;
	  dutotdsi = 0.;
	  ani      = sinphiw*cospsi;
	  bni      = sinphiw*sinpsi;
	  cni      = cosphiw;
	}
	else{
	  ddxddy2       = dutotdxi*dutotdxi+dutotdyi*dutotdyi;
	  ddxddyddz2    = ddxddy+dutotdzi*dutotdzi;
	  ddxddy        = sqrt(ddxddy2);
	  ddxddyddz     = sqrt(ddxddyddz2);
	  // ddxddy        = pow(ddxddy2,0.5);
	  // ddxddyddz     = pow(ddxddyddz2,0.5);
	  float ddxddyI = 1./ddxddy;;
	  if(wind.z>0.){
	    if(ddxddy>0.){
	      cospsi=dutotdxi*ddxddyI;
	      sinpsi=dutotdyi*ddxddyI;
	    }
	    else{
	      cospsi=1.;
	      sinpsi=0.;
	    }
	    alph1ij=0.;
	    alph2ij=sinpsi;
	    alph3ij=cospsi;
	    bet1ij=0.;
	    bet2ij=-cospsi;
	    bet3ij=sinpsi;
	    gam1ij=1.;
	    gam2ij=0.;
	    gam3ij=0.;
	    alphn1ij=0.;
	    alphn2ij=0.;
	    alphn3ij=1.;
	    betn1ij=sinpsi;
	    betn2ij=-cospsi;
	    betn3ij=0.;
	    gamn1ij=cospsi;
	    gamn2ij=sinpsi;
	    gamn3ij=0.;
	    // dutotdni=pow(dutotdxi*dutotdxi+dutotdyi*dutotdyi,0.5);
	    dutotdni = sqrt(dutotdxi*dutotdxi+dutotdyi*dutotdyi);
	    dutotdni=dutotdni;
	    if(dutotdni<1.e-12)dutotdni=1.e-12;
	    dutotdsi=dutotdzi;
	    ani=cospsi;
	    bni=sinpsi;
	    cni=0.;
	  }
	  else{
	    if(ddxddy>0.){
	      cospsi=dutotdxi*ddxddyI;
	      sinpsi=dutotdyi*ddxddyI;
	    }
	    else{
	      cospsi=1.;
	      sinpsi=0.;
	    }
	    alphn1ij=0.;
	    alphn2ij=0.;
	    alphn3ij=-1.;
	    betn1ij=-sinpsi;
	    betn2ij=cospsi;
	    betn3ij=0.;
	    gamn1ij=cospsi;
	    gamn2ij=sinpsi;
	    gamn3ij=0.;
	    alph1ij=0.;
	    alph2ij=-sinpsi;
	    alph3ij=cospsi;
	    bet1ij=0.;
	    bet2ij=cospsi;
	    bet3ij=sinpsi;
	    gam1ij=-1.;
	    gam2ij=0.;
	    gam3ij=0.;
	    dutotdni=ddxddy;
	    dutotdsi=-dutotdzi;
	  }
	}
      }

      bool loopdt = true;
      int ivrelch = 0;
      float dt,utot,taylor_microscale,tls,dudet,dvdet,dwdet,duran,dvran,dwran;
      int counter =0;
      while(loopdt){
	counter = counter + 1;//Balli : Can be used in above while loop to avoid infinite loop
	loopdt = false;
	dt      = timeStepSim;
            
	u1 = upPrev * alphn1ij + vpPrev * alphn2ij + wpPrev * alphn3ij;
	v1 = upPrev * betn1ij  + vpPrev * betn2ij  + wpPrev * betn3ij ;
	w1 = upPrev * gamn1ij  + vpPrev * gamn2ij  + wpPrev * gamn3ij ; 
            
	if(rcl<0 && pos.z<=.99*h){
	  coeps=5.7*(ustar3*(1.-.75*pos.z*rcl)*pow( (1.-.85*pos.z/h),(1.5))/(dwall*kkar));
	}
	else{
	  if(rcl<0.){
	    coeps=5.7*(ustar3*(1.-.75*.99*h*rcl)*pow( (1.-.85*.99),(1.5))/(dwall*kkar));
	  }
	  else{
	    if(pos.z<=.99*h){
	      coeps=5.7*(ustar3*(1.+3.7*pos.z*rcl)*pow( (1.-.85*pos.z/h),(1.5))/(dwall*kkar));
	    }
	    else{
	      coeps=5.7*(ustar3*(1.+3.7*.99*h*rcl)*pow( (1.-.85*.99),(1.5))/(dwall*kkar));
	    }
	  }
	} 
	utot    = tVel;
	taylor_microscale=0.;
                
	if(taylor_flag > 0){
	  // taylor_microscale=pow( (8.55e-3)*xnu*sigu*sigu/(coeps*utot*utot) ,0.5);
	  taylor_microscale = sqrt( (8.55e-3)*xnu*sigu*sigu/(coeps*utot*utot) );
	}
                
	tls=2.*sigw*sigw/(coeps+1.e-36);
	if(dt > 0.5*tls &&  dt > taylor_microscale){
	  //dt=0.5*tls; // Commented this out as the code freezes if I uncomment this ***!!!ATTENTION!!!***
	  timeStepSim = dt;
	}
	if((k+0.5*dz)<=.99*h){
	  tau13=ustar*ustar*pow( (1.-(k+0.5*dz)/h),(1.5));
	}
	else{
	  tau13=ustar*ustar*pow((1.-.99),(1.5));
	}
	lam11 = 1./(tau11-tau13*tau13/tau33);
	lam22 = 1./tau22;
	lam33 = 1./(tau33-tau13*tau13/tau11);
	lam13 = -tau13/(tau11*tau33-tau13*tau13);
                
	dutotdn  = dutotdni;
	dutotds  = dutotdsi;
                
	rTexCoord.s = rTexCoord.s + 1.0;
	rTexCoord.t = rTexCoord.t + 1.0;
	if (rTexCoord.s > float(random_texWidth))
	  rTexCoord.s = rTexCoord.s - float(random_texWidth);
	if (rTexCoord.t > float(random_texHeight))
	  rTexCoord.t = rTexCoord.t - float(random_texHeight);
                
	randn = vec3(texture2DRect(random, rTexCoord));
                
	dudet = -.5*coeps*(lam11*u1+lam13*w1)*dt+dutotdn*w1*dt 
	  +utot*dutotds*dt+u1*dutotds*dt+sigu/(1.3*2.)*dsigwdn*dt 
	  +sigu*dsigwdn
	  *(u1*(lam11*5./1.3+lam13/(1.25*1.3))+w1*(lam13*4./1.3 +lam33/1.3))
	  *(0.5*w1)*dt-w1*dutotdn*dt;
                
	// duran=pow(coeps,0.5)*randn.x*pow(dt,0.5);
	duran = sqrt(coeps) * randn.x * sqrt(dt);
                
	dvdet=-.5*coeps*lam22*v1*dt+(2./1.3)*sigv*dsigwdn*lam22*v1*w1*dt;
                
	// dvran=pow(coeps,.5)*randn.y*pow(dt,0.5);
	dvran = sqrt(coeps) * randn.y * sqrt(dt);
                
	dwdet=-.5*coeps*(lam13*u1+lam33*w1)*dt 
	  +sigw*dsigwdn*dt+((lam13+lam11/1.69)*u1+(lam13/1.69+lam33)*w1)*w1*sigw*dsigwdn*dt;
                
	// dwran=pow(coeps,.5)*randn.z*pow(dt,.5);
	dwran = sqrt(coeps) * randn.z * sqrt(dt);

	// float vrel2   = pow(((u1+dudet)*(u1+dudet)+(v1+dvdet)*(v1+dvdet)+(w1+dwdet)*(w1+dwdet)),.5);
	// float vrel1   = pow((u1*u1+v1*v1+w1*w1),.5);
	float vrel2   = sqrt(((u1+dudet)*(u1+dudet)+(v1+dvdet)*(v1+dvdet)+(w1+dwdet)*(w1+dwdet)));
	float vrel1   = sqrt((u1*u1+v1*v1+w1*w1));
	float vrel    = max(vrel2,vrel1);

                
	if(vrel>4.*sigu && ivrelch==0||vrel>100.){
	  rTexCoord.s = rTexCoord.s + 1.0;
	  rTexCoord.t = rTexCoord.t + 1.0;
	  if (rTexCoord.s > float(random_texWidth))
	    rTexCoord.s = rTexCoord.s - float(random_texWidth);
	  if (rTexCoord.t > float(random_texHeight))
	    rTexCoord.t = rTexCoord.t - float(random_texHeight);
                    
	  randn = vec3(texture2DRect(random, rTexCoord));
                    
	  u1      = sigu*randn.x;
	  v1      = sigv*randn.y;
	  w1      = sigw*randn.z;
	  upPrev  = u1 * alph1ij + v1 * alph2ij + w1 * alph3ij ;
	  vpPrev  = u1 * bet1ij  + v1 * bet2ij  + w1 * bet3ij  ;
	  wpPrev  = u1 * gam1ij  + v1 * gam2ij  + w1 * gam3ij  ;
	  ivrelch = 1;
	  loopdt=true;
	}
	float xrel = vrel*dt + pow( (wind.x*wind.x+wind.y*wind.y+wind.z*wind.z) , 0.5)*dt;
	if(xrel+absolute(wind.x)*dt>0.7*dx||5.*pow((coeps*dt),0.5)*dt>0.7*dx){
	  timeStepSim=min(dt/2.,timeStepRem);
	  loopdt=true;
	}
	if(xrel+absolute(wind.y)*dt>0.7*dy||5.*pow((coeps*dt),.5)*dt>0.7*dy){
	  timeStepSim=min(dt/2.,timeStepRem);
	  loopdt=true;
	}

	if(xrel+absolute(wind.z)*dt>0.7*dz||5.*pow( (coeps*dt),.5)*dt>0.7*dz){
	  timeStepSim=min(dt/2.,timeStepRem);
	  loopdt=true;
	}
	else{
	  loopdt=false;
	}
                
	if(timeStepSim<=epsilon){
	  loopdt = false;
	}
      }//while loopdt ends
      if(loopdt == false){
	du = dudet + duran;
	dv = dvdet + dvran;
	dw = dwdet + dwran;
                
	float uNewRan = (u1+du)*alph1ij + (v1+dv)*alph2ij + (w1+dw)*alph3ij ;
	float vNewRan = (u1+du)*bet1ij  + (v1+dv)*bet2ij  + (w1+dw)*bet3ij  ;
	float wNewRan = (u1+du)*gam1ij  + (v1+dv)*gam2ij  + (w1+dw)*gam3ij  ;

	prmCurr = vec3(uNewRan,vNewRan,wNewRan);

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

	  // This was in the code and seems redundant with the while loop below...
	  // if(ReturnCellType(i,j,k)==0)
	  // {   

	      // 
	      // Reflection using the celltypes of the building 
	      // 

	      // float isign,jsign,ksign;
	      // int imm,jmm,kmm;

	      // why is this bounded to cnt < 1 ?? -Pete  -- when this is set above 1, particles spread heavily 
	      // against the main flow
	      while((ReturnCellType(i,j,k)==0) && (cnt<2))
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
		     else */ if((s7 < smallestS) && (s7 >= -eps_S)){
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

		  i = int(floor(pos.x));
		  j = int(floor(pos.y));
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
	      // }//end reflection-if loop

	  timeStepUsed = timeStepUsed + timeStepSim;// stores total time used
	  timeStepRem = time_step - timeStepUsed; //stores time remaining in the time_step
	  timeStepSim = timeStepRem; // time stepfor next iteration

	  if(timeStepRem>epsilon){
	    loopThrough=true;
	    //I think we need to do this????
	    //update prevPos and prevPrime
	    prevPos = vec3(pos);
	    prmPrev = prmCurr;
	    // Commented out in bally's changes... not in prev version of LM
	    //upPrev=prmPrev.x;
	    //vpPrev=prmPrev.y;
	    //wpPrev=prmPrev.z;
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

      // if(color_advect_terms == 1){
	//Find largest advect term and set color
      //color = vec4(1.0,0.0,0.0,1.0);

	//float me = memory_term.x;//length(memory_term);
	//float dr = drift_term.x;//length(drift_term);
	//float ra = random_term.x;//length(random_term);

      //	float me = length(memory_term);
      //	float dr = length(drift_term);
      //	float ra = length(random_term);

      //	float largest = me;

      //	if(dr > largest){
      //	  largest = dr;
      //	  color = vec4(0.0,1.0,0.0,1.0);
      //	}
      //	if(ra > largest){
      //	  color = vec4(0.0,0.0,1.0,1.0);
      //        }
      // 

    }//if on domain check

  }//while loopthrough condition

  //  if(color_advect_terms == 1){
  //    color=poi;   
  //    if(pos.a <= 0.0 && (!(life_time <= 0.0))){
  //      gl_FragData[0] = vec4(100.0, 100.0, 100.0, life_time+1.0);
  //      gl_FragData[1] = vec4(prmCurr, 1.0);
  //      gl_FragData[2] = color;
  //    }
  //    else{
  //      gl_FragData[0] = pos;
  //      gl_FragData[1] = vec4(prmCurr, 1.0);
  //      gl_FragData[2] = color;
  //    }
  //  }
  //  else{
    if(pos.a <= 0.0 && (!(life_time <= 0.0))){
      gl_FragData[0] = vec4(100.0, 100.0, 100.0, life_time+1.0);
      gl_FragData[1] = vec4(prmCurr, 1.0);
    }
    else{
      gl_FragData[0] = pos;
      gl_FragData[1] = vec4(prmCurr, 1.0);
    }
    // }
}

//CFL Constant=1.4
//Other constant 2.0
