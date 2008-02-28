#include "particleEmitter.h"
#include <math.h>
#include <iostream>

ParticleEmitter::~ParticleEmitter(){

}
void ParticleEmitter::setVertices(){

}
int ParticleEmitter::EmitParticle(bool odd,GLuint pos0,GLuint pos1,
				  float time_step,GLuint prime0,GLuint prime1){
  int p_index;
  //float x;
  //float y;
  //float z;
  //float l;
  
    
 
  if(!indices->empty()){
    //THIS Method *seems* to work now!
    //Punch Hole method. Need to set drawbuffer and activate shader.
    if(Punch_Hole){

      if(odd)
	glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
      else 
	glDrawBuffer(GL_COLOR_ATTACHMENT1_EXT);

      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      gluOrtho2D(0, twidth, 0, theight);
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();     

    }
    	
    //Do this for each particle that is being emitted.
    for(int i = 0; i < numToEmit; i++){
      if(!indices->empty()){
	
      //First get available index
	p_index = indices->back();
	indices->pop_back();
	if(reuse){
	  pIndex newIndex;
	  newIndex.id = p_index;
	  newIndex.time = 0;
	  indicesInUse->push_back(newIndex);	   
	}	 

	//Determine the coordinates into the position texture
	s = (p_index%twidth);
	t = (p_index/twidth);	  

	/*
	//Get position value
	if(posCoord.empty())
	  std::cout << "Error: PosCoord empty" << std::endl;

	l = posCoord.back();
	posCoord.pop_back();
	z = posCoord.back();
	posCoord.pop_back();
	y = posCoord.back();
	posCoord.pop_back();
	x = posCoord.back();
	posCoord.pop_back();*/

	if(Punch_Hole){
	  glPointSize(1.0);

	  shader->activate();

       	  glViewport(s,t,1,1);
	  //std::cout << "particle num= " << p_index << "  s = " << s << "  t = " << t << std::endl;
	  glBegin(GL_POINTS);
	  {
	    glColor4f(xpos, ypos, zpos, lifeTime);
	    glVertex2f(0.5, 0.5);
	    //glVertex2f(s,t);
	  }
	  glEnd();
	  shader->deactivate();

	}
	else{
	  // Second mechanism to release particles.  Uses a texture
	  // copy which may likely be more expensive even though we're
	  // doing 1x1 pixels (but many times).  This operation
	  // appears to work consistently.
	  GLfloat value[4];
	  value[0] = xpos;
	  value[1] = ypos;
	  value[2] = zpos;
	  value[3] = lifeTime;

	  // write there via a glTexSubImage2D
	  if(odd)
	    // will read from this texture next
	    glBindTexture(GL_TEXTURE_RECTANGLE_ARB, pos0);
	  else 
	    // will read from this texture next
	    glBindTexture(GL_TEXTURE_RECTANGLE_ARB, pos1);

	  glTexSubImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, s, t, 1, 1, GL_RGBA, GL_FLOAT, value);
	}
  
      }
 
    }

    if(Punch_Hole){

      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      gluOrtho2D(-1, 1, -1, 1);
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();

    }
  
  }

  return numToEmit;

}
void ParticleEmitter::Draw(){

  glPointSize(4.0);
  glColor4f(0.0, 0.0, 1.0, 1.0);
  glBegin(GL_POINTS);
  {
    glVertex3f(xpos,ypos,zpos);
  }
  glEnd();

}
void ParticleEmitter::setPosition(float x, float y, float z){
  xpos = x;
  ypos = y;
  zpos = z;

}
void ParticleEmitter::getPosition(float* x, float*y, float*z){
  *x = xpos;
  *y = ypos;
  *z = zpos;
  //*r = 1.0;
}
void ParticleEmitter::getReleasedPosition(float*x,float*y,float*z){
  *x = xpos;
  *y = ypos;
  *z = zpos;
}
void ParticleEmitter::getIndex(int* x, int*y){
  *x = s;
  *y = t;
}

void ParticleEmitter::setParticleReuse(std::list<pIndex>* ind, float time){
  indicesInUse = ind;
 
  lifeTime = time;
  reuse = true;

}
void ParticleEmitter::setNumToEmit(float num){
  emitAmount = num;
}

//sets the value for numToEmit based on emission type
void ParticleEmitter::setEmitAmount(float time_step){
  switch(releaseType){
  case perSecond:
    if(!timeToEmit(time_step))
	numToEmit = 0;
    break;
  case onePerKeyPress:
    numToEmit = 1;
    emit = false;
    break;
  case instantaneous:
    numToEmit = (twidth*theight);
    emit = false;
    break;
  default:
    numToEmit = (int)emitAmount;
    numToEmit += addRemainder();
    //This will release per time step 
    break;
  }
}

int ParticleEmitter::addRemainder(){
  int rem = 0;

  remAmount += (emitAmount - floor(emitAmount));
  if(remAmount >= 1.0){
    rem = (int)floor(remAmount);
    remAmount -= floor(remAmount);
  }
  
  return rem;

}

//Determines whether or not it is time to emit a particle
bool ParticleEmitter::timeToEmit(float time_step){

  emitTime += time_step*releaseRate;

  if(emitTime >= 1.0){
    //Sets number of particles to emit
    numToEmit = (int)floor(emitTime);
    //Keeps track of the remainder value for emitTime.
    remTime += emitTime - floor(emitTime);

    //When the remainder gets greater than 1, release an extra particle.
    //This is an attempt to be more accurate in releasing particles,
    //instead of just "throwing away" the roundoff values
    if(remTime >= 1.0){
      numToEmit += (int)floor(remTime);
      remTime -= floor(remTime);
    }
    emitTime = 0;
    return true;
  }
  else return false;

}
