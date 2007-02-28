
#include "particleEmitter.h"
#include <math.h>

ParticleEmitter::ParticleEmitter(float x,float y,float z,float rate, int* w, 
		       int* h,std::list<int>* ind, GLSLObject* emit_shader){

  xpos = x;
  ypos = y;
  zpos = z;

  pps = rate;
  numToEmit = 1;

  twidth = *w;
  theight = *h;

  indices = ind;
  shader = emit_shader;
  
}
bool ParticleEmitter::timeToEmit(float time_step){
  emitTime += time_step*pps;

  if(emitTime >= 1.0){
    remTime += emitTime - floor(emitTime);

    if(remTime >= 1.0){
      numToEmit += (int)floor(remTime);
      remTime -= floor(remTime);
    }
    emitTime = 0;
    return true;
  }
  else return false;
}
void ParticleEmitter::EmitParticle(FramebufferObject* fbo, bool odd){
 
  int p_index;
  //Make sure there are available indices to emit particles.
    if(!indices->empty()){

      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      gluOrtho2D(0, twidth, 0, theight);
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
  
      if(odd)
	glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
      else 
	glDrawBuffer(GL_COLOR_ATTACHMENT1_EXT);

      
      //Do this for each particle that is being emitted.
      for(int i = 0; i < numToEmit; i++){
	if(!indices->empty()){
	
	  //First get available index
	  p_index = indices->back();
	  //std::cout << p_index <<std::endl;
	  indices->pop_back();

	  shader->activate();	

	  //Determine the coordinates into the position texture
	  int s = (p_index%twidth);
	  int t = (p_index/twidth);
	  //s = (s*(2.0/(float)twidth) - 1.0);
	  //t = (t*(2.0/(float)theight) - 1.0);
	  //std::cout << s << " " << t  <<std::endl;
	  glViewport(s,t,1,1);
       
	  glBegin(GL_POINTS);
	  {
	    glColor4f(xpos, ypos, zpos, 1.0);
	    glVertex2f(s, t);
	  }
	  glEnd();

	  shader->deactivate();
	}
      }
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      gluOrtho2D(-1, 1, -1, 1);
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();

    }
    numToEmit = 1;
}
