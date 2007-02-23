
#include "particleEmitter.h"

ParticleEmitter::ParticleEmitter(float x,float y,float z, int* w, int* h, 
				 std::list<int>* ind){

  xpos = x;
  ypos = y;
  zpos = z;

  numToEmit = 1;

  twidth = *w;
  theight = *h;

  indices = ind;
  
}
void ParticleEmitter::setNumToEmit(int num){
  numToEmit = num;
  
}
void ParticleEmitter::EmitParticle(FramebufferObject* fbo, GLSLObject emit_shader,
				   bool odd){
 
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

	  emit_shader.activate();	

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

	  emit_shader.deactivate();
	}
      }
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      gluOrtho2D(-1, 1, -1, 1);
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();

    }
}
