#include "pointEmitter.h"

PointEmitter::PointEmitter(float x,float y,float z,float rate, int* w, 
			   int* h,std::list<int>* ind, 
		       GLSLObject* emit_shader){

  xpos = x;
  ypos = y;
  zpos = z;

  reuse = false;
  lifeTime = -1.0;

  releaseRate = rate;
  numToEmit = 0;

  emitTime = 0;
  remTime = 0;

  twidth = *w;
  theight = *h;

  indices = ind; 

  shader = emit_shader;

}
PointEmitter::~PointEmitter(){}

int PointEmitter::EmitParticle(FramebufferObject* fbo, bool odd){
 
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
	  indices->pop_back();
	  if(reuse){
	    pIndex newIndex;
	    newIndex.id = p_index;
	    newIndex.time = 0;
	    indicesInUse->push_back(newIndex);	   
	  }

	  shader->activate();	

	  //Determine the coordinates into the position texture
	  int s = (p_index%twidth);
	  int t = (p_index/twidth);
	  //s = (s*(2.0/(float)twidth) - 1.0);
	  //t = (t*(2.0/(float)theight) - 1.0);
	  
	  glViewport(s,t,1,1);
       
	  glBegin(GL_POINTS);
	  {
	    glColor4f(xpos, ypos, zpos, lifeTime);
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
    //temp = numToEmit;
    //numToEmit = 0;
   
    return numToEmit;
}
