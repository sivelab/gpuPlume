#include <iostream>
#include "sphereEmitter.h"

#ifdef WIN32
// Rand functions
float rVal() { return (float)(rand()/(float)RAND_MAX); } 
#else
float rVal() { return drand48(); }
#endif

SphereEmitter::SphereEmitter(float x,float y,float z,float rate, float r,
	       int* w,int* h,std::list<int>* ind, 
	       GLSLObject* emit_shader){

  xpos = x;
  ypos = y;
  zpos = z;

  reuse = false;
  lifeTime = -1.0;

  pps = rate;
  radius = r;

  numToEmit = 1;
  emitTime = 1;
  remTime = 0;

  twidth = *w;
  theight = *h;

  indices = ind;
 
  shader = emit_shader;

}
SphereEmitter::~SphereEmitter(){}

void SphereEmitter::Draw(){
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  //glPointSize(2.0);
  glColor4f(0.0, 0.0, 1.0, 0.5);
  glPushMatrix();

  glTranslatef(xpos,ypos,zpos);
  glutSolidSphere(radius, 20, 16);

  glPopMatrix();

  glDisable(GL_BLEND);
  glDisable(GL_COLOR_MATERIAL);
  glDisable(GL_LIGHT0);
  glDisable(GL_LIGHTING);
}

int SphereEmitter::EmitParticle(FramebufferObject* fbo, bool odd){
 
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

	  offsetx = rVal()*2.0 - 1.0;
	  offsety = rVal()*2.0 - 1.0;
	  offsetz = rVal()*2.0 - 1.0;
	  float d = sqrt(offsetx*offsetx + offsety*offsety + offsetz*offsetz);
	  offsetx = (offsetx/d) * radius;
	  offsety = (offsety/d) * radius;
	  offsetz = (offsetz/d) * radius;
	
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
	  
	  glViewport(s,t,1,1);
       
	  glBegin(GL_POINTS);
	  {
	    glColor4f(xpos + offsetx, ypos + offsety, zpos + offsetz, lifeTime);
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
    temp = numToEmit;
    numToEmit = 1;
    return temp;
}
