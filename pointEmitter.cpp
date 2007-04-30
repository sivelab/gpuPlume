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

#if USE_VERTEX_PARTICLE_WRITE
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      gluOrtho2D(0, twidth, 0, theight);
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
  
      if(odd)
	glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
      else 
	glDrawBuffer(GL_COLOR_ATTACHMENT1_EXT);
#endif
      
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
	  s = (p_index%twidth);
	  t = (p_index/twidth);
	  //s = (s*(2.0/(float)twidth) - 1.0);
	  //t = (t*(2.0/(float)theight) - 1.0);
	  
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
	    glBindTexture(GL_TEXTURE_RECTANGLE_ARB, m_posTexID0);
	  else 
	    // will read from this texture next
	    glBindTexture(GL_TEXTURE_RECTANGLE_ARB, m_posTexID1);

	  glTexSubImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, s, t, 1, 1, GL_RGBA, GL_FLOAT, value);

#if USE_VERTEX_PARTICLE_WRITE
	  glViewport(s,t,1,1);
       
	  glBegin(GL_POINTS);
	  {
	    glColor4f(xpos, ypos, zpos, lifeTime);
	    glVertex2f(s, t);
	  }
	  glEnd();

	  shader->deactivate();
#endif
	}
	else 
	  {
	    std::cout << "No more indices!!!" << std::endl;
	  }
      }

#if USE_VERTEX_PARTICLE_WRITE
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      gluOrtho2D(-1, 1, -1, 1);
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
#endif

    }
    //temp = numToEmit;
    //numToEmit = 0;
   
    return numToEmit;
}
