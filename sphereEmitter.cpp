#include <iostream>
#include "sphereEmitter.h"


SphereEmitter::SphereEmitter(float x,float y,float z,float rate, float r,
	       int* w,int* h,std::list<int>* ind, 
	       GLSLObject* emit_shader){

  xpos = x;
  ypos = y;
  zpos = z;

  reuse = false;
  lifeTime = -1.0;

  releaseRate = rate;
  radius = r;

  numToEmit = 0;

  emitTime = 0;
  remTime = 0;

  twidth = *w;
  theight = *h;

  indices = ind;
 
  shader = emit_shader;

}
SphereEmitter::~SphereEmitter(){}

/*void SphereEmitter::getPosition(float*x, float*y, float*z){
  *x = xpos;
  *y = ypos;
  *z = zpos;
  //r = radius;

  }*/
void SphereEmitter::getReleasedPosition(float*x,float*y,float*z){
  *x = xpos + offsetx;
  *y = ypos + offsety;
  *z = zpos + offsetz;

}

void SphereEmitter::Draw(){
  glDisable(GL_TEXTURE_RECTANGLE_ARB);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  //glPointSize(2.0);
  
  glPushMatrix();
  glColor4f(0.0, 0.0, 1.0, 0.5);
  glTranslatef(xpos,ypos,zpos);
  glutSolidSphere(radius, 20, 16);

  glPopMatrix();

  glDisable(GL_BLEND);
  glDisable(GL_COLOR_MATERIAL);
  glDisable(GL_LIGHT0);
  glDisable(GL_LIGHTING);
  glEnable(GL_TEXTURE_RECTANGLE_ARB);
}

int SphereEmitter::EmitParticle(FramebufferObject* fbo, bool odd){
 
  int p_index;
  
  //Make sure there are available indices to emit particles.
    if(!indices->empty()){
      
#if USE_VERTEX_PARTICLE_WRITE
      if(odd)
	glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
      else 
	glDrawBuffer(GL_COLOR_ATTACHMENT1_EXT);
	        
      shader->activate();	
#endif

      //Do this for each particle that is being emitted.
      for(int i = 0; i < numToEmit; i++){

	if(!indices->empty()){

	  offsetx = Random::uniform()*2.0 - 1.0;
	  offsety = Random::uniform()*2.0 - 1.0;
	  offsetz = Random::uniform()*2.0 - 1.0;

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
	    newIndex.time = 0.0;
	    indicesInUse->push_back(newIndex);
	  }
	  
	  //Determine the coordinates into the position texture
	  s = (p_index%twidth);
	  t = (p_index/twidth);
	  
#if 0
	  std::cout << "\tp_index = " << p_index << ", (s,t) = (" << s << ", " << t << "), position = (" 
		    << xpos + offsetx << ", " 
		    << ypos + offsety << ", " 
		    << zpos + offsetz << ", "
		    << lifeTime << ")" << std::endl;
#endif

#if USE_VERTEX_PARTICLE_WRITE
	  glViewport(s,t,1,1);
	  
	  glMatrixMode(GL_PROJECTION);
	  glLoadIdentity();
	  gluOrtho2D(0-0.5, twidth+0.5, 0-0.5, theight+0.5);
	  glMatrixMode(GL_MODELVIEW);
	  glLoadIdentity();

	  // glPointSize(0.10);
	  // glPointSize(1.0);
	  glBegin(GL_POINTS);
	  {
	    glColor4f(xpos + offsetx, ypos + offsety, zpos + offsetz, lifeTime);
	    glVertex2f(s, t);
	  }
	  glEnd();
	  glFinish();
#endif

	  // Second mechanism to release particles.  Uses a texture
	  // copy which may likely be more expensive even though we're
	  // doing 1x1 pixels (but many times).  This operation
	  // appears to work consistently.
	  GLfloat value[4];
	  value[0] = xpos + offsetx;
	  value[1] = ypos + offsety;
	  value[2] = zpos + offsetz;
	  value[3] = lifeTime;

	  // write there via a glTexSubImage2D
	  if(odd)
	    // will read from this texture next
	    glBindTexture(GL_TEXTURE_RECTANGLE_ARB, m_posTexID0);
	  else 
	    // will read from this texture next
	    glBindTexture(GL_TEXTURE_RECTANGLE_ARB, m_posTexID1);

	  glTexSubImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, s, t, 1, 1, GL_RGBA, GL_FLOAT, value);
	}	
	else 
	  {
	    std::cout << "No more indices!!!" << std::endl;
	  }
      }

#if USE_VERTEX_PARTICLE_WRITE
	  shader->deactivate();

      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      gluOrtho2D(-1, 1, -1, 1);
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
#endif

      return numToEmit;
    }
    else 
      {
    //temp = numToEmit;
    //numToEmit = 0;
	return 0;
      }

}
