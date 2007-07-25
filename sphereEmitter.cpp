#include <iostream>
#include "sphereEmitter.h"
#include "Random.h"


SphereEmitter::SphereEmitter(float x,float y,float z,float rate, float r,
	       int w,int h,std::list<int>* ind, GLSLObject* emit_shader){

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

  twidth = w;
  theight = h;

  indices = ind;
 
  shader = emit_shader;

}
SphereEmitter::~SphereEmitter(){}

void SphereEmitter::getReleasedPosition(float*x,float*y,float*z){
  *x = xpos + offsetx;
  *y = ypos + offsety;
  *z = zpos + offsetz;
}

void SphereEmitter::Draw(){
  glDisable(GL_TEXTURE_RECTANGLE_ARB);
  //glEnable(GL_LIGHTING);
  //glEnable(GL_LIGHT0);
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
  //glDisable(GL_LIGHT0);
  //glDisable(GL_LIGHTING);
  glEnable(GL_TEXTURE_RECTANGLE_ARB);
}

void SphereEmitter::setVertices(){
  posCoord.clear();

  for(int i=0; i < numToEmit; i++){
    offsetx = Random::uniform()*2.0 - 1.0;
    offsety = Random::uniform()*2.0 - 1.0;
    offsetz = Random::uniform()*2.0 - 1.0;

    float d = sqrt(offsetx*offsetx + offsety*offsety + offsetz*offsetz);
    offsetx = (offsetx/d) * radius;
    offsety = (offsety/d) * radius;
    offsetz = (offsetz/d) * radius;
    
    posCoord.push_back(xpos+offsetx);
    posCoord.push_back(ypos+offsety);
    posCoord.push_back(zpos+offsetz);
    posCoord.push_back(lifeTime);

  }
  //int size = posCoord.size();
  //std::cout << "Stored " << size << " positions" << std::endl;

}
int SphereEmitter::EmitParticle(bool odd,GLuint pos0, GLuint pos1,
			       float time_step){
  int p_index;

  switch(releaseType){
  case perSecond:
    if(!timeToEmit(time_step))
	numToEmit = 0;
    break;
  case onePerKeyPress:
    setNumToEmit(1);
    emit = false;
    break;
  case instantaneous:
    setNumToEmit(twidth*theight);
    emit = false;
    break;
  default:
    //This will release per time step with numToEmit already defined!!!
    break;
  }

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

	offsetx = Random::uniform()*2.0 - 1.0;
	offsety = Random::uniform()*2.0 - 1.0;
	offsetz = Random::uniform()*2.0 - 1.0;

	float d = sqrt(offsetx*offsetx + offsety*offsety + offsetz*offsetz);
	offsetx = (offsetx/d) * radius;
	offsety = (offsety/d) * radius;
	offsetz = (offsetz/d) * radius;


	//Determine the coordinates into the position texture
	s = (p_index%twidth);
	t = (p_index/twidth);	  

	if(Punch_Hole){
	  glPointSize(1.0);

	  shader->activate();

	  glViewport(s,t,1,1);
	  //std::cout << "particle num= " << p_index << "  s = " << s << "  t = " << t << std::endl;
	  glBegin(GL_POINTS);
	  {
	    glColor4f(xpos+offsetx, ypos+offsety, zpos+offsetz, lifeTime);
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
	  value[0] = xpos+offsetx;
	  value[1] = ypos+offsety;
	  value[2] = zpos+offsetz;
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
