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
