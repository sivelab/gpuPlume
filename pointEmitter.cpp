#include "pointEmitter.h"

PointEmitter::PointEmitter(float x,float y,float z,float rate, int w, 
			   int h,std::list<int>* ind,GLSLObject* emit_shader){

  xpos = x;
  ypos = y;
  zpos = z;

  reuse = false;
  lifeTime = -1.0;

  releaseRate = rate;
  numToEmit = 0;

  emitTime = 0;
  remTime = 0;

  twidth = w;
  theight = h;

  indices = ind; 

  shader = emit_shader;

}
PointEmitter::~PointEmitter(){}

void PointEmitter::setVertices(){

  posCoord.clear();

  for(int i=0; i<numToEmit; i++){
    posCoord.push_back(xpos);
    posCoord.push_back(ypos);
    posCoord.push_back(zpos);
    posCoord.push_back(lifeTime);
    
  }
}
