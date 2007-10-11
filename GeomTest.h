#ifndef __GM_H__
#define __GM_H__

#include <list>
#include <vector>
#include "plumeControl.h"

typedef struct{
  //Index into position texture
  int x;
  int y;
  //Index into pathLine texture
  int s;
  int t;
}pathIdx;

class GeomTest : public PlumeControl{

 public:
  
  GeomTest(Util*);
  virtual void init(bool);
  virtual int display();
  virtual void setupTextures();
  //virtual void setupEmitters();
  
 protected:
  int pwidth,pheight;

  std::list<pathIdx> pathList;

  std::list<pathIdx>::iterator pathIter;
  
  virtual void initFBO();
  virtual ~GeomTest();

  GLint uniform_postPP, uniform_posPP, uniform_x,uniform_y;

  GLenum path_buffer,paths;

  GLSLObject geomShader,testAdvectShader, pathLineShader;

};

#endif 
