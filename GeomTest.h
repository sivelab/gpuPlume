#ifndef __GM_H__
#define __GM_H__

#include <list>
#include <vector>
#include "plumeControl.h"

class GeomTest : public PlumeControl{

 public:
  
  GeomTest(Util*);
  virtual void init(bool);
  virtual int display();
 
  
 protected:
  
  virtual ~GeomTest();

  GLint uniform_postPP, uniform_posPP, uniform_x,uniform_y;

  GLSLObject geomShader,testAdvectShader, pathLineShader;

};

#endif 
