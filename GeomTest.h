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
  virtual void initFBO();
  virtual void setupTextures();
 
  
 protected:
  
  virtual ~GeomTest();

  void readInTables();

  GLint uniform_postPP, uniform_posPP, uniform_x,uniform_y;
  GLint u_slice, uniform_tau, u_tau3D, u_case, u_edge;

  GLSLObject geomShader, testShader, pathLineShader;

};

#endif 
