#ifndef __GM_H__
#define __GM_H__

#include <list>
#include <vector>
#include "plumeControl.h"
#include "IsoSurface.h"

class GeomTest : public PlumeControl{

 public:
  
  GeomTest(Util*);
  virtual void init(bool);
  virtual int display();
  virtual void initFBO();
  virtual void setupTextures();
 
    
 protected:
  
  virtual ~GeomTest();

  //void readInTables();

  //IsoSurface* isoSurface;

  GLint uniform_postPP, uniform_posPP, uniform_x,uniform_y;
  //GLint u_slice, uniform_tau, u_tau3D, u_case, u_edge;
  //GLint u_dx,u_dy,u_dz,u_mesh;

  //GLSLObject geomShader, testShader, pathLineShader;
  //GLSLObject isoShader;

 

};

#endif 
