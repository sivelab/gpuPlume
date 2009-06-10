#ifndef __Gaussian_2shaders_Model_H__
#define __Gaussian_2shaders_Model_H__

#include "plumeControl.h"

class Gaussian_2shaders_Model : public PlumeControl{

 public:
  
  Gaussian_2shaders_Model(Util*);
  virtual void init(bool);
  virtual int display();
  virtual void setupTextures();
  //virtual void setupEmitters();
  
 protected:

  virtual ~Gaussian_2shaders_Model();



};

#endif //__Gaussian_2shaders_Model_H__
