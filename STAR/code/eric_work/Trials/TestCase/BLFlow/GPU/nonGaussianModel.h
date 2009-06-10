#ifndef __nonGaussianModel_H__
#define __nonGaussianModel_H__

#include "plumeControl.h"

class NonGaussianModel : public PlumeControl{

 public:
  
  NonGaussianModel(Util*);
  virtual void init(bool);
  virtual int display();
  virtual void setupTextures();
  //virtual void setupEmitters();
  
 protected:

  virtual void initFBO();
  virtual ~NonGaussianModel();



};

#endif //__nonGaussianModel_H__
