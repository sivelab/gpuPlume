#ifndef __GaussianModel_H__
#define __GaussianModel_H__

#include "plumeControl.h"

class GaussianModel : public PlumeControl{

 public:
  
  GaussianModel(Util*);
  virtual void init(bool);
  virtual int display();
  virtual void setupTextures();
  //virtual void setupEmitters();
  
 protected:

  virtual ~GaussianModel();



};

#endif //__GaussianModel_H__
