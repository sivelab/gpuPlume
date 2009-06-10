#ifndef __MBModel_H__
#define __MBModel_H__

#include "plumeControl.h"

class MultipleBuildingsModel : public PlumeControl{

 public:
  
  MultipleBuildingsModel(Util*);
  virtual void init(bool);
  virtual int display();
  virtual void setupTextures();
  //virtual void setupEmitters();
  
 protected:

 
  virtual void initFBO();
  virtual ~MultipleBuildingsModel();



};

#endif 
