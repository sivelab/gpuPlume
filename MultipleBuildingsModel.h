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

  //FramebufferObject* fbo2;

  virtual void initFBO();
  virtual ~MultipleBuildingsModel();



};

#endif //__MBModel_H__
