#ifndef __ReflectionModel_H__
#define __ReflectionModel_H__

#include "plumeControl.h"

class ReflectionModel : public PlumeControl{

 public:
  
  ReflectionModel(Util*);
  virtual void init(bool);
  virtual int display();
  virtual void setupTextures();
  virtual void setupEmitters();
  
 protected:

  virtual ~ReflectionModel();



};

#endif //__ReflectionModel_H__
