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
  
  // Shadow map setup is a simply a setup function for
  // the texture and frame buffer object that is used.
  void shadowMapSetup();

  // Generate shadow map will generate a shadow map
  // based on the position of the sun.
  // 
  // Currently, this function will only take snapshot
  // of what drawFeatures will draw. If additional 
  // geometry is added to the scene and not pushed
  // into that function, then this function must be
  // updated.
  void generateShadowMap();

  // genGridShadow calculates for each cell the percentage
  // viewable by the sun (or in other words the percentage
  // not in shadow). For this method to work, the shadow
  // map from the Sun's perspective must have been 
  // generated.
  void genGridShadow();
  
  // The following are angles used to calculate the Sun's position.
  // TODO: At some point these should be pushed into the Util class
  // and read in with the Settings file accociated with the .proj
  // file (or even the proj file it self).
  float sun_azimuth;
  float sun_altitude;

  // This is simply a static constant as to how far away
  // the sun is. At some point this value should either be
  // correctly set or dynamically calculated based on the size
  // of the scene.
  static const float SUN_DISTANCE = 1000;

 protected:
 
  virtual void initFBO();
  virtual ~MultipleBuildingsModel();

  Timer *mbaTimer;
  
  // shadowMap is the handle for the texture where the shadow 
  // (depth) map is stored on the graphics card.
  GLuint shadowMap;
  
  // shadowFBO is the frame buffer used to capture the shadow
  // map.
  FramebufferObject * shadowFBO;

};

#endif 
