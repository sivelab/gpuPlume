///////////////////////////////////////////////////////////
// This class takes care of drawing the axes for the domain,
// displaying the layers of the wind field, and drawing the
// frame rate to the screen.
///////////////////////////////////////////////////////////

#include <iostream>
#include <GL/glew.h>
#include "GLSL.h"

#include "Timer.h"

class DisplayControl{

 public:
  
  DisplayControl(int, int, int, GLenum);
  
  void drawVisuals(GLuint, GLuint, int, int, int);
  void drawAxes();
  void drawLayers(GLuint, int);
  void drawFeatures(void);
  void drawFrameRate(int, int);
  void OpenGLText(int, int, char*);
  void increaseVisualLayer();
  void decreaseVisualLayer();
  void setEyeValues(float);
  void setAzimuth(float, float);
  void setElevation(float, float);

  bool rotate_sphere, rotate_object, translate_view;
  bool frame_rate;

  bool draw_buildings;

  int visual_layer;

 private:
  int nx;
  int ny;
  int nz;

  GLSLObject render_shader;

  Timer *clock_timer;
  Timer_t graphics_time[2];

  GLfloat azimuth;
  GLfloat elevation;
  
  
  float eye_pos[3];
  float eye_gaze[3];


  GLenum texType;

};
