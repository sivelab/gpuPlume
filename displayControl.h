///////////////////////////////////////////////////////////
// This class takes care of drawing the axes for the domain,
// displaying the layers of the wind field, and drawing the
// frame rate to the screen.
///////////////////////////////////////////////////////////

#include <iostream>
#include <GL/glew.h>

#include "Timer.h"

class DisplayControl{

 public:
  
  DisplayControl(int, int, int, GLenum);
  
  void drawAxes();
  void drawLayers(int, GLuint, int);
  void drawFeatures(void);
  void drawFrameRate(int, int);
  void OpenGLText(int, int, char*);

 private:
  int nx;
  int ny;
  int nz;

  Timer *clock_timer;
  Timer_t graphics_time[2];

  GLenum texType;

};
