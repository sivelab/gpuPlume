#include <iostream>
#include <GL/glew.h>


class DisplayControl{

 public:
  
  DisplayControl(int, int, int, GLenum);
  
  void drawAxes();
  void drawLayers(int, GLuint, int);
  void drawFrameRate(int, int);
  void OpenGLText(int, int, char*);

 private:
  int nx;
  int ny;
  int nz;

  struct timeval startframe, endframe, diff;

  GLenum texType;

};
