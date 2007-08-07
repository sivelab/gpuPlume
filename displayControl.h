///////////////////////////////////////////////////////////
// This class takes care of drawing the axes for the domain,
// displaying the layers of the wind field, and drawing the
// frame rate to the screen.
///////////////////////////////////////////////////////////
#include <GL/glew.h>
#include <GL/glut.h>
#include "GLSL.h"
#include "Timer.h"

class DisplayControl{

 public:
  
  DisplayControl(int, int, int, GLenum);
  
  void drawVisuals(GLuint, GLuint, int, int, int);
  void drawAxes();
  void drawGround();
  void drawLayers(GLuint, int);
  void drawFeatures(void);
  void drawFrameRate(int, int);
  void OpenGLText(int, int, char*);
  void increaseVisualLayer();
  void decreaseVisualLayer();
  void moveForwardorBack(float);
  void setAzimuth(float, float);
  void setElevation(float, float);
  void setRotateAround(float, float,int,int);
  void initVars(int,double*,double*,double*,
		      double*,double*,double*);

  bool rotate_around, change_height, translate_view;
  bool frame_rate;

  bool draw_buildings;
  bool osgPlume;

  int visual_layer;
  float eye_pos[3];
  float eye_gaze[3];

 private:
  
  void createImageTex(GLuint, char*);
  GLubyte* readPPM(char*, int*, int*);
  
  int nx;
  int ny;
  int nz;

  //Feature Variables
  int numBuild;
  double* xfo;
  double* yfo;
  double* zfo;
  double* ht;
  double* wti;
  double* lti;
  

  GLSLObject render_shader;

  Timer *clock_timer;
  Timer_t graphics_time[2];

  GLfloat azimuth;
  GLfloat elevation;
  
  double tranx,trany,tranz;
  double angle;
  double xlook,ylook,zlook;

  GLuint displayTex[3];

  GLenum texType;

};
