///////////////////////////////////////////////////////////
// This class takes care of drawing the axes for the domain,
// displaying the layers of the wind field, and drawing the
// frame rate to the screen.
///////////////////////////////////////////////////////////
#include <GL/glew.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "GLSL.h"
#include "Timer.h"

class DisplayControl{

 public:
  
  DisplayControl(int, int, int, GLenum);
  
  void drawVisuals(GLuint, GLuint, GLuint, int, int, int);
  void drawAxes();
  void drawGrid();
  void drawGround();
  void drawSky();
  void drawLayers(GLuint, int);
  void drawTurbulenceLayers(GLuint, int);
  void drawFeatures(void);
  void drawFrameRate(int, int);
  void OpenGLText(int, int, char*);
  void increaseVisualLayer();
  void decreaseVisualLayer();
  void moveForwardorBack(float);
  void slideLeftorRight(float);
  void setAzimuth(float, float);
  void setElevation(float, float);
  void setRotateAround(float);
  void lookUporDown(float);
  void initVars(int,float*,float*,float*,
		      float*,float*,float*);

  bool rotate_around, change_height, change_look;
  bool frame_rate;

  bool draw_buildings;
  bool osgPlume;

  // Determine which visual state to draw the particles.  The default
  // should probably be point.
  enum ParticleVisualState
    {
      PARTICLE_POINT,
      PARTICLE_SPRITE
    };
  ParticleVisualState particle_visual_state;

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
  float* xfo;
  float* yfo;
  float* zfo;
  float* ht;
  float* wti;
  float* lti;
  

  GLSLObject render_shader, turbulence_shader;

  Timer *clock_timer;
  Timer_t graphics_time[2];
  Timer_t HUP_display_update_time[2];

  double estimated_rate, last_estimated_rate;

  GLfloat azimuth;
  GLfloat elevation;
  
  double tranx,trany,tranz;
  double angle,yangle;
  double xlook,ylook,zlook;
  double xslide,yslide;

  GLuint displayTex[3];

  GLenum texType;

  GLint uniform_vel_color, uniform_tauTex;

  //
  // Point Sprite visual data
  //
  GLuint point_sprite_textures[2];

  // * variable to reference the texture units for the point sprite
  // * and normal map textures
  GLint uniform_pointsprite_tex, uniform_normalmap_tex, uniform_visualization_tex, uniform_pointsprite_visuals, 
    uniform_nx, uniform_ny, uniform_nz, uniform_numInRow;

  // * function to create point sprite textures
  void createPointSpriteTextures();
};
