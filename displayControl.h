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
  
  void setupTurbulenceShader(float [], float [], float*, float*);
  void setupWindFieldShader(float [], float []);

  void drawVisuals(GLuint, GLuint, GLuint, int, int, int);
  void drawAxes();
  void drawGrid();
  void drawGround();
  void drawSky();
  void drawScale();
  void drawLayers(GLuint,GLuint, int);
  void drawTurbulenceLayers(GLuint, int);
  void drawFeatures(void);
  void drawFrameRate(int, int);
  void OpenGLText(int, int, char*);
  void increaseVisualLayer();
  void decreaseVisualLayer();
  void moveSliderUp();
  void moveSliderDown();
  void moveSlider(int);
  bool clickedSlider(int,int);
  void moveForwardorBack(float);
  void slideLeftorRight(float);
  void setAzimuth(float, float);
  void setElevation(float, float);
  void setRotateAround(float);
  void lookUporDown(float);
  void initVars(int,float*,float*,float*,
		      float*,float*,float*);

  bool rotate_around, change_height, change_look, move_slider;
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
  int visual_field;
  float eye_pos[3];
  float eye_gaze[3];
  //int tauValue; 

  int slider_x;
  bool localValues;

 private:
  
  void createImageTex(GLuint, char*);
  GLubyte* readPPM(char*, int*, int*);
  
  //Global Max and Min values of Tau
  float* TauMax;
  float* TauMin;
  float* tauLocalMax;
  float* tauLocalMin;
  float tauMin, tauMax;
  char* Taus[4];
  
  //Global Max and Min values of Wind Field
  float* WindMax;
  float* WindMin;
  float windMin, windMax;
  char* windTitle[4];

  float scale_xstart,scale_xend;
  float scale_ystart,scale_yend;

  //Controls what value the middle color is
  //represented at in the scale used for visualizing
  //the turbulence layers.
  float slider;

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
  

  GLSLObject render_shader, turbulence_shader, scale_shader;
  GLSLObject windField_shader;

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
  GLint uniform_max11,uniform_max22,uniform_max33,uniform_max13;
  GLint uniform_min11,uniform_min22,uniform_min33,uniform_min13;
  GLint uniform_controlTau;
  GLint uniform_xmax, uniform_xmin, uniform_tauMin, uniform_tauMax;
  GLint uniform_sliderScale, uniform_sliderTurb;
  GLint uniform_max_x,uniform_max_y,uniform_max_z,uniform_max_c;
  GLint uniform_min_x,uniform_min_y,uniform_min_z,uniform_min_c;
  GLint uniform_controlWind, uniform_sliderWind, uniform_windTex;


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
