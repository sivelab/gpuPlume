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
#include "particleEmitter.h"
#include "VisualPlane.h"
#include "framebufferObject.h"

enum tau_visual_type{draw_contours,draw_layers};

class DisplayControl{

 public:
  
  DisplayControl(int, int, int, GLenum,float,float,float);

  void drawVisuals(GLuint, GLuint, GLuint, int, int, int, GLuint postexid, GLuint veltexid);
  void drawAxes();
  void drawGrid();
  void drawGround();
  void drawSky();
  void drawLayers(GLuint, int, float);
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
  void initVars(int,int*,float*,float*,float*,
		      float*,float*,float*,float*);
  
  void setEmitter(ParticleEmitter*);
  void setVisualPlane(VisualPlane*);

  bool rotate_around, change_height, change_look;
  bool frame_rate;

  bool draw_buildings;
  bool osgPlume;

  float max_vel;
  // Determine which visual state to draw the particles.  The default
  // should probably be point.
  enum ParticleVisualState
    {
      PARTICLE_SPHERE,
      PARTICLE_SNOW,
    };
  ParticleVisualState particle_visual_state;

  int visual_layer;
  float eye_pos[3];
  float eye_gaze[3];
  
  tau_visual_type tau_visual;

  // sorting preference
  bool perform_cpu_sort;

  // TODO: DELETE ME, THIS IS FOR TESTING
  GLfloat inShadowData[30][30][30][4];
  void drawInShadowData();
  float sun_pos[3];
  
 private:

  void createImageTex(GLuint, char*);
  GLubyte* readPPM(char*, int*, int*);

  int nx;
  int ny;
  int nz;

  int nxdx;
  int nydy;
  int nzdz;

  float norm_x,norm_y,norm_z;
  float yaw,pitch,roll;

  void calculateNormal();

  //Feature Variables
  int numBuild;
  int* numSides;
  float* gamma;
  float* xfo;
  float* yfo;
  float* zfo;
  float* ht;
  float* wti;
  float* lti;
  
  
  //Copy of Particle Emitter for hand control
  ParticleEmitter* pe;

  //Copy of Visual Plane for hand control
  VisualPlane* plane;

  GLSLObject sphereParticle_shader;
  GLSLObject windField_shader;

  Timer *dTimer;

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
	GLuint skyBoxTex[6];

  GLenum texType;

  GLint uniform_vel_color; //, uniform_tauTex;
 
  //GLint uniform_max_x,uniform_max_y,uniform_max_z,uniform_max_c;
  //GLint uniform_min_x,uniform_min_y,uniform_min_z,uniform_min_c;
  //GLint uniform_controlWind, uniform_sliderWind, uniform_windTex;
  GLint uniform_windTex;
  GLint uniform_max_velocity;

  //
  // Point Sprite visual data
  //
  GLuint point_sprite_textures[2];

  // * variable to reference the texture units for the point sprite
  // * and normal map textures
  GLint uniform_pointsprite_tex, uniform_normalmap_tex, uniform_visualization_tex, uniform_pointsprite_visuals, 
    uniform_nx, uniform_ny, uniform_nz, uniform_numInRow;

  GLfloat uniform_max_vel;

  // * function to create point sprite textures
  void createPointSpriteTextures();

	void DrawSkyBox(float x, float y, float z,  float width, float height, float length);

};
