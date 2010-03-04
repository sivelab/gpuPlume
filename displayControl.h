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
#include "util.h"

/**
 *  The following is for synchronizing gpuPlume across
 *  multiple computers over a network and communicating
 *  with the TPAWT. Since, this code depends on Unix
 *  based network libraries it is being excluded from the
 *  Windows build.
 */
#if !WIN32
#include "NetworkManager.h"
#include "TreadportManager.h"
#include "graphicsUtil.h"
#include <vector>

struct screen_t {
  float asp_ratio;
  std::vector<float> lr, ll, ul, ur; // lowerRight, lowerLeft, upperLeft, upperRight == screen edges.
  std::vector<float> normal, center;
  std::vector<float> far1, far2, far3, far4;
  std::vector<float> proj;
  float phi;
};
#endif

enum tau_visual_type{draw_contours,draw_layers};

class DisplayControl{

 public:
  
  DisplayControl(int, int, int, GLenum, bool, Util *);
  ~DisplayControl();

  void initTreadport();
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
  float eye_offset[3];

  tau_visual_type tau_visual;

  // sorting preference
  bool perform_cpu_sort;
  
  /**
   * ViewingMode is the enum which describes the viewing state that Display
   * Control runs in. Note that some extra behavior is turned on when
   * different modes are enabled. For example, when TREADPORT is enabled
   * then information is communicated with the treadport system, otherwise
   * no communication occurs.
   */
  enum ViewingMode {
    STANDARD = 0,
    VR = 1,
    TREADPORT = 2,
    ORTHOGRAPHIC_TOP = 3
  };

  /**
   * initializeView creates the correct view based on the viewing mode that has
   * been specified (by viewingMode, and thus the in the config file). Note that
   * if we are using OpenSceneGraph this method does nothing.
   */
  void initializeView();

  /**
   * This method unitializes the view, which is needed because their a two
   * "views" for the treadport.
   */
  void deinitializeView();

  /**
   * syncTreadportData will syncs the data with the treadport system. Note that
   * this function will preform it's action only when the application is the
   * set to broadcast and the treadport view mode is enabled. IMPORTANT: This
   * sync should be done BEFORE the call to syncDataOverNetwork so that the
   * latest data reaches all of the computers.
   */
  void syncTreadportData();

  /**
   * syncDataOverNetwork will use the NetworkManager to sync the data contained
   * in the NetworkSyncData struct over the network (see the details of the
   * function as to what data is synced).
   */
  void syncDataOverNetwork();

  /**
   * getInPauseMode will retrieve the inPauseMode that has been stored within
   * DisplayControl. The value is stored here so that it can be syncronized
   * over the network to the other screens/eyes. Note that this value should
   * always be syncrhonized with the inPauseMode within which ever of
   * PlumeControl's sublcasses is being used (i.e. MultipleBuildingsModel).
   */
  bool getInPauseMode();
  
  /**
   * setInPauseMode will set the inPauseMode that has been stored within
   * DisplayControl. The value is set here so that it can be syncrhonized
   * over the network to the other screens/eyes. Note that this value should
   * always be syncrhonized with the inPauseMode within which ever of
   * PlumeControl's sublcasses is being used (i.e. MultipleBuildingsModel).
   *
   * @param newInPauseMode is a bool representing the a new inInPauseMode
   * variable to set set in DisplayControl and synchronized to the other
   * computers (if the correct network mode has been set).
   *
   */
  void setInPauseMode(bool newInPauseMode);
  
  // The following is for displaying the
  // cell shadow data.
  GLfloat * inShadowData;
  GLfloat * inShadowData2;
  void drawInShadowData();
  float sun_pos[3];
  bool drawISD;

  // The following is for visualizing the wind direction.
  GLfloat windDir[3];
  
 private:
  
#if !WIN32
  /**
   * network is the NetworkManager object that is used to synchronize data over
   * the local network to the other screens (or rather the other computers that
   * generate the image for the other screens).
   */

  NetworkManager network;

  /**
   * treadport is the TreadportManager object that is responsible for communicating
   * with the treadport system.
   */
  TreadportManager * treadport;
#endif

  /**
   * static_treadport_frustum is a flag that determines of a static or dynamic view
   * frustum should be used when we are in TREADPORT mode. This value should be set
   * within the Settings/input.txt file (or config file that has been passed in to
   * gpuPlume).
   */
  bool static_treadport_frustum;

  /**
   * inPauseMode is a duplicate of the inPauseMode flag which is set within
   * PlumeControl. This is done so that the value can be synchronized over 
   * the network when we are running more than one screen.
   */
  bool inPauseMode;
  
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
  
  /**
   * util is stored so that DisplayControl can access settings that are
   * retrieved from the config file that is read-in upon program start.
   */
  Util * util;

  /**
   * viewingMode is an enum which describes which viewing mode Display Control
   * will run in. Note that the default should be STANDARD or what has been
   * set within the input.txt config file.
   */
  ViewingMode viewingMode;

#if !WIN32
  /**
   * This is the mode in which the local communication will operate. The choices
   * are set in the config file (Settings/input.txt or passed in config file).
   */
  NetworkManager::Mode network_mode;

  // For treadport configuration.
  screen_t screen;
  std::vector<float> eye;
#endif
  
  /**
   * calcTreadportFrustum calculates and sets the frustum and opengl lookat in
   * order to set the view correctly for the treadport configuation.
   *
   * @param view is a charecter representing which screen the view should be
   * calculated for (either, l = left, c = center, or r = right).
   */
  void calcTreadportFrustum(char view);
  
  /**
   * blankSlides will black out the correct portion of the screen
   * when the treadport is in use.
   */
  void blankSides();
  
};
