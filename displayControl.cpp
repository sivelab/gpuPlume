#include <cstdlib>
#include <iostream>
#include <cstring>
#include <math.h>
#include <assert.h>
#include <new>
#include "displayControl.h"
#include "glErrorUtil.h"

static char text_buffer[128];
//static char number[128];
#ifdef WIN32
#define M_PI 3.141592654
#define M_PI_2 M_PI/2.0
#endif

// sorting tests /////////////////////////
float sort_eye[3];
double Distance2Eye(const float* fp)
{
  return sqrt( (sort_eye[0]-fp[0])*(sort_eye[0]-fp[0]) 
               + (sort_eye[1]-fp[1])*(sort_eye[1]-fp[1]) 
               + (sort_eye[2]-fp[2])*(sort_eye[2]-fp[2]) );
}

static int cmpVertices(const void *p1, const void *p2)
{
  const float* fp1 = (const float*)p1;
  const float* fp2 = (const float*)p2;

  float d1 = Distance2Eye(fp1);
  float d2 = Distance2Eye(fp2);

  if (d2 < d1)
    return -1;
  else 
    return 1;
}


DisplayControl::DisplayControl(int x, int y, int z, GLenum type, bool initialPauseMode, Util * newUtil)
{
  
  inPauseMode = initialPauseMode;
  
  util = newUtil;

  float dx = util->dx;
  float dy = util->dy;
  float dz = util->dz;
  
#if !WIN32
  network_mode = (NetworkManager::Mode)util->network_mode;
#endif

  viewingMode = (ViewingMode)util->viewing_mode;
  // viewingMode = ORTHOGRAPHIC_TOP;
  
  if(util->static_treadport_frustum == 1) {
    static_treadport_frustum = true;
  } else {
    static_treadport_frustum = false;
  }
  
  if(viewingMode == STANDARD) {
    std::cout << "Using STANDARD view (" << viewingMode << ")." << std::endl;
  } else if(viewingMode == VR) {
    std::cout << "Using VR view (" << viewingMode << ")." << std::endl;
  } else if(viewingMode == TREADPORT) {
    std::cout << "Using TREADPORT view (" << viewingMode << ")." << std::endl;
  }
  
  dTimer = new Timer(true);

  nx = x;
  ny = y;
  nz = z;

  nzdz = (int)(nz*(1.0/dz));
  nydy = (int)(ny*(1.0/dy));
  nxdx = (int)(nx*(1.0/dx));

  texType = type;

  // Adjust the default height, if we are using a system (such as TPAWT)
  // that adds in an actual measured height, then the default value
  // should be 0.
  if(viewingMode == TREADPORT) {
    eye_pos[0] = 10;
    eye_pos[1] = 0;
    eye_pos[2] = 0;
  } if (viewingMode == ORTHOGRAPHIC_TOP)
      {
	eye_pos[0] = 0;  
	eye_pos[1] = 0;
	eye_pos[2] = 20;
      }
  else {
    eye_pos[0] = nx+50;  
    eye_pos[1] = 0;
    eye_pos[2] = 5;
  }
  
  eye_gaze[0] = -1.0;
  eye_gaze[1] = 0;
  eye_gaze[2] = 0;

#if 0
  // starts the view closer to the particle release
  eye_pos[0] = 20;
  eye_pos[1] = -10;
  eye_pos[2] = 5;
  
  eye_gaze[0] = 0.0;
  eye_gaze[1] = 1.0;
  eye_gaze[2] = 0.0;
#endif

#if !WIN32
  // Set up the default eye (for the treadport system)
  eye.resize(3);
  eye[0] = 0;
  eye[1] = 0;
  eye[2] = 1.2;
#endif

  //New calculations for moving and looking around
  yaw = 0.0;
  pitch = 0.0;
  roll = 0.0;
  norm_x = -1.0;
  norm_y = 0.0;
  norm_z = 0.0;

  angle = M_PI;
  yangle = 0.0;

  xslide = 0.0;
  yslide = 1.0;

  change_height = false;
  change_look = false;
  rotate_around = false;
  azimuth = 0.0;
  elevation = 0.0;

  frame_rate = true;
  osgPlume = false;
  
  visual_layer = -1;

  
  // Set particle visual state to point based particles initially
  particle_visual_state = PARTICLE_SPHERE;

  //This shader is used to make final changes before rendering to the screen
  sphereParticle_shader.addShader("Shaders/sphereVisualize_vp.glsl", GLSLObject::VERTEX_SHADER);
  sphereParticle_shader.addShader("Shaders/sphereVisualize_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  sphereParticle_shader.createProgram();
  //  Field shader
  windField_shader.addShader("Shaders/windFieldLayer_vp.glsl", GLSLObject::VERTEX_SHADER);
  windField_shader.addShader("Shaders/windFieldLayer_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  windField_shader.createProgram();
  // Lighting shader
  lighting_shader.addShader("Shaders/lighting_vp.glsl", GLSLObject::VERTEX_SHADER);
  lighting_shader.addShader("Shaders/lighting_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  lighting_shader.createProgram();
  
  uniform_windTex = windField_shader.createUniform("Wind");
  uniform_max_velocity = windField_shader.createUniform("max_vel");
 
  // For the Sphere Particle Shader, we need the following data:
  // uniform variables for the texture units that hold the point
  // sprite and the normal map
  uniform_pointsprite_tex = sphereParticle_shader.createUniform("pointsprite_texunit");
  uniform_normalmap_tex = sphereParticle_shader.createUniform("pointspritenormal_texunit");
  uniform_visualization_tex = sphereParticle_shader.createUniform("visualization_texunit");

  // and determine which visual to use
  // uniform_pointsprite_visuals = sphereParticle_shader.createUniform("point_visuals");
  uniform_nx = sphereParticle_shader.createUniform("nx");
  uniform_ny = sphereParticle_shader.createUniform("ny");
  uniform_nz = sphereParticle_shader.createUniform("nz");
  uniform_numInRow = sphereParticle_shader.createUniform("numInRow");

  // POINT_SPRITE
  //
  // Initialize the texture data used for sprite based particle
  // visuals.  This currently means creating two textures: one
  // representing the base shape/color of the particle (white circle),
  // and a second to represent a normal map for a sphere to provide
  // Phong based lighting to the particle representations in a shader.
  createPointSpriteTextures();

  // float quadratic[] =  { 0.0f, 0.0f, 0.1f };
  // glPointParameterfvARB( GL_POINT_DISTANCE_ATTENUATION_ARB, quadratic );
  float maxSize = 0.0f, currSize = 1.0f;
  glGetFloatv(GL_POINT_SIZE_MAX_ARB, &maxSize);
  glGetFloatv(GL_POINT_SIZE, &currSize);
  glPointParameterfARB(GL_POINT_SIZE_MAX_ARB, currSize + 3.0);
  glPointParameterfARB(GL_POINT_SIZE_MIN_ARB, 1.0);
  
  // Create a high resolution clock timer - only works on Linux, x86
  // systems.  The basic timer works on Windows.  Setting the argument
  // to true will have no affect on windows implementations.
  clock_timer = new Timer(true);

  graphics_time[0] = clock_timer->tic();
  HUP_display_update_time[0] = clock_timer->tic();
  estimated_rate = 0.0;

  perform_cpu_sort = false;
  
  // Set and initialize the NetworkManager.
#if !WIN32
  network.setMode((NetworkManager::Mode)network_mode);
  network.init();
#endif
  
  // Set up for calculating and visualizing the percentage
  // of shadow for each cell.
  drawISD = false;
  inShadowData = new GLfloat[nz * nx * ny * 4];
  inShadowData2 = new GLfloat[nz * nx * ny * 4];
  
  // Setup the matrix to hold the camera model view matrix
  // (for the lighting shader).
  cameraModelviewMatrix = new GLfloat[16];
  
}

void DisplayControl::initTreadport() {
#if !WIN32
  // Check to make sure that we are in the right modes, currently
  // only the computer that is set to broadcast will communicate
  // with the treadport system.
  if(viewingMode == TREADPORT & network_mode == network.BROADCAST) {

    // Create the treadport object, this is done dynamically because
    // we attempt to connect to the treadport within the constructor
    // of TreadportManager.
    treadport = new TreadportManager();

    //
    // Set the initial values.
    //
    // Note that currently, their is no conversion nessisary to go
    // from the gpuPlume coordinate system to the treadport coordinate
    // system. However, if that were needed to be done it could be done
    // either here (and when we get the values back) OR it could be done
    // within the TreadportManager it self.
    
    std::vector<float> init_pos;
    init_pos.resize(3);
    init_pos[0] = eye_pos[0];
    init_pos[1] = eye_pos[1];
    init_pos[2] = eye_pos[2];
    
    std::vector<float> init_gaze;
    init_gaze.resize(3);
    init_gaze[0] = eye_gaze[0];
    init_gaze[1] = eye_gaze[1];
    init_gaze[2] = eye_gaze[2];

    std::vector<float> init_up;
    init_up.resize(3);
    init_up[0] = 0.0;
    init_up[1] = 0.0;
    init_up[2] = 1.0;

    // Call the initialization function within TreadportManager
    // to set up the initial values with the treadport system.
    treadport->init(init_pos, init_gaze, init_up);
  }
#endif
}

DisplayControl::~DisplayControl() {
  delete [] inShadowData;
  delete [] inShadowData2;
#if !WIN32
  delete treadport;
  delete [] cameraModelviewMatrix;
#endif
}

void DisplayControl::setEmitter(ParticleEmitter* p)
{
  pe = p;
}

void DisplayControl::setVisualPlane(VisualPlane* vp)
{
  plane = vp;
}

void DisplayControl::enableLightingShader() {
  
  glGetFloatv(GL_MODELVIEW_MATRIX, cameraModelviewMatrix);
  
  lighting_shader.activate();
  
  GLint shader_shadowTexture = lighting_shader.createUniform("shadowMap");
  glUniform1i(shader_shadowTexture, 1);
  
  GLint shader_boundTexture = lighting_shader.createUniform("boundTexture");
  glUniform1i(shader_boundTexture, 0);
  
  GLint shader_lightPosition = lighting_shader.createUniform("lightPosition");
  glUniform3f(shader_lightPosition, sun_pos[0],  sun_pos[1],  sun_pos[2]);
  
  GLint shaderID = lighting_shader.createUniform("lightModelviewMatrix");
  glUniformMatrix4fv(shaderID, 1, GL_FALSE, sunModelviewMatrix);
  
  shaderID = lighting_shader.createUniform("lightProjectionMatrix");
  glUniformMatrix4fv(shaderID, 1, GL_FALSE, sunProjectionMatrix);
    
  shaderID = lighting_shader.createUniform("lightScaleAndBiasMatrix");
  glUniformMatrix4fv(shaderID, 1, GL_FALSE, (GLfloat*)&sunScaleAndBiasMatrix);
  
  shaderID = lighting_shader.createUniform("cameraModelviewMatrix");
  glUniformMatrix4fv(shaderID, 1, GL_FALSE, cameraModelviewMatrix);
  
  CheckErrorsGL("DisplayControl::enableLightingShader: End of method.");
  
}

void DisplayControl::drawVisuals(GLuint vertex_buffer, GLuint texid3, GLuint color_buffer, 
				 int numInRow, int twidth, int theight, GLuint PositionTexId, GLuint VelTexId)
{

  // Timer_t displayStart = dTimer->tic();    
  
  DrawSkyBox(eye_pos[0], eye_pos[1], eye_pos[2], 50.0, 50.0, 50.0);
  
  drawAxes();
  
  if(draw_buildings) {
    drawFeatures();
  }
  
  enableLightingShader();
  drawGround();
  lighting_shader.deactivate();
  
  if(drawISD) {
    drawInShadowData();
  }
  
  // render the vertices in the VBO (the particle positions) as points in the domain
  
  if(color_buffer != 0) {
    glEnableClientState(GL_COLOR_ARRAY); 
    glBindBufferARB(GL_ARRAY_BUFFER, color_buffer);
    glColorPointer(4, GL_FLOAT, 0, 0);
  }
  
  glBindBufferARB(GL_ARRAY_BUFFER, vertex_buffer);
  glVertexPointer(4, GL_FLOAT, 0, 0);
  
  glEnableClientState(GL_VERTEX_ARRAY); 
  
  // see what the cost is for sorting the vertex buffer elements based on distance from the eye
  if (perform_cpu_sort)
  {

    glBindBufferARB(GL_ARRAY_BUFFER, vertex_buffer);
    GLfloat* vertex_data = (GLfloat*)glMapBuffer(GL_ARRAY_BUFFER, GL_READ_WRITE);

    if (vertex_data != (GLfloat*)NULL)
    {
      memcpy(&sort_eye, &eye_pos, sizeof(float)*3);
      qsort(vertex_data, twidth*theight, sizeof(GLfloat)*4, cmpVertices);
    }

    glUnmapBuffer(GL_ARRAY_BUFFER);
    
    // shut it back off
    // perform_cpu_sort = false;
  }

  sphereParticle_shader.activate();  

  glPointSize(6.0);
  glUniform1iARB(uniform_nx, nx);
  glUniform1iARB(uniform_ny, ny);
  glUniform1iARB(uniform_nz, nz);
  glUniform1iARB(uniform_numInRow, numInRow);

  glEnable(GL_TEXTURE_2D);
  glEnable(GL_TEXTURE_RECTANGLE_ARB);

  glActiveTextureARB(GL_TEXTURE0_ARB);
  glUniform1iARB(uniform_pointsprite_tex, 0);
  glBindTexture(GL_TEXTURE_2D, displayTex[0]); // point_sprite_textures[0]);

  glActiveTextureARB(GL_TEXTURE1_ARB);
  glUniform1iARB(uniform_normalmap_tex, 1);
  glBindTexture(GL_TEXTURE_2D, point_sprite_textures[1]);

  glActiveTexture(GL_TEXTURE2);
  glUniform1iARB(uniform_visualization_tex, 2);
  glBindTexture(GL_TEXTURE_RECTANGLE_ARB, texid3);

  glActiveTextureARB(GL_TEXTURE0_ARB);

  // all of our particle rendering methods now use the imposters
  glTexEnvf(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);
  glEnable(GL_POINT_SPRITE_ARB);
  glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

  if(color_buffer == 0)
    glColor4f(1.0,1.0,1.0,1.0);

  glDrawArrays(GL_POINTS, 0, twidth*theight);

  // all of our particle rendering methods now use the imposters
  glDisable(GL_POINT_SPRITE_ARB);
  glDisable(GL_TEXTURE_2D);

  sphereParticle_shader.deactivate();

#if 0
  // Sample code to draw the particle positions...
    {
      glMatrixMode(GL_PROJECTION);
      glPushMatrix();
      glLoadIdentity();
      gluOrtho2D(-1, 1, -1, 1);

      glMatrixMode(GL_MODELVIEW);
      glPushMatrix();
      glLoadIdentity();

      glEnable(texType);

      glUniform1iARB(uniform_PosTexSampler, 0);
      glUniform1iARB(uniform_DomainX, nx);
      glUniform1iARB(uniform_DomainY, ny);
      glUniform1iARB(uniform_DomainZ, nz);
      glUniform1iARB(uniform_DoNorm, 0);

      glBindTexture(texType, PositionTexId);
      glBegin(GL_QUADS);
      {
	glTexCoord2f(0, 0);			glVertex3f(-1, -1, -0.5f);
	glTexCoord2f(twidth, 0);		glVertex3f( 0, -1, -0.5f);
	glTexCoord2f(twidth, theight);          glVertex3f( 0,  1, -0.5f);
	glTexCoord2f(0, theight);		glVertex3f(-1,  1, -0.5f);
      }
      glEnd();

      glUniform1iARB(uniform_DomainX, 1);
      glUniform1iARB(uniform_DomainY, 1);
      glUniform1iARB(uniform_DomainZ, 1);
      glUniform1iARB(uniform_DoNorm, 1);

      glBindTexture(texType, VelTexId);
      glBegin(GL_QUADS);
      {
	glTexCoord2f(0, 0);			glVertex3f(0, -1, -0.5f);
	glTexCoord2f(twidth, 0);		glVertex3f( 1, -1, -0.5f);
	glTexCoord2f(twidth, theight);          glVertex3f( 1,  1, -0.5f);
	glTexCoord2f(0, theight);		glVertex3f(0,  1, -0.5f);
      }
      glEnd();

      glDisable(texType);

      glMatrixMode(GL_PROJECTION);
      glPopMatrix();

      glMatrixMode(GL_MODELVIEW);
      glPopMatrix();
    }

#endif

#if 0
  // spit out frame rate
  if(!osgPlume){
    if (frame_rate){
      drawFrameRate(twidth, theight);
    }
  }
#endif

  // Timer_t displayEnd = dTimer->tic();      
  // std::cout << "DC Display Time: " << dTimer->deltau(displayStart, displayEnd) << " us." << std::endl;  

#if 0
  // Display wind data for testing  
  sprintf(text_buffer, "Wind: %f %f %f", windDir[0], windDir[1], windDir[2]);
  OpenGLText(5, 5, text_buffer);
#endif

}

void DisplayControl::increaseVisualLayer()
{
  visual_layer++;
  if(visual_layer > nz) visual_layer = nz;

}

void DisplayControl::decreaseVisualLayer()
{
  visual_layer--;
  if(visual_layer < -1) visual_layer = -1;
}

void DisplayControl::moveForwardorBack(float change)
{
  //eye_pos[0] = eye_pos[0] + change;
  //eye_gaze[0] = eye_pos[0] - 5.0;
  if(change < 0.0){
    eye_pos[0] -= eye_gaze[0];
    eye_pos[1] -= eye_gaze[1];
    eye_pos[2] -= eye_gaze[2];
  }
  else{
    eye_pos[0] += eye_gaze[0];
    eye_pos[1] += eye_gaze[1];
    eye_pos[2] += eye_gaze[2];
  }
}

void DisplayControl::slideLeftorRight(float direction)
{

  if(direction < 0.0){
    eye_pos[0] -= xslide;
    eye_pos[1] -= yslide;
  }
  else{
    eye_pos[0] += xslide;
    eye_pos[1] += yslide;
  }

}

void DisplayControl::lookUporDown(float change)
{

  if(change < 0.0){
    pitch -= M_PI_2/90.0;
    //yaw = 0.0;
    //roll = 0.0;
  }
  else{
    pitch += M_PI_2/90.0;
    //yaw = 0.0;
    //roll = 0.0;
  }

  /*if(change < 0.0)
    yangle = yangle + (M_PI/90.0);
  else
    yangle = yangle - (M_PI/90.0);

  if(yangle > M_PI/2.0)
    yangle = yangle - (M_PI/90.0);
  else if(yangle < -M_PI/2.0)
    yangle = yangle + (M_PI/90.0);

    eye_gaze[2] = sin(yangle);*/

  if(pitch > M_PI_2)
    pitch = pitch - M_PI_2/90.0;
  else if(pitch < -M_PI_2)
    pitch = pitch + M_PI_2/90.0;

  calculateNormal();

  eye_gaze[0] = norm_x;
  eye_gaze[1] = norm_y;
  eye_gaze[2] = norm_z;

  //eye_gaze[1] = sin(angle);
  //eye_gaze[0] = cos(angle);

}

void DisplayControl::setAzimuth(float change, float rate)
{
  azimuth = azimuth + change*rate;
}

void DisplayControl::setElevation(float change, float rate)
{
  //elevation = elevation + change*rate;
  eye_pos[2] = eye_pos[2] + change*rate;
  //eye_gaze[2] = eye_gaze[2] + change*rate;
}

void DisplayControl::setRotateAround(float change)
{

  if(change < 0){
    yaw -= M_PI/90.0;
    //pitch = 0.0;
    //roll = 0.0;
  }
  else{
    yaw += M_PI/90.0;
    //pitch = 0.0;
    //roll = 0.0;
  }

  calculateNormal();

  /*if(change < 0)
    angle = angle + (M_PI/90.0);
  else 
    angle = angle - (M_PI/90.0);

  if(angle > 2*M_PI)
    angle = 0.0;
  else if(angle < 0.0)
    angle = (2*M_PI - M_PI/90.0);

  eye_gaze[0] = cos(angle);
  eye_gaze[1] = sin(angle);*/

  eye_gaze[0] = norm_x;
  eye_gaze[1] = norm_y;
  eye_gaze[2] = norm_z;

  //xslide = cos(angle-M_PI/2.0);
  // yslide = sin(angle-M_PI/2.0);
  

}

void DisplayControl::calculateNormal()
{

  float a11;//,a12,a13;
  float a21;//,a22,a23;
  float a31;//,a32,a33;
  
  //Transformation matrix
  a11 = cos(pitch)*cos(yaw);
  //a12 = sin(roll)*sin(pitch)*cos(yaw) + cos(roll)*sin(yaw);
  //a13 = (-cos(roll))*sin(pitch)*cos(yaw) + sin(roll)*sin(yaw);
  a21 = (-cos(pitch))*sin(yaw);
  //a22 = (-sin(roll))*sin(pitch)*sin(yaw) + cos(roll)*cos(yaw);
  //a23 = cos(roll)*sin(pitch)*sin(yaw) + sin(roll)*cos(yaw);
  a31 = sin(pitch);
  //a32 = (-sin(roll))*cos(pitch);
  //a33 = cos(roll)*cos(pitch);

  norm_x = a11*(-1.0);
  norm_y = a21*(-1.0);
  norm_z = a31*(-1.0);
  
  //norm_x = a11*norm_x + a12*norm_y + a13*norm_z;
  //norm_y = a21*norm_x + a22*norm_y + a23*norm_z;
  //norm_z = a31*norm_x + a32*norm_y + a33*norm_z;
  //normalize the normal
  float length;
  length = sqrt((norm_x*norm_x)+(norm_y*norm_y)+(norm_z*norm_z));
   
  norm_x = norm_x/length;
  norm_y = norm_y/length;
  norm_z = norm_z/length;

  float b11;//,a12,a13;
  float b21;//,a22,a23;
  //float b31;//,a32,a33;
  
  //Transformation matrix
  b11 = cos(pitch)*cos(yaw-M_PI_2);
  //a12 = sin(roll)*sin(pitch)*cos(yaw) + cos(roll)*sin(yaw);
  //a13 = (-cos(roll))*sin(pitch)*cos(yaw) + sin(roll)*sin(yaw);
  b21 = (-cos(pitch))*sin(yaw-M_PI_2);
  //a22 = (-sin(roll))*sin(pitch)*sin(yaw) + cos(roll)*cos(yaw);
  //a23 = cos(roll)*sin(pitch)*sin(yaw) + if(numSides[qi]==5)sin(roll)*cos(yaw);
  //b31 = sin(pitch);
  //a32 = (-sin(roll))*cos(pitch);
  //a33 = cos(roll)*cos(pitch);

  xslide = b11*(1.0);
  yslide = b21*(1.0);
    
  //norm_x = a11*norm_x + a12*norm_y + a13*norm_z;
  //norm_y = a21*norm_x + a22*norm_y + a23*norm_z;
  //norm_z = a31*norm_x + a32*norm_y + a33*norm_z;
  //normalize the normal
  
  length = sqrt((xslide*xslide)+(yslide*yslide));
   
  xslide = xslide/length;
  yslide = yslide/length;
  
}

void DisplayControl::drawSky()
{
  // double skyr = 0.4, skyg = 0.5, skyb = 0.8;
  double skyr = 0.5, skyg = 0.5, skyb = 0.6;
  glDisable(texType);

  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glDisable ( GL_TEXTURE_2D );
  glDepthMask(GL_FALSE);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
 
  //Put sky in background
  glBegin(GL_POLYGON);
  glColor4f(0.9,0.9,1.0, 1.0); //mid color
  glVertex3f(-1,-0.5,-1);
  glColor4f(0.9,0.9,1.0, 1.0); //white
  glVertex3f(-1,-1,-1);
  glVertex3f(1,-1,-1);
  glColor4f(0.9,0.9,1.0, 1.0); //mid color
  glVertex3f(1,-.5,-1);
  glColor4f(skyr, skyg, skyb, 1.0);
  glVertex3f(1,1,-1);
  glVertex3f(-1,1,-1);
  glEnd();
 
  glDepthMask(GL_TRUE);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glPopAttrib();

  glEnable(texType);

}

void DisplayControl::drawGround()
{

  glDisable(texType);
  glEnable(GL_TEXTURE_2D);
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D, shadowMap);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, displayTex[0]);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
  
  glColor4f(0.5,0.5,0.5,1.0);
  glBegin(GL_QUADS);
  {
    /*glTexCoord2f(0,0);       glVertex3f(0.0,0.0,-0.05);
    glTexCoord2f(1,0);       glVertex3f(nx,0.0,-0.05);
    glTexCoord2f(1,1);       glVertex3f(nx,ny,-0.05);
    glTexCoord2f(0,1);       glVertex3f(0.0,ny,-0.05); */

		glTexCoord2f(0,0);       glVertex3f(-nx*5,-ny*5,-0.05);
    glTexCoord2f(1,0);       glVertex3f(nx*5,-ny*5,-0.05);
    glTexCoord2f(1,1);       glVertex3f(nx*5,ny*5,-0.05);
    glTexCoord2f(0,1);       glVertex3f(-nx*5,ny*3,-0.05);
  }
  glEnd();

  lighting_shader.deactivate();

  glBindTexture(GL_TEXTURE_2D, 0);
  glDisable(GL_TEXTURE_2D);
  glEnable(texType);

}

void DisplayControl::drawGrid()
{
  GLint lwidth;

  glGetIntegerv(GL_LINE_WIDTH, &lwidth);

  glDisable(texType);

  glLineWidth(0.5);
  glBegin(GL_LINES);
  {
    glColor3f(0.4,0.4,0.4);
    for(int i=-ny; i <= ny*2; i+=5){
      glVertex3f(-nx,i,-0.15);
      glVertex3f(nx*2,i,-0.15);
    }
    for(int i=-nx; i <= nx*2; i+=5){
      glVertex3f(i,-ny,-0.15);
      glVertex3f(i,ny*2,-0.15);
    }

  
  }
  glEnd();

  glLineWidth(lwidth);
 
  glEnable(texType);
}

void DisplayControl::drawAxes()
{
  // query the current line width so we can set it back at the end of
  // the function
  //int numLines = ny + 10;
  glDisable(texType);
  GLint lwidth;
  glGetIntegerv(GL_LINE_WIDTH, &lwidth);
  glLineWidth(3.0);

  // Modified colors of axes to represent the manner in which
  // direction is being visualized.
  // 
  // Which makes me think we might want to make the axes' colors
  // dependent on what is being visualized.

  glBegin(GL_LINES);
  glColor3f(0.0, 0.0, 1.0);
  glVertex3f(0.0, 0.0, 0.0);
  glColor3f(1.0, 1.0, 0.0);
  glVertex3f(nx, 0.0, 0.0);

  glColor3f(0.0, 1.0, 0.0);
  glVertex3f(0.0, 0.0, 0.0);
  glColor3f(1.0, 0.0, 0.0);
  glVertex3f(0.0, ny, 0.0);

  glColor3f(0.0, 0.0, 0.0);
  glVertex3f(0.0, 0.0, 0.0);
  glColor3f(1.0, 1.0, 1.0);
  glVertex3f(0.0, 0.0, nz);
  glEnd();

  //This draws a label for the axis by drawing a textured quad.
  /*
  glDisable(texType);
  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D, axisLabel[0]);
  glPushMatrix();
  glRotatef(180.0, 0.0, 0.0, 1.0);

  glBegin(GL_QUADS);
  {
    glTexCoord2f(0, 0);      glVertex3f(-0.5, -0.5, nz+1.0);
    glTexCoord2f(1, 0);      glVertex3f(0.5, -0.5, nz+1.0);
    glTexCoord2f(1, 1);      glVertex3f(0.5, 0.5, nz+1.0);
    glTexCoord2f(0, 1);      glVertex3f(-0.5, 0.5, nz+1.0);
  }
  

  glEnd();

  glPopMatrix();
  glDisable(GL_TEXTURE_2D);
  glEnable(texType);*/

  // set the line width back to what it was
  glLineWidth(lwidth);
  glEnable(texType);

}

void DisplayControl::drawLayers(GLuint texId, int numInRow, float maxVel)
{
  if (visual_layer >= 0 && visual_layer < nz)
    {
      glPushMatrix();
      //glEnable(GL_LIGHTING);
      //glEnable(GL_LIGHT0);
      glEnable(GL_COLOR_MATERIAL);
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);     

      glEnable(texType);
      glBindTexture(texType, texId);
      glTexParameteri(texType, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(texType, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      
      // The s and t parameters reference pixels on the UVW texture
      // map.  Because we are using texture rectangles, the s and t
      // parameters do not need to be in normalized device coordinates
      // so the values represent actual pixels positions in the
      // texture map.
      int s = 0;
      int t = 0;

      // The texture we use here is the packed texture containing all
      // cells of the 3D uvw wind field data.  It is packed because we
      // flatten the 3D structure into a series of discrete 2D
      // elements. The variable numInRow is the number of these
      // discrete 2D layers that can fit in each row the texture. 

      // The coordinate frame we use is with Y up, so X and Z (at Y=0)
      // forms the ground plane.
      //

      // s (or the value in the x dimension of the texture) can be
      // determined with a mod of the layer by the number of layers
      // that can be held in each row of the texutre. [ COMPLETE DESCRIPTION ]
      s = (int)(visual_layer % numInRow) * nxdx;

      // t (or the value in the y dimension of the texture) can be 
      // calculated by the floor of the layer to be visualized divided
      // by the number of layers that can be held in each row of
      // the texture. 
      t = (int)(floor(visual_layer/(float)numInRow) * nydy);

      windField_shader.activate();

      glUniform1fARB(uniform_max_velocity, maxVel);
      //glUniform1iARB(uniform_controlWind, visual_field);
      //glUniform1fARB(uniform_sliderWind, slider);

      glUniform1iARB(uniform_windTex, 0); 

      // Create a quad at this layer with 50% transparency
      glColor4f(1.0, 1.0, 1.0, 0.8);
      glBegin(GL_QUADS);
      {
	glNormal3f(0.0, 0.0, 1.0);
	glTexCoord2f(s, t);         glVertex3f(0, 0, visual_layer);
	glTexCoord2f(s+nxdx, t);      glVertex3f(nx, 0,visual_layer);
	glTexCoord2f(s+nxdx, t+nydy);   glVertex3f(nx, ny,visual_layer);
	glTexCoord2f(s, t+nydy);      glVertex3f(0, ny,visual_layer);
      }
      glEnd();

      windField_shader.deactivate();

      // set back to nearest
      glTexParameteri(texType, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(texType, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

      //glBindTexture(texType, 0);
      glDisable(texType);
      glDisable(GL_BLEND);
      glDisable(GL_COLOR_MATERIAL);
      //glDisable(GL_LIGHT0);
      //glDisable(GL_LIGHTING);

      glPopMatrix();
    }
  
}

void instanceCube()
{
}

void DisplayControl::initVars(int nb,int* ns,float* x, float* y, float* z,
				    float* h, float* w, float* l,float* g)
{
  numBuild = nb;
  numSides = ns;
  xfo = x;
  yfo = y;
  zfo = z;
  ht = h;
  wti = w;
  lti = l;
  gamma=g;
  glDisable(texType);
  glEnable(GL_TEXTURE_2D);
  glGenTextures(3,displayTex);
  
  createImageTex(displayTex[0], (char *)"concrete.ppm");
  createImageTex(displayTex[1], (char *)"building.ppm");
  createImageTex(displayTex[2], (char *)"buildingRoof.ppm");

  glGenTextures(6,skyBoxTex);
/*  createImageTex(skyBoxTex[0], (char *)"SkyBox/front.ppm");
  createImageTex(skyBoxTex[1], (char *)"SkyBox/left.ppm");
  createImageTex(skyBoxTex[2], (char *)"SkyBox/back.ppm");
  createImageTex(skyBoxTex[3], (char *)"SkyBox/right.ppm");
  createImageTex(skyBoxTex[4], (char *)"SkyBox/up.ppm");
  createImageTex(skyBoxTex[5], (char *)"SkyBox/down.ppm"); */
	/*
		for the mystic skybox. 
	createImageTex(skyBoxTex[0], (char *)"SkyBox/mystic/mystic_east.ppm"); //front
  createImageTex(skyBoxTex[1], (char *)"SkyBox/mystic/mystic_north.ppm"); //left
  createImageTex(skyBoxTex[2], (char *)"SkyBox/mystic/mystic_west.ppm"); //back
  createImageTex(skyBoxTex[3], (char *)"SkyBox/mystic/mystic_south.ppm");//right
  createImageTex(skyBoxTex[4], (char *)"SkyBox/mystic/mystic_up.ppm");
  createImageTex(skyBoxTex[5], (char *)"SkyBox/mystic/mystic_down.ppm");
  */
	/* for skybox with the cloudy reef
	createImageTex(skyBoxTex[0], (char *)"SkyBox/clouds/reef_east.ppm"); //front
  createImageTex(skyBoxTex[1], (char *)"SkyBox/clouds/reef_north.ppm"); //left
  createImageTex(skyBoxTex[2], (char *)"SkyBox/clouds/reef_west.ppm"); //back
  createImageTex(skyBoxTex[3], (char *)"SkyBox/clouds/reef_south.ppm");//right
  createImageTex(skyBoxTex[4], (char *)"SkyBox/clouds/reef_up.ppm");
  createImageTex(skyBoxTex[5], (char *)"SkyBox/clouds/reef_down.ppm"); */

  createImageTex(skyBoxTex[0], (char *)"SkyBox/rays/rays_east.ppm"); //front
  createImageTex(skyBoxTex[1], (char *)"SkyBox/rays/rays_north.ppm"); //left
  createImageTex(skyBoxTex[2], (char *)"SkyBox/rays/rays_west.ppm"); //back
  createImageTex(skyBoxTex[3], (char *)"SkyBox/rays/rays_south.ppm");//right
  createImageTex(skyBoxTex[4], (char *)"SkyBox/rays/rays_up.ppm");
  createImageTex(skyBoxTex[5], (char *)"SkyBox/rays/rays_down.ppm");

 
  glDisable(GL_TEXTURE_2D);
  glEnable(texType);

}

void DisplayControl::DrawSkyBox(float x, float y, float z,  float width, float height, float length)
{
  glDisable(GL_DEPTH_TEST);
  glDisable(texType);
  glEnable(GL_TEXTURE_2D);
  
  // center the box around the postion paramaters
  x = x - width  / 2;
  y = y - height / 2;
  z = z - length / 2;
  
  glColor4f(1.0,1.0,1.0,1.0f); // set color to white
  
  // Draw the Front (in front of you)
  glBindTexture(GL_TEXTURE_2D, skyBoxTex[0]);
  glBegin(GL_QUADS);
    glTexCoord2f(0.0f, 1.0f); glVertex3f(x, y, z); //bottom left
    glTexCoord2f(0.0f, 0.0f); glVertex3f(x, y, z+height);
    glTexCoord2f(1.0f, 0.0f); glVertex3f(x, y+width, z+height);
    glTexCoord2f(1.0f, 1.0f); glVertex3f(x, y+width, z);
  glEnd();

  // Draw Right side
  glBindTexture(GL_TEXTURE_2D, skyBoxTex[3]);
  glBegin(GL_QUADS);		
    glTexCoord2f(0.0f, 0.0f); glVertex3f(x, y+width, z+height);
    glTexCoord2f(1.0f, 0.0f); glVertex3f(x+length, y+width, z+height);
    glTexCoord2f(1.0f, 1.0f); glVertex3f(x+length, y+width, z);
    glTexCoord2f(0.0f, 1.0f); glVertex3f(x, y+width, z);
  glEnd();

  // Draw Back (behind you)
  glBindTexture(GL_TEXTURE_2D, skyBoxTex[2]);
  glBegin(GL_QUADS);
    glTexCoord2f(1.0f, 0.0f); glVertex3f(x+length, y, z+height);
    glTexCoord2f(0.0f, 0.0f); glVertex3f(x+length, y+width, z+height);
    glTexCoord2f(0.0f, 1.0f); glVertex3f(x+length, y+width, z);
    glTexCoord2f(1.0f, 1.0f); glVertex3f(x+length, y, z); //from behind you (facing away)- bottom left.
  glEnd();
  
  // Draw the left side
  glBindTexture(GL_TEXTURE_2D, skyBoxTex[1]);
  glBegin(GL_QUADS);
    glTexCoord2f(0.0f, 0.0f); glVertex3f(x+length, y, z+height);
    glTexCoord2f(1.0f, 0.0f); glVertex3f(x, y, z+height);
    glTexCoord2f(1.0f, 1.0f); glVertex3f(x, y, z);
    glTexCoord2f(0.0f, 1.0f); glVertex3f(x+length, y, z);
  glEnd();

  // Draw the top
  glBindTexture(GL_TEXTURE_2D, skyBoxTex[4]);
  glBegin(GL_QUADS);
    glTexCoord2f(1.0f, 1.0f); glVertex3f(x+length, y+width, z+height); //back right
    glTexCoord2f(1.0f, 0.0f); glVertex3f(x+length, y, z+height); //back left
    glTexCoord2f(0.0f, 0.0f); glVertex3f(x, y, z+height);  //front left
    glTexCoord2f(0.0f, 1.0f); glVertex3f(x, y+width, z+height); ///front right
  glEnd();

  //draw the bottom
  glBindTexture(GL_TEXTURE_2D, skyBoxTex[5]);
  glBegin(GL_QUADS);
    glTexCoord2f(0.0f, 0.0f); glVertex3f(x, y+width, z); //front right
    glTexCoord2f(1.0f, 0.0f); glVertex3f(x+length, y+width, z); //back rigth
    glTexCoord2f(1.0f, 1.0f); glVertex3f(x+length, y, z); //back left
    glTexCoord2f(0.0f, 1.0f); glVertex3f(x, y, z);// front left
  glEnd();
	
  glDisable(GL_TEXTURE_2D);
  glEnable(texType);
  
  if(drawISD) {
    glPushMatrix();
      glTranslatef(eye_pos[0] + sun_pos[0], eye_pos[1] + sun_pos[1], eye_pos[2] + sun_pos[2]);
      glColor4f(1.0f, 1.0f, 0.0f, 1.0f);
      glutSolidSphere(5, 20, 20);
    glPopMatrix();
  }
  
  glEnable(GL_DEPTH_TEST);
} // end DrawSkyBox

void DisplayControl::createImageTex(GLuint texture, char* filename)
{
  GLubyte* testImage;
  int w, h;

  testImage = readPPM(filename, &w, &h);
  if(testImage == 0){
    std::cout << "Didn't Load Texture File" << std::endl;
  }
  glBindTexture(GL_TEXTURE_2D, texture); 
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, testImage);
}


void DisplayControl::drawFeatures(void)
{
  //float grid_scale = 1.0;  // currently, just 1 but likely needs to come from QUICPLUME
  glDisable(texType);
  glEnable(GL_TEXTURE_2D);
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D, shadowMap);
  glActiveTexture(GL_TEXTURE0);
  // glEnable(GL_COLOR_MATERIAL);
  // glEnable(GL_LIGHTING);
  // glEnable(GL_LIGHT0);
  
  // Draw the building
  for (int qi=0; qi<numBuild; qi++) 
    {
      
      if(numSides[qi]==4)
      {
      glPushMatrix();
      glTranslatef(xfo[qi],yfo[qi],zfo[qi]);
      glRotatef((GLfloat)(gamma[qi]), 0.0, 0.0, 1.0);
      glTranslatef(-xfo[qi],-yfo[qi],-zfo[qi]);
      glColor4f(1.0, 1.0, 1.0,1.0);

      glBindTexture(GL_TEXTURE_2D, displayTex[1]);

      glBegin(GL_QUADS);
      {
	glTexCoord2f(0,0);       glVertex3f(xfo[qi],yfo[qi]-(wti[qi]/2.0),zfo[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]+lti[qi],yfo[qi]-(wti[qi]/2.0),zfo[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]+lti[qi],yfo[qi]-(wti[qi]/2.0),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi],yfo[qi]-(wti[qi]/2.0),zfo[qi]+ht[qi]);
	
	glTexCoord2f(0,0);       glVertex3f(xfo[qi],yfo[qi]+(wti[qi]/2.0),zfo[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]+lti[qi],yfo[qi]+(wti[qi]/2.0),zfo[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]+lti[qi],yfo[qi]+(wti[qi]/2.0),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi],yfo[qi]+(wti[qi]/2.0),zfo[qi]+ht[qi]);
	
	glTexCoord2f(0,0);       glVertex3f(xfo[qi],yfo[qi]+(wti[qi]/2.0),zfo[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi],yfo[qi]-(wti[qi]/2.0),zfo[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi],yfo[qi]-(wti[qi]/2.0),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi],yfo[qi]+(wti[qi]/2.0),zfo[qi]+ht[qi]);

	glTexCoord2f(0,0);       glVertex3f(xfo[qi]+lti[qi],yfo[qi]+(wti[qi]/2.0),zfo[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]+lti[qi],yfo[qi]-(wti[qi]/2.0),zfo[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]+lti[qi],yfo[qi]-(wti[qi]/2.0),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi]+lti[qi],yfo[qi]+(wti[qi]/2.0),zfo[qi]+ht[qi]);
      }
      glEnd();

      glBindTexture(GL_TEXTURE_2D,displayTex[2]);
      glBegin(GL_QUADS);
      {
	glTexCoord2f(0,0);	glVertex3f(xfo[qi],yfo[qi]-(wti[qi]/2.0),zfo[qi]+ht[qi]);
	glTexCoord2f(1,0);	glVertex3f(xfo[qi]+lti[qi],yfo[qi]-(wti[qi]/2.0),zfo[qi]+ht[qi]);
	glTexCoord2f(1,1);	glVertex3f(xfo[qi]+lti[qi],yfo[qi]+(wti[qi]/2.0),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);	glVertex3f(xfo[qi],yfo[qi]+(wti[qi]/2.0),zfo[qi]+ht[qi]);
		
      }
      glEnd();
      glPopMatrix();
     }

    //drawing a pentagon
     else if(numSides[qi]==5)
     {
        glPushMatrix();
        glTranslatef(xfo[qi],yfo[qi],zfo[qi]);
      glRotatef((GLfloat)(gamma[qi]), 0.0, 0.0, 1.0);
      glTranslatef(-xfo[qi],-yfo[qi],-zfo[qi]);
        glBindTexture(GL_TEXTURE_2D, displayTex[1]);
        //glTranslatef(3,0,0);

      glBegin(GL_QUADS);
      {
	glTexCoord2f(0,0);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((18.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((18.0*3.14)/180.0)),zfo[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((-54.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((-54.0*3.14)/180.0)),zfo[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((-54.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((-54.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((18.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((18.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	
	glTexCoord2f(0,0);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((-54.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((-54.0*3.14)/180.0)),zfo[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((-126.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((-126.0*3.14)/180.0)),zfo[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((-126.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((-126.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((-54.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((-54.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	
	glTexCoord2f(0,0);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((-126.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((-126.0*3.14)/180.0)),zfo[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((162.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((162.0*3.14)/180.0)),zfo[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((162.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((162.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((-126.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((-126.0*3.14)/180.0)),zfo[qi]+ht[qi]);

	glTexCoord2f(0,0);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((162.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((162.0*3.14)/180.0)),zfo[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((90.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((90.0*3.14)/180.0)),zfo[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((90.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((90.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((162.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((162.0*3.14)/180.0)),zfo[qi]+ht[qi]);

        glTexCoord2f(0,0);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((90.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((90.0*3.14)/180.0)),zfo[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((18.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((18.0*3.14)/180.0)),zfo[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((18.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((18.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((90.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((90.0*3.14)/180.0)),zfo[qi]+ht[qi]);

      }
      glEnd();

      glBegin(GL_QUADS);
      {
	glTexCoord2f(0,0);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((18.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((18.0*3.14)/180.0)),zfo[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((-54.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((-54.0*3.14)/180.0)),zfo[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((-54.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((-54.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((18.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((18.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	
	glTexCoord2f(0,0);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((-54.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((-54.0*3.14)/180.0)),zfo[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((-126.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((-126.0*3.14)/180.0)),zfo[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((-126.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((-126.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((-54.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((-54.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	
	glTexCoord2f(0,0);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((-126.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((-126.0*3.14)/180.0)),zfo[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((162.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((162.0*3.14)/180.0)),zfo[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((162.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((162.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((-126.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((-126.0*3.14)/180.0)),zfo[qi]+ht[qi]);

	glTexCoord2f(0,0);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((162.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((162.0*3.14)/180.0)),zfo[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((90.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((90.0*3.14)/180.0)),zfo[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((90.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((90.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((162.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((162.0*3.14)/180.0)),zfo[qi]+ht[qi]);

        glTexCoord2f(0,0);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((90.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((90.0*3.14)/180.0)),zfo[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((18.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((18.0*3.14)/180.0)),zfo[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((18.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((18.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((90.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((90.0*3.14)/180.0)),zfo[qi]+ht[qi]);

      }
      glEnd();
 

      glBindTexture(GL_TEXTURE_2D, displayTex[2]);
      glBegin(GL_QUADS);
      {
	
	glTexCoord2f(0,0);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((18.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((18.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((18.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((18.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((-54.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((-54.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((-54.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((-54.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	
	glTexCoord2f(0,0);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((-54.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((-54.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((-54.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((-54.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((-126.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((-126.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((-126.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((-126.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	
	glTexCoord2f(0,0);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((-126.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((-126.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((-126.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((-126.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((162.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((162.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((162.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((162.0*3.14)/180.0)),zfo[qi]+ht[qi]);

	glTexCoord2f(0,0);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((162.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((162.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((162.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((162.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((90.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((90.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((90.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((90.0*3.14)/180.0)),zfo[qi]+ht[qi]);

        glTexCoord2f(0,0);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((90.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((90.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((90.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((90.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]+((wti[qi]/5.0)*cos((18.0*3.14)/180.0)),yfo[qi]+((wti[qi]/5.0)*sin((18.0*3.14)/180.0)),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi]+((wti[qi]/2.0)*cos((18.0*3.14)/180.0)),yfo[qi]+((wti[qi]/2.0)*sin((18.0*3.14)/180.0)),zfo[qi]+ht[qi]);

      }
      glEnd();
      glPopMatrix();
      //glTranslatef(-3,0,0);

     }

     //drawing a hexagon
     else if(numSides[qi]==6)
     {
        glPushMatrix();
        glTranslatef(xfo[qi],yfo[qi],zfo[qi]);
      glRotatef((GLfloat)(gamma[qi]), 0.0, 0.0, 1.0);
      glTranslatef(-xfo[qi],-yfo[qi],-zfo[qi]);
        glBindTexture(GL_TEXTURE_2D, displayTex[1]);
        //glTranslatef(3,0,0);

      glBegin(GL_QUADS);
      {
	glTexCoord2f(0,0);       glVertex3f(xfo[qi]-(lti[qi]/4.0),yfo[qi]-(wti[qi]/2.0),zfo[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]+(lti[qi]/4.0),yfo[qi]-(wti[qi]/2.0),zfo[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]+(lti[qi]/4.0),yfo[qi]-(wti[qi]/2.0),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi]-(lti[qi]/4.0),yfo[qi]-(wti[qi]/2.0),zfo[qi]+ht[qi]);
	
	glTexCoord2f(0,0);       glVertex3f(xfo[qi]+(lti[qi]/4.0),yfo[qi]-(wti[qi]/2.0),zfo[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]+(lti[qi]/2.0),yfo[qi],zfo[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]+(lti[qi]/2.0),yfo[qi],zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi]+(lti[qi]/4.0),yfo[qi]-(wti[qi]/2.0),zfo[qi]+ht[qi]);
	
	glTexCoord2f(0,0);       glVertex3f(xfo[qi]+(lti[qi]/2.0),yfo[qi],zfo[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]+(lti[qi]/4.0),yfo[qi]+(wti[qi]/2.0),zfo[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]+(lti[qi]/4.0),yfo[qi]+(wti[qi]/2.0),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi]+(lti[qi]/2.0),yfo[qi],zfo[qi]+ht[qi]);

	glTexCoord2f(0,0);       glVertex3f(xfo[qi]+(lti[qi]/4.0),yfo[qi]+(wti[qi]/2.0),zfo[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]-(lti[qi]/4.0),yfo[qi]+(wti[qi]/2.0),zfo[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]-(lti[qi]/4.0),yfo[qi]+(wti[qi]/2.0),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi]+(lti[qi]/4.0),yfo[qi]+(wti[qi]/2.0),zfo[qi]+ht[qi]);

        glTexCoord2f(0,0);       glVertex3f(xfo[qi]-(lti[qi]/4.0),yfo[qi]+(wti[qi]/2.0),zfo[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]-(lti[qi]/2.0),yfo[qi],zfo[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]-(lti[qi]/2.0),yfo[qi],zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi]-(lti[qi]/4.0),yfo[qi]+(wti[qi]/2.0),zfo[qi]+ht[qi]);

        glTexCoord2f(0,0);       glVertex3f(xfo[qi]-(lti[qi]/2.0),yfo[qi],zfo[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]-(lti[qi]/4.0),yfo[qi]-(wti[qi]/2.0),zfo[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]-(lti[qi]/4.0),yfo[qi]-(wti[qi]/2.0),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi]-(lti[qi]/2.0),yfo[qi],zfo[qi]+ht[qi]);
	
      }
      glEnd();

      glBindTexture(GL_TEXTURE_2D, displayTex[2]);
      glBegin(GL_POLYGON);
      {
	
	glTexCoord2f(0.16,0);       glVertex3f(xfo[qi]-(lti[qi]/4.0),yfo[qi]-(wti[qi]/2.0),zfo[qi]+ht[qi]);

	glTexCoord2f(0.83,0);       glVertex3f(xfo[qi]+(lti[qi]/4.0),yfo[qi]-(wti[qi]/2.0),zfo[qi]+ht[qi]);
	
	glTexCoord2f(1,0.5);       glVertex3f(xfo[qi]+(lti[qi]/2.0),yfo[qi],zfo[qi]+ht[qi]);
	
	glTexCoord2f(0.83,1);        glVertex3f(xfo[qi]+(lti[qi]/4.0),yfo[qi]+(wti[qi]/2.0),zfo[qi]+ht[qi]);

	glTexCoord2f(0.16,1);       glVertex3f(xfo[qi]-(lti[qi]/4.0),yfo[qi]+(wti[qi]/2.0),zfo[qi]+ht[qi]);

        glTexCoord2f(0,0.5);       glVertex3f(xfo[qi]-(lti[qi]/2.0),yfo[qi],zfo[qi]+ht[qi]);

      }
      glEnd();
      glPopMatrix();
      //glTranslatef(-3,0,0);

     }
    //drawing a cylinder 
    else
    {
         
         float t1,t2;
         glPushMatrix();
         glTranslatef(xfo[qi],yfo[qi],zfo[qi]);
      glRotatef((GLfloat)(gamma[qi]), 0.0, 0.0, 1.0);
      glTranslatef(-xfo[qi],-yfo[qi],-zfo[qi]);
         glBindTexture(GL_TEXTURE_2D, displayTex[1]);
         glBegin(GL_QUADS);
         {
	      float incr = 10.0;
              for(float Theta=0.0;Theta<=360.0;Theta+=incr)
              {
                 t1=float(Theta/360.0);
                 t2=float((Theta+incr)/360.0);
                 //printf(" \n %f %f %f %f\n",xfo[qi]-((lti[qi]/2.0)*cos(Theta)),yfo[qi]-((lti[qi]/2.0)*sin(Theta)),xfo[qi]-((lti[qi]/2.0)*cos(Theta+incr)),yfo[qi]-((lti[qi]/2.0)*sin(Theta+incr)));
                 //printf(" \n %f %f\n",cos(Theta+incr),sin(Theta+incr));
                 glTexCoord2f(t1,0);       glVertex3f(xfo[qi]-((lti[qi]/2.0)*cos(Theta*(3.14/180.0))-lti[qi]/2.0),yfo[qi]-((wti[qi]/2.0)*sin(Theta*(3.14/180.0))),zfo[qi]);
	         glTexCoord2f(t2,0);       glVertex3f(xfo[qi]-((lti[qi]/2.0)*cos((Theta+incr)*(3.14/180.0))-lti[qi]/2.0),yfo[qi]-((wti[qi]/2.0)*sin((Theta+incr)*(3.14/180.0))),zfo[qi]);
	         glTexCoord2f(t2,1);       glVertex3f(xfo[qi]-((lti[qi]/2.0)*cos((Theta+incr)*(3.14/180.0))-lti[qi]/2.0),yfo[qi]-((wti[qi]/2.0)*sin((Theta+incr)*(3.14/180.0))),zfo[qi]+ht[qi]); 
	         glTexCoord2f(t1,1);       glVertex3f(xfo[qi]-((lti[qi]/2.0)*cos(Theta*(3.14/180.0))-lti[qi]/2.0),yfo[qi]-((wti[qi]/2.0)*sin(Theta*(3.14/180.0))),zfo[qi]+ht[qi]);               
              }
         }
         glEnd();

         glBindTexture(GL_TEXTURE_2D, displayTex[2]);
         glBegin(GL_TRIANGLES);
         {
              float incr = 10.0;
              for(float Theta=0;Theta<=360;Theta+=incr)
              {
                 glTexCoord2f(0,0);       glVertex3f(xfo[qi]-((lti[qi]/2.0)*cos(Theta*(3.14/180.0))-lti[qi]/2.0),yfo[qi]-((wti[qi]/2.0)*sin(Theta*(3.14/180.0))),zfo[qi]+ht[qi]);
	         glTexCoord2f(1,0);       glVertex3f(xfo[qi],yfo[qi],zfo[qi]+ht[qi]);
	         glTexCoord2f(0.5,0.5);   glVertex3f(xfo[qi]-((lti[qi]/2.0)*cos((Theta+incr)*(3.14/180.0))-lti[qi]/2.0),yfo[qi]-((wti[qi]/2.0)*sin((Theta+incr)*(3.14/180.0))),zfo[qi]+ht[qi]); 
	                       
              }
         }
         glEnd();
         glPopMatrix();
    }
  
  // glPopMatrix();

    }
	
  // glDisable(GL_LIGHT0);
  // glDisable(GL_LIGHTING);
  glDisable(GL_TEXTURE_2D);
	
  glEnable(texType);
	
}

//text: draws a string of text with an 18 point helvetica bitmap font
// at position (x,y) in window space(bottom left corner is (0,0).

void DisplayControl::drawFrameRate(int twidth, int theight)
{
  static bool do_first_time_only = true;
  double alpha = 0.125;
  double sampled_rate;

  // record end clock time
  graphics_time[1] = clock_timer->tic();

  sampled_rate = clock_timer->deltam( graphics_time[0], graphics_time[1] );

  if (do_first_time_only) 
    {
      estimated_rate = sampled_rate;
      last_estimated_rate = sampled_rate;
      do_first_time_only = false;
    }
  else 
    estimated_rate = (1.0 - alpha) * estimated_rate + alpha * sampled_rate;

  //
  // update the time on the screen every so often, but sample it more often
  //
  HUP_display_update_time[1] = clock_timer->tic();
  if (clock_timer->deltas( HUP_display_update_time[0], HUP_display_update_time[1] ) > 0.25)
    {
      last_estimated_rate = estimated_rate;
      HUP_display_update_time[0] = clock_timer->tic();
    }

  sprintf(text_buffer, "%d particles, Sim Step Time: %.2f ms", twidth*theight, last_estimated_rate);
  OpenGLText(5, 5, text_buffer);

  // record start clock time
  graphics_time[0] = clock_timer->tic();
}

void DisplayControl::OpenGLText(int x, int y, char* s)
{
  int lines;
  char* p;
  GLint* vp = new GLint[4];
  glGetIntegerv(GL_VIEWPORT,vp);

  // glDisable(GL_LIGHTING);
  glDisable(texType);
  glDisable(GL_DEPTH_TEST);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0, vp[2], 
	  0, vp[3], -1, 1);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glColor3ub(0, 0, 0);
  glRasterPos2i(x+1, y-1);
  for (p=s, lines=0; *p; p++) {
    if (*p == '\n') {
      lines++;
      glRasterPos2i(x+1, y-1-(lines*18));
    }
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *p);
  }
  glColor3ub(255, 255, 0);
  glRasterPos2i(x, y);
  for (p=s, lines=0; *p; p++) {
    if (*p == '\n') {
      lines++;
      glRasterPos2i(x, y-(lines*18));
    }
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *p);
  }
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glEnable(GL_DEPTH_TEST);
  glEnable(texType);
  // glEnable(GL_LIGHTING);

}

GLubyte* DisplayControl::readPPM(char* filename, int* width, int* height)
{
    FILE* fp;
    int i, w, h, d;
    unsigned char* image;
    char head[70];          /* max line <= 70 in PPM (per spec). */
    
    fp = fopen(filename, "rb");
    if (!fp) {
        perror(filename);
        return NULL;
    }
    
    /* grab first two chars of the file and make sure that it has the
       correct magic cookie for a raw PPM file. */
    fgets(head, 70, fp);
    if (strncmp(head, "P6", 2)) {
        fprintf(stderr, "%s: Not a raw PPM file\n", filename);
        return NULL;
    }
    
    /* grab the three elements in the header (width, height, maxval). */
    i = 0;
    while(i < 3) {
        fgets(head, 70, fp);
        if (head[0] == '#')     /* skip comments. */
            continue;
        if (i == 0)
            i += sscanf(head, "%d %d %d", &w, &h, &d);
        else if (i == 1)
            i += sscanf(head, "%d %d", &h, &d);
        else if (i == 2)
            i += sscanf(head, "%d", &d);
    }
    
    /* grab all the image data in one fell swoop. */
    image = (unsigned char*)malloc(sizeof(unsigned char)*w*h*3);
    fread(image, sizeof(unsigned char), w*h*3, fp);
    fclose(fp);
    
    *width = w;
    *height = h;
    return image;
}

void DisplayControl::createPointSpriteTextures()
{
  unsigned int radius = 32;
  unsigned int width = 128;
  unsigned int height = 128;

  // for this, I want width & height to be the same
  assert(width == height);
  GLfloat *data = new GLfloat[ width * height * 4 ];

  glEnable(GL_TEXTURE_2D);
  glGenTextures(2, point_sprite_textures);

  glBindTexture(GL_TEXTURE_2D, point_sprite_textures[0]);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  // make the shape/intensity map
  for (unsigned int y=0; y<height; y++)
    for (unsigned int x=0; x<width; x++)
      {
	int idx = y*width*4 + x*4;

	// convert to sphere coordinates with radius half of width/height
	double r = (double)radius;
	double xr = x - (double)width/2.0;
	double yr = y - (double)height/2.0;

	// if coordinate is on or inside the circle, keep going
	double mag = sqrt(xr*xr + yr*yr);
	if (mag <= r)
	  {
	    // White for now... may use later
	    data[idx] = 1.0; data[idx+1] = 1.0; data[idx+2] = 1.0; data[idx+3] = 1.0;
	  }
	else if ((mag >= r+(radius*0.20)) && (mag <= r+(radius*0.40)))
	  {
	    // store transparent white in image space
	    data[idx] = 1.0; data[idx+1] = data[idx+2] = 0.0;
	    data[idx+3] = 1.0;
	  }
	else 
	  {
	    // store transparent white in image space
	    data[idx] = data[idx+1] = data[idx+2] = 1.0;
	    data[idx+3] = 0.0;
	  }
      }

  gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, width, height, GL_RGBA, GL_FLOAT, data);
  // glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_FLOAT, data);
  
  // 
  // make the normal map
  //
  glBindTexture(GL_TEXTURE_2D, point_sprite_textures[1]);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  for (unsigned int y=0; y<height; y++)
    for (unsigned int x=0; x<width; x++)
      {
	int idx = y*width*4 + x*4;

	// convert to sphere coordinates with radius half of width/height
	double r = (double)radius;
	double xr = x - (double)width/2.0;
	double yr = y - (double)height/2.0;

	// if coordinate is on or inside the circle, keep going
	double mag = sqrt(xr*xr + yr*yr);
	if (mag <= r)
	  {
	    // calculate what z would be for this x, y, radius combination
	    double z = sqrt( r*r - xr*xr - yr*yr );

	    // convert to positive half space since some x and y will
	    // be negative (bad for RGBA values)
	    xr /= r;
	    xr = (xr + 1.0) / 2.0;

	    yr /= r;
	    yr = (yr + 1.0) / 2.0;

	    z /= r;
	    z = (z + 1.0) / 2.0;

	    // store in image space
	    data[idx] = xr;
	    data[idx+1] = yr;
	    data[idx+2] = z;
	    data[idx+3] = 1.0;
	  }
	else if ((mag >= r+(radius*0.20)) && (mag <= r+(radius*0.40)))
	  {
	    // store transparent white in image space
	    data[idx] = 0.0; data[idx+1] = 0.0; data[idx+2] = 1.0;
	    data[idx+3] = 1.0;
	  }
	else 
	  {
	    // store transparent white in image space
	    data[idx] = data[idx+1] = data[idx+2] = 1.0;
	    data[idx+3] = 0.0;
	  }
      }
	
  // glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_FLOAT, data);
  gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, width, height, GL_RGBA, GL_FLOAT, data);

  glDisable(GL_TEXTURE_2D);

  delete [] data;
}


void DisplayControl::syncTreadportData() {
#if !WIN32
  // Note, we need to sync the treadport information before
  // we can sync it to the other computers.

  // We should only sync if we are in the treadport viewing mode
  // and this is the broadcasting computer.

  if(viewingMode == TREADPORT & network_mode == network.BROADCAST) {
    
    // Sync data.
    treadport->sync();
    
    // Set the retrieved values. Note, if a conversion matrix is
    // needed to switch coordinate systems from the treadport
    // coordinate system to the application coordinate system
    // this is where that could be done (or within the 
    // TreadportManager it self).

    std::vector<float> new_eye_pos = treadport->getEyePosition();
    std::vector<float> new_eye_gaze = treadport->getEyeGaze();
    std::vector<float> new_eye = treadport->getEyeOffset();

    eye_pos[0] = new_eye_pos[0];
    eye_pos[1] = new_eye_pos[1];
    eye_pos[2] = new_eye_pos[2];

    eye_gaze[0] = new_eye_gaze[0];
    eye_gaze[1] = new_eye_gaze[1];
    eye_gaze[2] = new_eye_gaze[2];

    eye[0] = new_eye[0];
    eye[1] = new_eye[1];
    eye[2] = new_eye[2];

  }
#endif
}

void DisplayControl::syncDataOverNetwork() {
#if !WIN32
   NetworkManager::NetworkSyncData data;

   memcpy(data.eye_pos, eye_pos, sizeof(float) * 3);
   memcpy(data.eye_gaze, eye_gaze, sizeof(float) * 3);

   float tmp[3];
   tmp[0] = eye[0];
   tmp[1] = eye[1];
   tmp[2] = eye[2];

   memcpy(data.eye_offset, tmp, sizeof(float) * 3);

   memcpy(&data.inPauseMode, &inPauseMode, sizeof(bool));

   // The following data has not yet been implemented in this version of gpuPlume
   // and thus will not be synced.

   // memcpy(data.eye_up, eye_up, sizeof(float) * 3);
   // memcpy(&data.wim_scale, &wim_scale_factor, sizeof(float));
   // memcpy(&data.gesture_id, &gesture_id, sizeof(int));
   // memcpy(data.hand_pos, current_hand_pos, sizeof(float) * 3);
   // memcpy(data.hand_up, current_hand_up, sizeof(float) * 3);

   network.sync(data);

   memcpy(eye_pos, data.eye_pos, sizeof(float) * 3);
   memcpy(eye_gaze, data.eye_gaze, sizeof(float) * 3);
   memcpy(tmp, data.eye_offset, sizeof(float) * 3);
   memcpy(&inPauseMode, &data.inPauseMode, sizeof(bool));

   eye[0] = tmp[0];
   eye[1] = tmp[1];
   eye[2] = tmp[2];

   // The following data has not yet been implemented in this version of gpuPlume
   // and thus will not be synced.

   // memcpy(eye_up, data.eye_up, sizeof(float) * 3);
   // memcpy(&wim_scale_factor, &data.wim_scale, sizeof(float));
   // memcpy(&gesture_id, &data.gesture_id, sizeof(int));
   // memcpy(current_hand_pos, data.hand_pos, sizeof(float) * 3);
   // memcpy(current_hand_up, data.hand_up, sizeof(float) * 3);
   // memcpy(data.hand_matrix, hand_matrix, sizeof(float) * 16);
#endif
}

void DisplayControl::initializeView() {

  if(osgPlume) {
    // Don't do anything.
    return;
  }

  // Set up the viewing frustum according to what view mode is set.
  if(viewingMode == STANDARD) {

    // Modify the view frustum. Note that the legacy code for OpenScene Graph is
    // in which ever model you are running (i.e. MultipleBuildingsModel).
    glViewport(0, 0, glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0, glutGet(GLUT_WINDOW_WIDTH)/float(glutGet(GLUT_WINDOW_HEIGHT)), 1.0, 250.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glClearColor(util->bcolor[0], util->bcolor[1], util->bcolor[2], 1.0);
      gluLookAt( eye_pos[0], eye_pos[1], eye_pos[2],
                 eye_gaze[0]+eye_pos[0], eye_gaze[1]+eye_pos[1], eye_gaze[2]+eye_pos[2],
                 0, 0, 1 );
  }
  else if (viewingMode == ORTHOGRAPHIC_TOP)
    {
#if 0
    // Modify the view frustum. Note that the legacy code for OpenScene Graph is
    // in which ever model you are running (i.e. MultipleBuildingsModel).
    glViewport(0, 0, glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    // compute a perspective view that encompasses the domain
    float d = sqrt(nx*nx + ny*ny);
    float ht = tan(30.0 * M_PI/180.0) * d;
    
    Vector3 pt2(nx, ny, ht);
    Vector3 mpt1(nx, 0, pt2[2]/2.0);
    Vector3 mpt2(0, ny, pt2[2]/2.0);
    Vector3 mdpt3 = (mpt2 + mpt1)/2.0;

    float b1 = (mdpt3 - pt2).norm();
    float hyp = (mpt1 - pt2).norm();
    float horizTheta = asin( b1/hyp ) * 2.0 * 180.0/M_PI;
    

    float aspectRatio = glutGet(GLUT_WINDOW_WIDTH)/float(glutGet(GLUT_WINDOW_HEIGHT));
    gluPerspective(horizTheta/aspectRatio, aspectRatio, 0.5, 500.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    gluLookAt( pt2[0], pt2[1], pt2[2],
	       0.0, 0.0, 0.0,
	       0, 0, 1 );

    glClearColor(util->bcolor[0], util->bcolor[1], util->bcolor[2], 1.0);
#endif

      /// to be completed orthographic views
      
      // Modify the view frustum. Note that the legacy code for OpenScene Graph is
      // in which ever model you are running (i.e. MultipleBuildingsModel).
    glViewport(0, 0, glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1.0, nx+1.0, -1.0, ny+1.0, -nz, nz*2.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glClearColor(util->bcolor[0], util->bcolor[1], util->bcolor[2], 1.0);

    // glTranslatef(nx/2.0, ny/2.0, 0.0);
    // glRotatef(180.0, 0.0, 1.0, 0.0);
    // glRotatef(-90.0, 0.0, 0.0, 1.0);
    glTranslatef(eye_pos[0], 0, eye_pos[1]);

    } 
  else if(viewingMode == VR) {
    // TODO: Handle the left and right eye; at the moment we are just doing
    // the same as standard.

    // Modify the view frustum. Note that the legacy code for OpenScene Graph is
    // in which ever model you are running (i.e. MultipleBuildingsModel).
    glViewport(0, 0, glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0, glutGet(GLUT_WINDOW_WIDTH)/float(glutGet(GLUT_WINDOW_HEIGHT)), 1.0, 250.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glClearColor(util->bcolor[0], util->bcolor[1], util->bcolor[2], 1.0);
      gluLookAt( eye_pos[0], eye_pos[1], eye_pos[2],
                 eye_gaze[0]+eye_pos[0], eye_gaze[1]+eye_pos[1], eye_gaze[2]+eye_pos[2],
                 0, 0, 1 );
  } else if(viewingMode == TREADPORT) {
#if !WIN32
    // Initialize 1st treadport view
    calcTreadportFrustum(util->treadport_view);

    // Now push so we can initialize the gpuPlume view
    glPushMatrix();
    // glLoadIdentity();
    // glRotatef(90.0f, 1.0f, 0.0f, 0.0f);
    // glRotatef(180.0f, 0.0f, 0.0f, 0.0f);
#endif
  }
}

void DisplayControl::deinitializeView() {
  if(viewingMode == TREADPORT) {
#if !WIN32
    glPopMatrix();
#endif
  }
}

void DisplayControl::calcTreadportFrustum(char view) {
#if !WIN32
  // TODO: Change to use frustum call parameters from frustumCalc app w/ 't'
  //       option.

  // Parameters for use with original treadport setup.
  float W = 3.048; // 10ft wide "effective" screens
  float H = W/2.0;  // half of effective screen width
  float trueH = (W - 0.6096) / 2.0;  // get back to half of an 8 foot screen - the true width
  float f = 0.0; // no f in this equation
  float theta = 30.0 * M_PI/180.0;  // convert 30.0 degrees to radians
  float overscan = (1.0 - 2.4384/W) / 2.0;
  float screen_height = 2.4384; // 8.0ft
  float offset_from_bottom = 0.0 * 0.0254; // inches converted to meters
  float zbot = offset_from_bottom;
  float ztop = offset_from_bottom + screen_height; // relative height of screen top (to feet) not including frame
  float D;
  float a, b, c;
  a = trueH * cos(theta);
  c = trueH * sin(theta);
  b = (trueH + c) * tan(theta);
  D = a + b;

  Vector4 coord, transformed( 1.0, 1.0, 1.0, 1.0 );

  Matrix4x4 rotation;
  rotation.set2Identity();

  float lrll[3], lrur[3], lrll_mag, lrur_mag;

  screen.lr.resize(3);
  screen.ll.resize(3);
  screen.ul.resize(3);
  screen.ur.resize(3);

  // Set the initial position of the eye if we are using a static frusstum.
  // This is currently disabled because it is set when we sync
  // the data from the treadport.
  if(static_treadport_frustum) {
    eye[0] = 0.0;
    eye[1] = 0.0;
    eye[2] = screen_height / 2.0;
  }
  
  if(view == 'l') { // left screen.

    //
    // Left Screen
    //

    rotation.setRotationAboutZ( 60.0 * M_PI / 180.0 ); // 60.0 degrees in radians

    // transform ll

    coord.set( -H + f, D, zbot, 1.0 );
    transformed =  rotation * coord;
    screen.ll[0] = transformed[0];       screen.ll[1] = transformed[1];     screen.ll[2] = transformed[2];

    coord.set( -H + f, D, ztop, 1.0 );
    transformed =  rotation * coord;
    screen.ul[0] = transformed[0];       screen.ul[1] = transformed[1];     screen.ul[2] = transformed[2];

    coord.set( H - f, D, zbot, 1.0 );
    transformed =  rotation * coord;
    screen.lr[0] = transformed[0];       screen.lr[1] = transformed[1];     screen.lr[2] = transformed[2];

    coord.set( H - f, D, ztop, 1.0 );
    transformed =  rotation * coord;
    screen.ur[0] = transformed[0];       screen.ur[1] = transformed[1];     screen.ur[2] = transformed[2];

    screen.phi = theta - (90.0 * M_PI/180.0);

    screen.far1.resize(3);
    screen.far2.resize(3);
    screen.far3.resize(3);
    screen.far4.resize(3);
    screen.proj.resize(3);
    screen.normal.resize(3);

    lrll[0] = screen.ll[0] - screen.lr[0];
    lrll[1] = screen.ll[1] - screen.lr[1];
    lrll[2] = screen.ll[2] - screen.lr[2];
    lrll_mag = sqrt( lrll[0]*lrll[0] + lrll[1]*lrll[1] + lrll[2]*lrll[2] );
    lrll[0] /= lrll_mag;
    lrll[1] /= lrll_mag;
    lrll[2] /= lrll_mag;

    lrur[0] = screen.ur[0] - screen.lr[0];
    lrur[1] = screen.ur[1] - screen.lr[1];
    lrur[2] = screen.ur[2] - screen.lr[2];
    lrur_mag = sqrt( lrur[0]*lrur[0] + lrur[1]*lrur[1] + lrur[2]*lrur[2] );
    lrur[0] /= lrur_mag;
    lrur[1] /= lrur_mag;
    lrur[2] /= lrur_mag;

    screen.asp_ratio = lrur_mag / lrll_mag;

    screen.normal[0] = lrur[1] * lrll[2] - lrur[2] * lrll[1];
    screen.normal[1] = lrur[2] * lrll[0] - lrur[0] * lrll[2];
    screen.normal[2] = lrur[0] * lrll[1] - lrur[1] * lrll[0];

    screen.center.resize(3);
    screen.center[0] = (screen.lr[0] + screen.ll[0] + screen.ul[0] + screen.ur[0]) / 4.0;
    screen.center[1] = (screen.lr[1] + screen.ll[1] + screen.ul[1] + screen.ur[1]) / 4.0;
    screen.center[2] = (screen.lr[2] + screen.ll[2] + screen.ul[2] + screen.ur[2]) / 4.0;

    /*
    cout << "Left Screen Config:" << endl;
    cout << "\tll=[" << screen.ll[0] << ", " << screen.ll[1] << ", " << screen.ll[2] << "]" << endl;
    cout << "\tul=[" << screen.ul[0] << ", " << screen.ul[1] << ", " << screen.ul[2] << "]" << endl;
    cout << "\tlr=[" << screen.lr[0] << ", " << screen.lr[1] << ", " << screen.lr[2] << "]" << endl;
    cout << "\tur=[" << screen.ur[0] << ", " << screen.ur[1] << ", " << screen.ur[2] << "]" << endl;
    cout << "\tnormal=" << screen.normal[0] << ", " << screen.normal[1] << ", " << screen.normal[2] << endl;
    */

  } else if(view == 'c') { // center screen.

    //
    // Center Screen
    //
    screen.ll[0] = -H + f;       screen.ll[1] = D;     screen.ll[2] = zbot;
    screen.ul[0] = -H + f;       screen.ul[1] = D;     screen.ul[2] = ztop;
    screen.lr[0] = H - f;        screen.lr[1] = D;     screen.lr[2] = zbot;
    screen.ur[0] = H - f;        screen.ur[1] = D;     screen.ur[2] = ztop;

    screen.phi = 0.0; // was 0.0

    screen.far1.resize(3);
    screen.far2.resize(3);
    screen.far3.resize(3);
    screen.far4.resize(3);
    screen.proj.resize(3);
    screen.normal.resize(3);

    lrll[0] = screen.ll[0] - screen.lr[0];
    lrll[1] = screen.ll[1] - screen.lr[1];
    lrll[2] = screen.ll[2] - screen.lr[2];
    lrll_mag = sqrt( lrll[0]*lrll[0] + lrll[1]*lrll[1] + lrll[2]*lrll[2] );
    lrll[0] /= lrll_mag;
    lrll[1] /= lrll_mag;
    lrll[2] /= lrll_mag;

    lrur[0] = screen.ur[0] - screen.lr[0];
    lrur[1] = screen.ur[1] - screen.lr[1];
    lrur[2] = screen.ur[2] - screen.lr[2];
    lrur_mag = sqrt( lrur[0]*lrur[0] + lrur[1]*lrur[1] + lrur[2]*lrur[2] );
    lrur[0] /= lrur_mag;
    lrur[1] /= lrur_mag;
    lrur[2] /= lrur_mag;

    screen.asp_ratio = lrur_mag / lrll_mag;

    screen.normal[0] = lrur[1] * lrll[2] - lrur[2] * lrll[1];
    screen.normal[1] = lrur[2] * lrll[0] - lrur[0] * lrll[2];
    screen.normal[2] = lrur[0] * lrll[1] - lrur[1] * lrll[0];

    screen.center.resize(3);
    screen.center[0] = (screen.lr[0] + screen.ll[0] + screen.ul[0] + screen.ur[0]) / 4.0;
    screen.center[1] = (screen.lr[1] + screen.ll[1] + screen.ul[1] + screen.ur[1]) / 4.0;
    screen.center[2] = (screen.lr[2] + screen.ll[2] + screen.ul[2] + screen.ur[2]) / 4.0;

    /*
    cout << "Center Screen Config:" << endl;
    cout << "\tll=[" << screen.ll[0] << ", " << screen.ll[1] << ", " << screen.ll[2] << "]" << endl;
    cout << "\tul=[" << screen.ul[0] << ", " << screen.ul[1] << ", " << screen.ul[2] << "]" << endl;
    cout << "\tlr=[" << screen.lr[0] << ", " << screen.lr[1] << ", " << screen.lr[2] << "]" << endl;
    cout << "\tur=[" << screen.ur[0] << ", " << screen.ur[1] << ", " << screen.ur[2] << "]" << endl;
    cout << "\tnormal=" << screen.normal[0] << ", " << screen.normal[1] << ", " << screen.normal[2] << endl;
    */

  } else if(view == 'r') { // right screen.

    //
    // Right Screen
    //
    rotation.setRotationAboutZ( -60.0 * M_PI / 180.0 );

    // transform ll

    coord.set( -H + f, D, zbot, 1.0 );
    transformed =  rotation * coord;
    screen.ll[0] = transformed[0];       screen.ll[1] = transformed[1];     screen.ll[2] = transformed[2];

    coord.set( -H + f, D, ztop, 1.0 );
    transformed =  rotation * coord;
    screen.ul[0] = transformed[0];       screen.ul[1] = transformed[1];     screen.ul[2] = transformed[2];

    coord.set( H - f, D, zbot, 1.0 );
    transformed =  rotation * coord;
    screen.lr[0] = transformed[0];       screen.lr[1] = transformed[1];     screen.lr[2] = transformed[2];

    coord.set( H - f, D, ztop, 1.0 );
    transformed =  rotation * coord;
    screen.ur[0] = transformed[0];       screen.ur[1] = transformed[1];     screen.ur[2] = transformed[2];

    screen.phi = (90.0 * M_PI/180.0) - theta;

    screen.far1.resize(3);
    screen.far2.resize(3);
    screen.far3.resize(3);
    screen.far4.resize(3);
    screen.proj.resize(3);
    screen.normal.resize(3);

    float lrll[3], lrur[3], lrll_mag, lrur_mag;
    lrll[0] = screen.ll[0] - screen.lr[0];
    lrll[1] = screen.ll[1] - screen.lr[1];
    lrll[2] = screen.ll[2] - screen.lr[2];
    lrll_mag = sqrt( lrll[0]*lrll[0] + lrll[1]*lrll[1] + lrll[2]*lrll[2] );
    lrll[0] /= lrll_mag;
    lrll[1] /= lrll_mag;
    lrll[2] /= lrll_mag;

    lrur[0] = screen.ur[0] - screen.lr[0];
    lrur[1] = screen.ur[1] - screen.lr[1];
    lrur[2] = screen.ur[2] - screen.lr[2];
    lrur_mag = sqrt( lrur[0]*lrur[0] + lrur[1]*lrur[1] + lrur[2]*lrur[2] );
    lrur[0] /= lrur_mag;
    lrur[1] /= lrur_mag;
    lrur[2] /= lrur_mag;

    screen.asp_ratio = lrur_mag / lrll_mag;

    screen.normal[0] = lrur[1] * lrll[2] - lrur[2] * lrll[1];
    screen.normal[1] = lrur[2] * lrll[0] - lrur[0] * lrll[2];
    screen.normal[2] = lrur[0] * lrll[1] - lrur[1] * lrll[0];

    screen.center.resize(3);
    screen.center[0] = (screen.lr[0] + screen.ll[0] + screen.ul[0] + screen.ur[0]) / 4.0;
    screen.center[1] = (screen.lr[1] + screen.ll[1] + screen.ul[1] + screen.ur[1]) / 4.0;
    screen.center[2] = (screen.lr[2] + screen.ll[2] + screen.ul[2] + screen.ur[2]) / 4.0;

    /*
    cout << "Right Screen Config:" << endl;
    cout << "\tll=[" << screen.ll[0] << ", " << screen.ll[1] << ", " << screen.ll[2] << "]" << endl;
    cout << "\tul=[" << screen.ul[0] << ", " << screen.ul[1] << ", " << screen.ul[2] << "]" << endl;
    cout << "\tlr=[" << screen.lr[0] << ", " << screen.lr[1] << ", " << screen.lr[2] << "]" << endl;
    cout << "\tur=[" << screen.ur[0] << ", " << screen.ur[1] << ", " << screen.ur[2] << "]" << endl;
    cout << "\tnormal=" << screen.normal[0] << ", " << screen.normal[1] << ", " << screen.normal[2] << endl;
    */

  }

  GLfloat left, right, bottom, top, near, far;
  std::vector<float> eyeproj;

  Matrix4x4 translation;
  Vector4 Cs, Ct;
  rotation.set2Identity();
  translation.set2Identity();

  rotation.setRotationAboutZ(screen.phi);
  translation.setTranslationV( Vector3( -eye[0], -eye[1], -eye[2] ) );

  // transform ll
  Cs =  rotation * translation * Vector4( screen.ll[0], screen.ll[1], screen.ll[2], 1.0 );
  left = Cs[0];
  bottom = Cs[2];

  // tranform ur
  Cs = rotation * translation * Vector4( screen.ur[0], screen.ur[1], screen.ur[2], 1.0 );
  right = Cs[0];
  top = Cs[2];
  near = Cs[1];
  far = near + 250.0; // was + 20.0

  glViewport(0, 0, glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glFrustum(left, right, bottom, top, near, far);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glClearColor(util->bcolor[0], util->bcolor[1], util->bcolor[2], 1.0);

  float dotprod;
  float center_eye[3], tmll[3], pt[3], mag;

  //
  // Begin transforming the screens and the eye into the world space.
  //
  // Note, that when calculating the frustum, the x & y coordinates were swapped
  // but z is still up. TODO: Fix this, so everything is in the same coordinate
  // system.

  // Adjust the screen center & normals so we are facing the right direction.
  rotation.set2Identity();
  rotation.setRotationAboutZ(M_PI); // 180 degrees.

  Vector4 adjusted_screen_center = rotation * Vector4(screen.center[0], screen.center[1], screen.center[2], 1.0);;
  Vector4 adjusted_screen_normal = rotation * Vector4(screen.normal[0], screen.normal[1], screen.normal[2], 1.0);;

  screen.center[0] = adjusted_screen_center[0];
  screen.center[1] = adjusted_screen_center[1];
  screen.center[2] = adjusted_screen_center[2];

  screen.normal[0] = adjusted_screen_normal[0];
  screen.normal[1] = adjusted_screen_normal[1];
  screen.normal[2] = adjusted_screen_normal[2];


  // Transform the eye (on the treadport, or "local" space) into the gpuPlume
  // world (or "global" or "world" space).
  eye[0] += eye_pos[1]; // <-- --> (side to side)
  eye[2] += eye_pos[2]; // ^ (UP)
  eye[1] += eye_pos[0]; // (forward and backward)

  // Transform center of screen.
  screen.center[0] += eye_pos[1];
  screen.center[2] += eye_pos[2];
  screen.center[1] += eye_pos[0];

  //
  // Done transforming the screens and the eye into the world space.
  //

  //
  // Begin calculate projection of eye point onto each screen plane
  //

  // this can be done by dotting the screen normal with the center_eye
  // vector and then subtracting this distance along the normal vector
  // from the eye point.

  center_eye[0] = eye[0] - screen.center[0];
  center_eye[1] = eye[1] - screen.center[1];
  center_eye[2] = eye[2] - screen.center[2];

  dotprod = screen.normal[0] * center_eye[0] +
            screen.normal[1] * center_eye[1] +
            screen.normal[2] * center_eye[2];

  screen.proj[0] = eye[0] - screen.normal[0] * dotprod;
  screen.proj[1] = eye[1] - screen.normal[1] * dotprod;
  screen.proj[2] = eye[2] - screen.normal[2] * dotprod;

  // Calculate the rotation matrix that needs to be applied to the screen
  // projection to rotate the projection of the eye on the screen into the
  // opengl "world" current gaze direction.

  // First create and normalize the "starting" or "default" direction vector
  // and the current "world" or "global" gaze direction vector.

  Vector3 default_gaze = Vector3(0.0, 0.0, -1.0); // TODO: Remove hard coding (this is the default "gaze").
  Vector3 world_gaze = Vector3(eye_gaze[1], eye_gaze[2], eye_gaze[0]);

  // Now normalize.
  default_gaze.normalize();
  world_gaze.normalize();

  // Then take the cross-product of those two vectors to get the axis on which
  // we need to rotate on (this was Josh's "duh" moment, *sigh*). Since OBVIOUSLY
  // you want to rotate along the axis that is orthogonal to the two vectors.
  Vector3 rotation_axis = default_gaze.cross(world_gaze);

  float rotationMag = rotation_axis.norm();

  // Get the angle of rotation.
  float angle_of_rotation = acosf(default_gaze.dot(world_gaze) / (default_gaze.norm() * world_gaze.norm()));

  // Normalize.
  rotation_axis.normalize();

  float cos_theta = cos(angle_of_rotation);
  float sin_theta = sin(angle_of_rotation);

  float x = rotation_axis[0];
  float y = rotation_axis[1];
  float z = rotation_axis[2];

  Matrix4x4 rotation_matrix = Matrix4x4(
      (double)(cos_theta + ( 1 - cos_theta) * (x * x)), (double)((1 - cos_theta) * x * y - z * sin_theta), (double)((1 - cos_theta) * x * z + y * sin_theta), 0,
      (double)((1 - cos_theta) * x * y + z * sin_theta), (double)(cos_theta + (1 - cos_theta) * (y * y)), (double)((1 - cos_theta) * y * z - x * sin_theta), 0,
      (double)((1 - cos_theta) * x * z - y * sin_theta), (double)((1 - cos_theta) * y * z + x * sin_theta), (double)(cos_theta + (1 - cos_theta) * (z * z)), 0,
      0, 0, 0, 1
  );

  // We are now in TREADPORT unit space.
  Vector4 orig_point = Vector4(screen.proj[0] - eye[0], screen.proj[1] - eye[1], screen.proj[2] - eye[2], 1.0);
  rotation.set2Identity();
  rotation.setRotationAboutZ(screen.phi); // Note this isn't negative because we already rotated 180 degress (so we are backwards in treadport space).
  Vector4 rotated_look_at = rotation * orig_point;

  // We are now in NORMAL unit space.
  Vector4 aligned_look_at = rotation_matrix * Vector4(rotated_look_at[0], rotated_look_at[2], rotated_look_at[1], 1.0);

  Vector3 standard_up = Vector3(0.0, 1.0, 0.0); // Again this is hard coded, TODO: Fix the hard coded "default" up axis.

  Vector4 aligned_axis_temp = rotation_matrix * Vector4(standard_up[0], standard_up[1], standard_up[2], 1.0);
  Vector3 aligned_axis = Vector3(aligned_axis_temp[0], aligned_axis_temp[1], aligned_axis_temp[2]);
  aligned_axis.normalize();

  // cout << "Aligned Axis: " << aligned_axis[0] << " " << aligned_axis[1] << " " << aligned_axis[2] << "\n" << endl;

  // Rotate the screen "back"

  // Note, we are reusing values, so be careful...

  cos_theta = cos(-screen.phi);// Note this is negative because we already rotated 180 degress (so we are backwards in treadport space).
  sin_theta = sin(-screen.phi);

  x = aligned_axis[0];
  y = aligned_axis[1];
  z = aligned_axis[2];

  rotation_matrix = Matrix4x4(
      (double)(cos_theta + ( 1 - cos_theta) * (x * x)), (double)((1 - cos_theta) * x * y - z * sin_theta), (double)((1 - cos_theta) * x * z + y * sin_theta), 0,
      (double)((1 - cos_theta) * x * y + z * sin_theta), (double)(cos_theta + (1 - cos_theta) * (y * y)), (double)((1 - cos_theta) * y * z - x * sin_theta), 0,
      (double)((1 - cos_theta) * x * z - y * sin_theta), (double)((1 - cos_theta) * y * z + x * sin_theta), (double)(cos_theta + (1 - cos_theta) * (z * z)), 0,
      0, 0, 0, 1
  );

  Vector4 realigned_look_at = rotation_matrix * aligned_look_at;

  // Done Rotating the screen "back"

  // We are in TREADPORT unit space.

  //
  // Done calculating the projection of the eye point onto the screen frame.
  //

  // cout << "-----------------------------------------\n" << endl;


  // This one works (rotation does not work but looking up and down does).
  gluLookAt(eye[1], eye[0], eye[2],
      realigned_look_at[2] + eye[1], realigned_look_at[0] + eye[0], realigned_look_at[1] + eye[2],
      aligned_axis[2], aligned_axis[0], aligned_axis[1]);


  /* This one also works (sides do not work)
  gluLookAt(eye[1], eye[0], eye[2],
      realigned_look_at[2] + eye[1], realigned_look_at[0] + eye[0], realigned_look_at[1] + eye[2],
      0, 0, 1);
  */

  /* This one works (can not look around at all)
  gluLookAt(eye[1], eye[0], eye[2],
            screen.proj[1], screen.proj[0], screen.proj[2],
            0, 0, 1);
  */
#endif
}

void DisplayControl::blankSides()
{
  
  int lines;
  char* p;
  GLint* vp = new GLint[4];
  glGetIntegerv(GL_VIEWPORT,vp);

  // glDisable(GL_LIGHTING);
  glDisable(texType);
  glDisable(GL_DEPTH_TEST);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0, vp[2],
    0, vp[3], -1, 1);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  
  // glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT)

  // Blank out left 10%
  glColor3ub(0, 0, 0);
  glBegin(GL_POLYGON);
    glVertex2f(0, glutGet(GLUT_WINDOW_HEIGHT));
    glVertex2f(glutGet(GLUT_WINDOW_WIDTH) * 0.10, glutGet(GLUT_WINDOW_HEIGHT));
    glVertex2f(glutGet(GLUT_WINDOW_WIDTH) * 0.10, 0);
    glVertex2f(0, 0);
  glEnd();

  // Blank out right 10%
  glColor3ub(0, 0, 0);
  glBegin(GL_POLYGON);
    glVertex2f(glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
    glVertex2f(glutGet(GLUT_WINDOW_WIDTH) - glutGet(GLUT_WINDOW_WIDTH) * 0.10, glutGet(GLUT_WINDOW_HEIGHT));
    glVertex2f(glutGet(GLUT_WINDOW_WIDTH) - glutGet(GLUT_WINDOW_WIDTH) * 0.10, 0);
    glVertex2f(glutGet(GLUT_WINDOW_WIDTH), 0);
    glEnd();

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glEnable(GL_DEPTH_TEST);
  glEnable(texType);
  // glEnable(GL_LIGHTING);
  
}

bool DisplayControl::getInPauseMode() {
  return inPauseMode;
}

void DisplayControl::setInPauseMode(bool newInPauseMode) {
  inPauseMode = newInPauseMode;
}

void DisplayControl::drawInShadowData() {

  for(int i = 0; i < nz; i++) {
    for(int j = 0; j < ny; j++) {
      for(int k = 0; k < nx; k++) {
				int index = i*ny*nx*4 + j*nx*4 + k*4;

        /* Center
        if(inShadowData[index] < 1) {
					glPushMatrix();
					glTranslatef(k * util->dx + 0.5 * util->dx, j * util->dy + 0.5 * util->dy, i * util->dz + 0.5 * util->dz);
					glColor4f(inShadowData2[index], inShadowData2[index], inShadowData2[index], 1);
					glutSolidCube(0.08);
					glPopMatrix();
				}*/

        
				// Bottom
				if(inShadowData[index] < 1) {
					glPushMatrix();
					glTranslatef(k * util->dx + 0.5 * util->dx, j * util->dy + 0.5 * util->dy, i * util->dz);
					glColor4f(inShadowData[index], inShadowData[index], inShadowData[index], 1);
					// std::cout << inShadowData[i][j][0] << " " << inShadowData[i][j][1] << " " << inShadowData[i][j][2] << "  ";
					glutSolidCube(0.05);
					glPopMatrix();
				}
				
				// North
				if(inShadowData[index + 1] < 1) {
					glPushMatrix();
					glTranslatef(k * util->dx, j * util->dy + 0.5 * util->dy, i * util->dz + 0.5 * util->dz);
					glColor4f(inShadowData[index + 1], inShadowData[index + 1], inShadowData[index + 1], 1);
					// std::cout << inShadowData[i][j][0] << " " << inShadowData[i][j][1] << " " << inShadowData[i][j][2] << "  ";
					glutSolidCube(0.05);
					glPopMatrix();
				}
				
				// West
				if(inShadowData[index + 2] < 1) {
					glPushMatrix();
					glTranslatef(k * util->dx + 0.5 * util->dx, j * util->dy, i * util->dz + 0.5 * util->dz);
					glColor4f(inShadowData[index + 2], inShadowData[index + 2], inShadowData[index + 2], 1);
					// std::cout << inShadowData[i][j][0] << " " << inShadowData[i][j][1] << " " << inShadowData[i][j][2] << "  ";
					glutSolidCube(0.05);
					glPopMatrix();
				}

        // Top
				if(inShadowData2[index] < 1) {
					glPushMatrix();
					glTranslatef(k * util->dx + 0.5 * util->dx, j * util->dy + 0.5 * util->dy, i * util->dz + util->dz);
					glColor4f(inShadowData2[index], inShadowData2[index], inShadowData2[index], 1);
					// std::cout << inShadowData[i][j][0] << " " << inShadowData[i][j][1] << " " << inShadowData[i][j][2] << "  ";
					glutSolidCube(0.05);
					glPopMatrix();
				}
				
				// South
				if(inShadowData2[index + 1] < 1) {
					glPushMatrix();
					glTranslatef(k * util->dx + util->dx, j * util->dy + 0.5 * util->dy, i * util->dz + 0.5 * util->dz);
					glColor4f(inShadowData2[index + 1], inShadowData2[index + 1], inShadowData2[index + 1], 1);
					// std::cout << inShadowData[i][j][0] << " " << inShadowData[i][j][1] << " " << inShadowData[i][j][2] << "  ";
					glutSolidCube(0.05);
					glPopMatrix();
				}
				
				// East
				if(inShadowData2[index + 2] < 1) {
					glPushMatrix();
					glTranslatef(k * util->dx + 0.5 * util->dx, j * util->dy + util->dy, i * util->dz + 0.5 * util->dz);
					glColor4f(inShadowData2[index + 2], inShadowData2[index + 2], inShadowData2[index + 2], 1);
					// std::cout << inShadowData[i][j][0] << " " << inShadowData[i][j][1] << " " << inShadowData[i][j][2] << "  ";
					glutSolidCube(0.05);
					glPopMatrix();
				}

      }
    }
    // std::cout << std::endl;
  }
  
}


