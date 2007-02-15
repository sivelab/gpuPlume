#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <list>

#include <GL/glew.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "particleControl.h"
#include "displayControl.h"
#include "framebufferObject.h"
#include "renderbuffer.h"
#include "GLSL.h"
#include "glErrorUtil.h"

#ifdef WIN32
#include <windows.h>
#include <stdio.h>
#include <conio.h>
#include <tchar.h>

// Rand functions
float randVal() { return (float)(rand()/(float)RAND_MAX); } 
#else
float randVal() { return drand48(); }
#endif

//These values hold the domain of the 3D area
int nx;
int ny; //This is up for our orientation
int nz;
float time_step = 0.0012; //This is the time step for the movement of particles

std::list<int> indices;

ParticleControl* pc;
DisplayControl* dc;
int numInRow;   // represents the number of layers that are stored in one row in a 2D texture
int visual_layer = -1;    // the layer we're visualizing in the display (if == -1, do not display)

bool frame_rate = true;
bool dump_contents = false;
bool emit = false;

void init(void);
void initFBO(void);
void idle();
void reshape(int w, int h);
void display(void);
void keyboard_cb(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);

bool show_particle_visuals = true;

int winid;
float eye_z = 5.0;

int winwidth = 512, winheight = 512;

//
// The values here determine the number of particles
//
int twidth = 128, theight = 128;

static bool rotate_sphere = false;
static bool rotate_object = false;
static bool translate_view = false;
static GLfloat azimuth = 0.0;
static GLfloat elevation = 0.0;
static float eye_pos[3];
static float eye_gaze[3];

GLuint texid[6];
GLint  uniform_postex, uniform_wind, uniform_timeStep;
GLSLObject init_shader, pass1_shader, render_shader, emit_shader;
GLuint vertex_buffer;

GLenum texType = GL_TEXTURE_RECTANGLE_ARB;
GLenum int_format = GL_RGBA32F_ARB;
GLenum int_format_init = GL_RGBA;

int main(int argc, char** argv)
{
  if (argc == 2)
    {
      twidth = atoi(argv[1]);
      theight = twidth;
    }

#ifdef WIN32
  TCHAR buffer[MAX_PATH];
  DWORD dwRet;

  // dwRet = GetCurrentDirectory(MAX_PATH, buffer);

  // Set the current working directory back a level so shader access is uniform across platforms
  if (!SetCurrentDirectory(_T("..")))
  {
	  std::cerr << "SetCurrentDirectory failed (" << GetLastError() << ")" << std::endl;
  }
  else 
  {
	   dwRet = GetCurrentDirectory(MAX_PATH, buffer);
	   std::cout << "Current directory set to " << buffer << std::endl;
	   system("pause");
  }

#endif

#ifndef WIN32
  srand48( time(0) % getpid() );
#else
  srand(2);
#endif

  glutInitDisplayMode( GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE );
  glutInit(&argc, argv);
  glutInitWindowSize(winwidth, winheight);
  winid = glutCreateWindow("gpuPLUME");
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutIdleFunc(idle);
  glutKeyboardFunc(keyboard_cb);
  glutMotionFunc(motion);
  glutMouseFunc(mouse);

  GLenum err = glewInit();
  if (GLEW_OK != err) 
    {
      std::cout << "Error: " << glewGetErrorString(err) << std::endl;
    }

  if (GL_ARB_vertex_buffer_object) 
    {
      std::cout << "GL_ARB_vertex_buffer_object available!" << std::endl;
    }
  else 
    {
      std::cout << "GL_ARB_vertex_buffer_object is NOT available!  Exiting!" << std::endl;
      exit(-1);
    }

  init();

  glutMainLoop();
  return 0;
}

// GLUT reshape function
void reshape(int w, int h)
{
    if (h == 0) h = 1;

    glViewport(0, 0, w, h);

    // GPGPU CONCEPT 3b: One-to-one Pixel to Texel Mapping: An Orthographic
    //                   Projection.
    // This code sets the projection matrix to orthographic with a range of
    // [-1,1] in the X and Y dimensions. This allows a trivial mapping of
    // pixels to texels.
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(-1, 1, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

int sz = 4;
FramebufferObject *fbo;
Renderbuffer *rb;

void initFBO(void){
  // Get the Framebuffer Object ready
  fbo = new FramebufferObject();
  fbo->Bind();
      
  rb = new Renderbuffer();
  rb->Set(GL_DEPTH_COMPONENT24, twidth, theight);
  fbo->AttachRenderBuffer(GL_DEPTH_ATTACHMENT_EXT, rb->GetId() );

  //Attach textures to framebuffer object
  fbo->AttachTexture(GL_COLOR_ATTACHMENT0_EXT, texType, texid[0]);
  fbo->AttachTexture(GL_COLOR_ATTACHMENT1_EXT, texType, texid[1]);

  fbo->IsValid();
  FramebufferObject::Disable();
}

void init(void)
{
  glEnable(GL_DEPTH_TEST);

  pc = new ParticleControl(texType);
  pc->getDomain(&nx,&ny,&nz);

  dc = new DisplayControl(nx, ny, nz, texType);

  for(int i = twidth*theight-1; i >= 0; i--)
    indices.push_back(i);

  eye_pos[0] = 0;
  eye_pos[1] = 0;
  eye_pos[2] = 0;
  
  eye_gaze[0] = 0;
  eye_gaze[1] = 0;
  eye_gaze[2] = -5;

  // load up textures to hold position
  CheckErrorsGL("BEGIN : Creating textures");

  glEnable(texType);
  glGenTextures(6, texid); // create (reference to) a new texture
  GLfloat *data = new GLfloat[ twidth * theight * sz];
  
  for (int j=0; j<theight; j++)
    for (int i=0; i<twidth; i++)
      {
	int idx = j*twidth*sz + i*sz;
	
	//
	// Generate random positions for the particles within the
	// domain.  Currently, the domain is positive.
	//
	// With floating point textures, we have to create the inital
	// values between 0 and 1 and then use an initial shader to
	// transform the normalized coordinates to the correct domain.
	
	data[idx] = randVal();
	data[idx+1] = randVal();
	data[idx+2] = randVal();
	data[idx+3] = randVal();
      }
  pc->createTexture(texid[2], int_format_init, twidth, theight, data);

  // Creates wind field data texture
  pc->initWindTex(texid[3], &numInRow);
  
  for (int j=0; j<theight; j++)
    for (int i=0; i<twidth; i++)
      {
	int idx = j*twidth*sz + i*sz;
	data[idx] = data[idx] * (nx-1) + 100;
	data[idx+1] = data[idx+1] * ny + 100;
	data[idx+2] = data[idx+2] * (nz-1) + 100;
      }
  
  // create the base texture with inital vertex positions
  pc->createTexture(texid[0], int_format, twidth, theight, data);

  // create a second texture to double buffer the vertex positions
  pc->createTexture(texid[1], int_format, twidth, theight, NULL);

  delete [] data;

  CheckErrorsGL("END : Creating textures");

  //
  // set up vertex buffer
  // 
  glGenBuffersARB(1, &vertex_buffer);
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, vertex_buffer);
  glBufferDataARB(GL_ARRAY_BUFFER_ARB, twidth*theight*4*sizeof(GLfloat), 0, GL_STREAM_COPY);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  // Load up the shader programs
  //This shader is used to move the particles
  pass1_shader.addShader("Shaders/plumeAdvect_vp.glsl", GLSLObject::VERTEX_SHADER);
  pass1_shader.addShader("Shaders/plumeAdvect_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  pass1_shader.createProgram();

  // Get location of the sampler uniform
  uniform_postex = pass1_shader.createUniform("pos_texunit");
  uniform_wind = pass1_shader.createUniform("wind_texunit");
  uniform_timeStep = pass1_shader.createUniform("time_step");
  GLint unx = pass1_shader.createUniform("nx");
  GLint uny = pass1_shader.createUniform("ny");
  GLint unz = pass1_shader.createUniform("nz");
  GLint uNumInRow = pass1_shader.createUniform("numInRow");
  pass1_shader.activate();
  glUniform1fARB(unx, nx);
  glUniform1fARB(uny, ny);
  glUniform1fARB(unz, nz);
  glUniform1fARB(uNumInRow, numInRow);
  glUniform1fARB(uniform_timeStep, time_step);
  pass1_shader.deactivate();

  //This shader is used to initialize the particle positions
  init_shader.addShader("Shaders/initialize_vp.glsl", GLSLObject::VERTEX_SHADER);
  init_shader.addShader("Shaders/initialize_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  init_shader.createProgram();

  //This shader is used to make final changes before rendering to the screen
  render_shader.addShader("Shaders/particleVisualize_vp.glsl", GLSLObject::VERTEX_SHADER);
  render_shader.addShader("Shaders/particleVisualize_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  render_shader.createProgram();

  //This shader is used to emmit particles
  emit_shader.addShader("Shaders/emitParticle_vp.glsl", GLSLObject::VERTEX_SHADER);
  emit_shader.addShader("Shaders/emitParticle_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  emit_shader.createProgram();

  // ///////////////////////////////////////////// 
  initFBO(); //Sets up the framebuffer object  

  //Initializes the particle positions in the domain
  //We don't need to do this anymore since we now can initialize data values
  //past 0-1 range in a RGBA32F texture without a shader.
  //pc->initParticlePositions(fbo, twidth, theight, init_shader, texid[2]); 
  

  CheckErrorsGL("END of init");
}

bool odd = true;
int p_index;
int numToEmit = 1; //Number of particles to emit in each pass
int offset = 0;

void display(void)
{
  GLint draw_buffer;
  glGetIntegerv(GL_DRAW_BUFFER, &draw_buffer);
  
  // bind the framebuffer object so we can render to the 2nd texture
  fbo->Bind();

  ////////////////////////////////////////////////////////////
  // Emit Particles
  ///////////////////////////////////////////////////////////
  if(emit){

    //Make sure there are available indices to emit particles.
    if(!indices.empty()){

      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      gluOrtho2D(0, twidth, 0, theight);
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
  
      if(odd)
	glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
      else 
	glDrawBuffer(GL_COLOR_ATTACHMENT1_EXT);

      
      //Do this for each particle that is being emitted.
      for(int i = 0; i < numToEmit; i++){
	if(!indices.empty()){
	
	  //First get available index
	  p_index = indices.back();
	  std::cout << p_index <<std::endl;
	  indices.pop_back();

	  emit_shader.activate();	

	  //Determine the coordinates into the position texture
	  int s = (p_index%twidth);
	  int t = (p_index/twidth);
	  //s = (s*(2.0/(float)twidth) - 1.0);
	  //t = (t*(2.0/(float)theight) - 1.0);
	  std::cout << s << " " << t << " " << offset <<std::endl;
	  glViewport(s,t,1,1);
       
	  glBegin(GL_POINTS);
	  {
	    glColor4f(10.0, 10.0, 10.0, 1.0);
	    glVertex2f(s, t);
	  }
	  glEnd();

	  emit_shader.deactivate();
	}
      }
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      gluOrtho2D(-1, 1, -1, 1);
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();

      //offset++;
    }
    /*if(offset == (twidth*theight)){     
      if(odd)
	glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
      else 
	glReadBuffer(GL_COLOR_ATTACHMENT1_EXT);

      pc->dumpContents(twidth, theight);
      //emit = false;
      offset++;
      }*/
  
  }
  emit = false;
  
  ////////////////////////////////////////////////////////////
  // Update Particle Positions 
  ///////////////////////////////////////////////////////////
  if (odd)
    glDrawBuffer(GL_COLOR_ATTACHMENT1_EXT);
  else 
    glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glViewport(0, 0, twidth, theight);
  
  
  glEnable(texType);
  pass1_shader.activate();
  
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(texType, texid[3]);
  glUniform1iARB(uniform_wind, 1);
  
  glActiveTexture(GL_TEXTURE0);
  glUniform1iARB(uniform_postex, 0);

  if (odd)
    glBindTexture(texType, texid[0]);  // read from texture 0
  else 
    glBindTexture(texType, texid[1]);  // read from texture 1

  glBegin(GL_QUADS);
  {
    glTexCoord2f(0, 0);            glVertex3f(-1, -1, -0.5f);
    glTexCoord2f(twidth, 0);       glVertex3f( 1, -1, -0.5f);
    glTexCoord2f(twidth, theight); glVertex3f( 1,  1, -0.5f);
    glTexCoord2f(0, theight);      glVertex3f(-1,  1, -0.5f);
  }
  glEnd();
  
  pass1_shader.deactivate();
 
  glBindTexture(texType, 0);

  CheckErrorsGL("END : after 1st pass");
  
  //Switches the frame buffer and binding texture
  odd = !odd;

  // We only need to do PASS 2 (copy to VBO) and PASS 3 (visualize) if
  // we actually want to render to the screen.  Rendering to the
  // screen will make the simulation run more slowly. This feature is
  // mainly included to allow some idea of how much faster the
  // simulation can run if left to run on the GPU.
  if (show_particle_visuals)
    {
      
      // //////////////////////////////////////////////////////////////
      // PASS 2 - copy the contents of the 2nd texture (the new positions)
      // into the vertex buffer
      // //////////////////////////////////////////////////////////////

      // In some circumstances, we may want to dump the contents of
      // the FBO to a file.
      if (dump_contents)
	{
	  pc->dumpContents(twidth, theight);
	  dump_contents = false;
	}
      
      glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB, vertex_buffer);
      glReadPixels(0, 0, twidth, theight, GL_RGBA, GL_FLOAT, 0);
      glBindBuffer(GL_PIXEL_PACK_BUFFER_ARB, 0);
      CheckErrorsGL("after glReadPixels");
      
      // Disable the framebuffer object
      FramebufferObject::Disable();
      glDrawBuffer(draw_buffer); // send it to the original buffer
      CheckErrorsGL("END : after 2nd pass");

      // //////////////////////////////////////////////////////////////
      // PASS 3 - draw the vertices; This represents the visualization
      // of the PLUME particle field.
      // //////////////////////////////////////////////////////////////

      // clear the color and depth buffer before drawing the scene, and
      // set the viewport to the window dimensions
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      glViewport(0, 0, glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
      
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      gluPerspective(60.0, 1.0, 1.0, 250.0);
      
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
      gluLookAt( eye_pos[0], eye_pos[1], eye_pos[2],
		 eye_gaze[0], eye_gaze[1], eye_gaze[2],
		 0, 1, 0 );

      if (!rotate_sphere)
	{
	  // allow rotation of this object
	  glRotatef(elevation, 1,0,0);
	  glRotatef(azimuth, 0,1,0);
	  // glTranslatef(0,0,5.0);
	}

      // render the vertices in the VBO (the particle positions) as points in the domain
      glBindBufferARB(GL_ARRAY_BUFFER, vertex_buffer);
      glEnableClientState(GL_VERTEX_ARRAY);
      glColor3f(1.0, 1.0, 1.0);
      glVertexPointer(4, GL_FLOAT, 0, 0);
      glPointSize(2.0);
      render_shader.activate();
      glDrawArrays(GL_POINTS, 0, twidth*theight); 
      render_shader.deactivate();
  
      glDisable(texType);
      CheckErrorsGL("END : visualization");
      
      // Draw axes that represent the domain dimension
      dc->drawAxes();
      
      // Draw the wind texture layer if necessary
      dc->drawLayers(visual_layer, texid[3], numInRow);
      
#ifndef WIN32
      // spit out frame rate
      if (frame_rate){
	dc->drawFrameRate(twidth, theight);
      }
#endif

      // If we've chose to display the 3D particle domain, we need to
      // set the projection and modelview matrices back to what is
      // needed for the particle advection step
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      gluOrtho2D(-1, 1, -1, 1);
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
      
      // Finally, swap the front and back buffers to display the
      // particle field to the monitor
      glutSwapBuffers();
    }
}

void idle()
{
    glutPostRedisplay();
}

void keyboard_cb(unsigned char key, int x, int y)
{
  if (key == 'z')
    eye_z -= 0.5;
  else if (key == 'Z')
    eye_z += 0.5;

  if (key == 'k') 
    {
      visual_layer--;
      if (visual_layer < -1) visual_layer = -1;
    }
  else if (key == 'K')
    {
      visual_layer++;
      if (visual_layer > ny) visual_layer = ny;
    }

  else if (key == 'd')
    {
      // toggle whether to display output
      show_particle_visuals = !show_particle_visuals;
    }

  else if (key == 27)
    {
      glutDestroyWindow(winid);
      exit(0);
    }
  else if (key == 'r')
    dump_contents = true;
  else if( key == 'e')
    emit = true;

  glutPostRedisplay();
}

static int last_x, last_y;
void mouse(int button, int state, int x, int y)
{
  last_x = x;
  last_y = y;

  if (state == GLUT_DOWN && button == GLUT_LEFT_BUTTON)
    rotate_object = true;
  else // state == GLUT_UP
    rotate_object = false;

  if (state == GLUT_DOWN && button == GLUT_RIGHT_BUTTON)
    translate_view = true;
  else // state == GLUT_UP
    translate_view = false;

  glutPostRedisplay();
}


void motion(int x, int y)
{
  if (translate_view) 
    {
      // pan view around gaze center...
      // since y is up, move eye in z only to take it into and out of the screen
      float change = y - last_y;
      eye_pos[2] = eye_pos[2] + change;
      eye_gaze[2] = eye_pos[2] - 5.0;
    }

    if (rotate_object) 
    {
	// since y is up, move eye in z only to take it into and out of the screen
	float change = x - last_x;
	float rate = 0.1;

	change = x - last_x;
	rate = 0.1;
	azimuth = azimuth + change * rate;

	change = y - last_y;
	rate = 0.1;
	elevation = elevation + change * rate;
    }

    last_x = x;
    last_y = y;

    glutPostRedisplay();
}
