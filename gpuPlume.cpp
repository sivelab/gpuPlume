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
#include "particleEmitter.h"
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
ParticleEmitter* pe;

int numInRow;   // represents the number of layers that are stored in one row in a 2D texture

bool dump_contents = false;
bool emit = false;
bool show_particle_visuals = true;

void init(void);
void initFBO(void);
void idle();
void reshape(int w, int h);
void display(void);
void keyboard_cb(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);

int winid;

int winwidth = 512, winheight = 512;
//
// The values here determine the number of particles
//
int twidth = 128, theight = 128;

GLuint texid[8];
GLSLObject emit_shader;
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

  pc = new ParticleControl(texType,twidth,theight);
  pc->getDomain(&nx,&ny,&nz);

  dc = new DisplayControl(nx, ny, nz, texType);

  //Create a particleEmitter with position 10,10,10
  pe = new ParticleEmitter(10.0, 10.0, 10.0, &twidth, &theight, &indices);

  for(int i = twidth*theight-1; i >= 0; i--)
    indices.push_back(i);

  // load up textures to hold position
  CheckErrorsGL("BEGIN : Creating textures");

  glEnable(texType);
  glGenTextures(8, texid); // create (reference to) a new texture
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

  //
  // create random texture for use with particle simulation and turbulence
  //
  for (int j=0; j<theight; j++)
    for (int i=0; i<twidth; i++)
      {
	int idx = j*twidth*sz + i*sz;
	
	//
	// Generate random values should of normal distribution with zero mean and standard deviation of one.
	// Need to pull classes from sim_fast that handle this... 
	// For now, generate random values between -1 and 1.... shader subtracts 1.0
	//
	data[idx] = randVal() * 2.0 - 1.0;
	data[idx+1] = randVal() * 2.0 - 1.0;
	data[idx+2] = randVal() * 2.0 - 1.0;
	data[idx+3] = 0.0;
      }
  pc->createWrappedTexture(texid[4], int_format, twidth, theight, data);

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
  //This shader is used to advect the particles using the windfield
  pc->setupAdvectShader(&time_step, &numInRow);

  //This shader is used to emmit particles
  emit_shader.addShader("Shaders/emitParticle_vp.glsl", GLSLObject::VERTEX_SHADER);
  emit_shader.addShader("Shaders/emitParticle_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  emit_shader.createProgram();

  // ///////////////////////////////////////////// 
  initFBO(); //Sets up the framebuffer object  

  //Initializes the particle positions in the domain
  //We don't need to do this anymore since we now can initialize data values
  //past 0-1 range in a RGBA32F texture without a shader.
  //pc->initParticlePositions(fbo, texid[2]); 
  

  CheckErrorsGL("END of init");
}

bool odd = true;
//int numToEmit = 1; //Number of particles to emit in each pass

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

      pe->EmitParticle(fbo, emit_shader, odd);
  
  }
  emit = false;
  ////////////////////////////////////////////////////////////
  // Update Particle Positions 
  ////////////////////////////////////////////////////////////

  pc->advect(fbo, odd, texid[4], texid[3], texid[0], texid[1]);

  ////////////////////////////////////////////////////////////

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
	  pc->dumpContents();
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

      dc->drawVisuals(vertex_buffer, texid[3], numInRow, twidth, theight);
        
      glDisable(texType);
      CheckErrorsGL("END : visualization");
      
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
  if (key == 'k') 
    {
      dc->visual_layer--;
      if (dc->visual_layer < -1) dc->visual_layer = -1;
    }
  else if (key == 'K')
    {
      dc->visual_layer++;
      if (dc->visual_layer > ny) dc->visual_layer = ny;
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
    dc->rotate_object = true;
  else // state == GLUT_UP
    dc->rotate_object = false;

  if (state == GLUT_DOWN && button == GLUT_RIGHT_BUTTON)
    dc->translate_view = true;
  else // state == GLUT_UP
    dc->translate_view = false;

  glutPostRedisplay();
}


void motion(int x, int y)
{
  if (dc->translate_view) 
    {
      // pan view around gaze center...
      // since y is up, move eye in z only to take it into and out of the screen
      float change = y - last_y;
      dc->setEyeValues(change);
      
    }

    if (dc->rotate_object) 
    {
	// since y is up, move eye in z only to take it into and out of the screen
	float change = x - last_x;
	float rate = 0.1;

	change = x - last_x;
	rate = 0.1;
	dc->setAzimuth(change,rate);

	change = y - last_y;
	rate = 0.1;
	dc->setElevation(change,rate);
    }

    last_x = x;
    last_y = y;

    glutPostRedisplay();
}
