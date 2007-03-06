#include <iostream>
#include <stdlib.h>

#include <GL/glew.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "plumeControl.h"
#include "framebufferObject.h"
#include "renderbuffer.h"
#include "GLSL.h"
#include "glErrorUtil.h"

Timer * display_clock;
PlumeControl* plume;

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

GLuint vertex_buffer;

int main(int argc, char** argv)
{
  int w, h, t;
  if (argc == 2)
    {
      w = atoi(argv[1]);
      h = w;
      t = 3;
    }
  else if(argc == 3)
    {
      w = atoi(argv[1]);
      h = w;
      t = atoi(argv[2]);
    }
  else{ w = 128; h = 128; t = 3;}

  plume = new PlumeControl(w,h, t);

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

FramebufferObject *fbo;
Renderbuffer *rb;

void initFBO(void){
  // Get the Framebuffer Object ready
  fbo = new FramebufferObject();
  fbo->Bind();
      
  rb = new Renderbuffer();
  rb->Set(GL_DEPTH_COMPONENT24, plume->twidth, plume->theight);
  fbo->AttachRenderBuffer(GL_DEPTH_ATTACHMENT_EXT, rb->GetId() );

  //Attach textures to framebuffer object
  fbo->AttachTexture(GL_COLOR_ATTACHMENT0_EXT, plume->texType, plume->texid[0]);
  fbo->AttachTexture(GL_COLOR_ATTACHMENT1_EXT, plume->texType, plume->texid[1]);

  fbo->IsValid();
  FramebufferObject::Disable();
}

void init(void)
{
  glEnable(GL_DEPTH_TEST);
  
  plume->init();

  display_clock = new Timer(true);

  //
  // set up vertex buffer
  // 
  glGenBuffersARB(1, &vertex_buffer);
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, vertex_buffer);
  glBufferDataARB(GL_ARRAY_BUFFER_ARB, plume->twidth*plume->theight*4*sizeof(GLfloat), 0, GL_STREAM_COPY);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  // ///////////////////////////////////////////// 
  initFBO(); //Sets up the framebuffer object  

  //Initializes the particle positions in the domain
  //We don't need to do this anymore since we now can initialize data values
  //past 0-1 range in a RGBA32F texture without a shader.
  //plume->pc->initParticlePositions(fbo, plume->texid[2]); 
  

  CheckErrorsGL("END of init");
}

bool odd = true;
Timer_t display_time[2];

void display(void)
{
  GLint draw_buffer;
  glGetIntegerv(GL_DRAW_BUFFER, &draw_buffer);
  
  // bind the framebuffer object so we can render to the 2nd texture
  fbo->Bind();

  ////////////////////////////////////////////////////////////
  // Emit Particles
  ///////////////////////////////////////////////////////////
  //record end time
  display_time[1] = display_clock->tic();
  
  if(emit){
    float runtime = display_clock->deltas(display_time[0],display_time[1]);
    plume->time_step = runtime;
    //std::cout << runtime << std::endl;
    plume->injectParticles(fbo, odd);
  
  }
  //record start time
  display_time[0] = display_clock->tic();
  //emit = false;
  ////////////////////////////////////////////////////////////
  // Update Particle Positions 
  ////////////////////////////////////////////////////////////
  plume->advectParticles(fbo,odd);

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
	  plume->pc->dumpContents();
	  dump_contents = false;
	}
      
      glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB, vertex_buffer);
      glReadPixels(0, 0, plume->twidth, plume->theight, GL_RGBA, GL_FLOAT, 0);
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

      plume->displayVisual(vertex_buffer);
     
      glDisable(plume->texType);
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
      plume->dc->decreaseVisualLayer();
    }
  else if (key == 'K')
    {
      plume->dc->increaseVisualLayer();
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
    plume->dc->rotate_object = true;
  else // state == GLUT_UP
    plume->dc->rotate_object = false;

  if (state == GLUT_DOWN && button == GLUT_RIGHT_BUTTON)
    plume->dc->translate_view = true;
  else // state == GLUT_UP
    plume->dc->translate_view = false;

  glutPostRedisplay();
}


void motion(int x, int y)
{
  if (plume->dc->translate_view) 
    {
      // pan view around gaze center...
      // since y is up, move eye in z only to take it into and out of the screen
      float change = y - last_y;
      plume->dc->setEyeValues(change);
      
    }

    if (plume->dc->rotate_object) 
    {
	// since y is up, move eye in z only to take it into and out of the screen
	float change = x - last_x;
	float rate = 0.1;

	change = x - last_x;
	rate = 0.1;
	plume->dc->setAzimuth(change,rate);

	change = y - last_y;
	rate = 0.1;
	plume->dc->setElevation(change,rate);
    }

    last_x = x;
    last_y = y;

    glutPostRedisplay();
}
