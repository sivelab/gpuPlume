#include <iostream>
#include <stdlib.h>

#include <GL/glew.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#ifdef WIN32
#include <windows.h>
#include <stdio.h>
#include <conio.h>
#include <tchar.h>
#endif

#include "plumeControl.h"

PlumeControl* plume;

void init(void);
void idle();
void reshape(int w, int h);
void display(void);
void keyboard_cb(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);

int winid;

int winwidth = 512, winheight = 512;

int main(int argc, char** argv)
{
  int w, h, t;
  if (argc == 2)
    {
      w = atoi(argv[1]);
      h = w;
      t = 4;
    }
  else if(argc == 3)
    {
      w = atoi(argv[1]);
      h = w;
      t = atoi(argv[2]);
    }
  else{ w = 128; h = 128; t = 4;}

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

void init(void)
{
  glEnable(GL_DEPTH_TEST);
  
  plume->init(); 
 
}

void display(void)
{
  plume->display();
 
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
      plume->show_particle_visuals = !plume->show_particle_visuals;
    }

  else if (key == 27)
    {
      glutDestroyWindow(winid);
      exit(0);
    }
  else if (key == 'r')
    {
      plume->dump_contents = true;
    }
  else if( key == 'e')
    {
      plume->emit = !plume->emit;
    }
  /*else if(key == 'c')
    {
      plume->collectionBoxes = !plume->collectionBoxes;
    }
  else if(key == 'b')
    {
      plume->output_CollectionBox = true;
      }*/

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
