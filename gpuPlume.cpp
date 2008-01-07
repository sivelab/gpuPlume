#include <iostream>

#include <GL/glew.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <stdlib.h>
#include <GL/glut.h>
#endif

#ifdef WIN32
#include <windows.h>
#include <stdio.h>
#include <conio.h>
#include <tchar.h>
#endif

#include "plumeControl.h"
#include "nonGaussianModel.h"
#include "GaussianModel.h"
#include "Gaussian_2shaders_Model.h"
#include "ReflectionModel.h"
#include "MultipleBuildingsModel.h"
#include "GeomTest.h"

#include "Timer.h"

// Variables to hold timing array and record timings
long timing_count = 0;
bool compute_timings = false;
long timing_N;
double *timing_array;
Timer_t total_timer[2];
Timer_t plume_timer[2];
Timer *plume_clock;

Util* util;

PlumeControl* plume;
int curr;

void init(void);
void idle();
void reshape(int w, int h);
void display(void);
void keyboard_cb(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);

int winwidth = 512, winheight = 512;

//int last_x, last_y;

int main(int argc, char** argv)
{
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
  
  // allocate memory for the timing values
  // keep 1000 values
  timing_N = 1000;
  timing_array = new double[timing_N];

  // Create a timer to use with timing calculations
  // & query the start time so we have a value in it
  // 
  // Attempt to create a high resolution timer.
  plume_clock = new Timer(true);
  plume_timer[0] = plume_clock->tic();

  util = new Util();
  util->readInput("Settings/input.txt");

  //plume = new PlumeControl(util);
  //Options are:
  //  nonGaussianModel();
  //  GaussianModel();
  switch(util->advectChoice){
  case 0:
    //util->windFieldData = 4;
    plume = new Gaussian_2shaders_Model(util);
    break;
  case 1:
    //util->windFieldData = 4;
    plume = new GaussianModel(util);
    break;
  case 2:
    //util->windFieldData = 5;
    plume = new NonGaussianModel(util);
    break;
  case 3:
    plume = new ReflectionModel(util);
    break;
  case 4:
    plume = new MultipleBuildingsModel(util);
    break;
  case 5:
    plume = new GeomTest(util);
    break;
  default:
    std::cout << "Error in advection Choice in Settings file!" << std::endl;
  }

  Random random_gen(2);

  glutInitDisplayMode( GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE );
  glutInit(&argc, argv);

  // A variable to allow us (eventually) to control full screen, hopefully, SLI rendering.
  bool use_game_mode = false;
  if (use_game_mode) 
    {
      glutGameModeString("1280x1024");
      std::cout << "Game Mode Width: " << glutGameModeGet(GLUT_GAME_MODE_WIDTH) << std::endl;
      std::cout << "Game Mode Height: " << glutGameModeGet(GLUT_GAME_MODE_HEIGHT) << std::endl;
      glutEnterGameMode();
    }
  else 
    {
      glutInitWindowSize(winwidth, winheight);
      plume->winid = glutCreateWindow("gpuPLUME");
    }

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

  // We really should place a function call here to do a scan and
  // check over the various OpenGL extensions and states that we
  // require for this code to function.
  if (!GL_ARB_vertex_buffer_object) 
    {
      std::cout << "GL_ARB_vertex_buffer_object is NOT available!  Exiting!" << std::endl;
      exit(-1);
    }
/*
  if(GL_EXT_geometry_shader4){
    std::cout << "Ready for geom shader!" << std::endl;
    int numV = 0;
    glGetIntegerv(GL_MAX_GEOMETRY_OUTPUT_VERTICES_EXT, &numV);
    std::cout << numV << " number of output vertices." << std::endl;
  }
  else
    std::cout << "NOT ready for geom shader." << std::endl;
  */

  init();

  // record the start time
  total_timer[0] = plume_clock->tic();  

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
  curr = 0;
  glEnable(GL_DEPTH_TEST);
  
  plume->init(false); 
  plume->paused = false;
  plume->inPauseMode = util->pauseMode;

}

void display(void)
{
  if (compute_timings)
    {
      plume_timer[0] = plume_clock->tic();      
    }

  // if quitting the simulation, 0 is returned
  int quitSimulation = 1;
  quitSimulation = plume->display();
  
  if (quitSimulation == 0)
    {
      total_timer[1] = plume_clock->tic();  
      std::cout << "Total Simulation Time: " << plume_clock->deltas(total_timer[0], total_timer[1]) << " seconds." << std::endl;

#ifdef WIN32
	  system("pause");
#endif

      // system("pause");
      glutDestroyWindow(plume->winid);
      exit(0);      
    }
 
  // stats computation
  if (compute_timings)
    {
      plume_timer[1] = plume_clock->tic();      
      if (timing_count >= 0 && timing_count < timing_N)
	{
	  // calculate the difference between start and end, returns end - start
	  timing_array[timing_count] = plume_clock->deltas(plume_timer[0], plume_timer[1]);
	}
      
      // increment timing value
      timing_count++;
	  
      if (timing_count == timing_N) 
	{
	  std::cout << "Timing complete, restoring visuals..." << std::endl;
	  compute_timings = false;
	  util->show_particle_visuals = true;
	}
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
      plume->contours->decreaseContourLayer();
    }
  else if (key == 'K')
    {
      plume->dc->increaseVisualLayer();
      plume->contours->increaseContourLayer();
    }
  else if (key == 'l')
    {
      plume->dc->visual_field++;
      if(plume->dc->visual_field > 4)
	plume->dc->visual_field = 0;

      //value of 1 displays Tau11
      //value of 2 displays Tau22
      //value of 3 displays Tau33
      //value of 4 displays Tau13
      plume->contours->tauValue = plume->dc->visual_field -1;
      

    }
  else if (key == '<')
    {
      plume->dc->moveSliderDown();
    }
  else if (key == '>')
    {
      plume->dc->moveSliderUp();
    }
  else if (key == 't')
    {
      // toggle whether to display output
      util->show_particle_visuals = !util->show_particle_visuals;
    }

  else if (key == 'b')
    {
      // Turn on timing - shut off visuals and allow a ten cycle count before timings start
      util->show_particle_visuals = false;
      timing_count = -10;
      compute_timings = true;
    }
  else if (key == 27)
    {
      total_timer[1] = plume_clock->tic();  
      std::cout << "Total Simulation Time: " << plume_clock->deltas(total_timer[0], total_timer[1]) << " seconds." << std::endl;

      // Before we exit, write out timing values (if collected) to a file.
      std::cout << "Exiting... writing out timing values..." << std::endl;
      std::ofstream outfile("gpuplumetimes.m");
      outfile << "plumetimes = [" << std::endl;
      for (long i=0; i<timing_N; i++)
	{
	  outfile << timing_array[i] << ";" << std::endl;
	}
      outfile << "];" << std::endl;

      // 
      // Before deleting the timing array, perform the mean, std, and
      // var calculations.
      //
      double mean = 0.0;
      for (long i=0; i<timing_N; i++)
	mean += timing_array[i];
      mean = mean/(double)timing_N;
      std::cout << "Mean Advection Step Time: " << mean << " seconds (" 
		<< mean * 1000.0 << " ms, " 
		<< mean * 1000000.0 << " us)" << std::endl;

      double isum = 0.0;
      for (long i=0; i<timing_N; i++)
	isum += ((timing_array[i] - mean) * (timing_array[i] - mean));
      double sigma = sqrt(1.0/(double)timing_N * isum);
      std::cout << "\tStandard Deviation: " << sigma << ", Variance = " << sigma*sigma << std::endl;

      // Delete the timing array
      delete [] timing_array;

#ifdef WIN32
	  system("pause");
#endif

      // Next, destroy the glut window and exit
      glutDestroyWindow(plume->winid);
      exit(0);
    }
  else if (key == 'r')
    {
      plume->dump_contents = true;
    }
  else if( key == 'e')
    {
      plume->pe[curr]->emit = !plume->pe[curr]->emit;
    }
  else if( key == 'o')
    {
      plume->pe[curr]->releaseType = onePerKeyPress;
    }
  else if(key == 'p')
    {
      plume->pe[curr]->releaseType = perSecond;
    }
  else if (key == 'W')
    {
      plume->pe[curr]->setPosition(plume->pe[curr]->xpos, plume->pe[curr]->ypos-1.0, plume->pe[curr]->zpos);
    }
  else if (key == 'S')
    {
      plume->pe[curr]->setPosition(plume->pe[curr]->xpos, plume->pe[curr]->ypos+1.0, plume->pe[curr]->zpos);
    }
   else if (key == 'A')
    {
      plume->pe[curr]->setPosition(plume->pe[curr]->xpos-1.0, plume->pe[curr]->ypos, plume->pe[curr]->zpos);
    }
  else if (key == 'D')
    {
      plume->pe[curr]->setPosition(plume->pe[curr]->xpos+1.0, plume->pe[curr]->ypos, plume->pe[curr]->zpos);
    }
  else if (key == '>')
    {
      curr++;
      if(curr == util->numOfPE)
	curr = 0;

    }
  else if(key == 'c')
    {
      //plume->stream->clear();   
      plume->pathLines->clear();
    }
  else if(key == 'g')
    {
      plume->pc->outputPrime = true;
    }
  else if(key == 'h')
    {
      plume->createImages = true;
    }
  else if(key == 'm')
    {
      plume->print_MeanVel = true;
    }
  else if(key == 'd')
    {
      plume->dc->slideLeftorRight(1.0);
    }
  else if(key == 'a')
    {
      plume->dc->slideLeftorRight(-1.0);
    }
  else if(key == ' ')
    {
      plume->paused = false;
    }
  else if(key == 'z')
    {
      plume->inPauseMode = !plume->inPauseMode;
    }
  else if(key == 'w')
    {
      plume->dc->moveForwardorBack(1.0);
    }
  else if(key == 's')
    {
      plume->dc->moveForwardorBack(-1.0);
    }
  else if(key == '1')
    {
      plume->draw_contours = true;
      plume->draw_layers = false;
    }
  else if(key == '2')
    {
      plume->draw_contours = false;
      plume->draw_layers = true;
    }

  glutPostRedisplay();
}

static int last_x, last_y;
void mouse(int button, int state,int x, int y)
{
  last_x = x;
  last_y = y;

  if (state == GLUT_DOWN && button == GLUT_MIDDLE_BUTTON)
    plume->dc->change_height = true;
  else // state == GLUT_UP
    plume->dc->change_height = false;

  if (state == GLUT_DOWN && button == GLUT_RIGHT_BUTTON)
    plume->dc->change_look = true;
  else // state == GLUT_UP
    plume->dc->change_look = false;

  if (state == GLUT_DOWN && button == GLUT_LEFT_BUTTON) {
    plume->dc->clickedRangeButton(x,y);

    if(plume->dc->clickedSlider(x,y)){
      plume->dc->move_slider = true;
      plume->dc->rotate_around = false;
    }
    else{
      plume->dc->rotate_around = true;
      plume->dc->move_slider = false;
    }
   
  }
  else{
    plume->dc->rotate_around = false;
    plume->dc->move_slider = false;
  }

  glutPostRedisplay();
}


void motion(int x, int y)
{
  
  float change_y = y - last_y;
  float change_x = x - last_x;

  //std::cout << "last: " << last_y << "   current: " << y << std::endl;

  //plume->dc->lookUporDown(change_y);
  //plume->dc->setRotateAround(change_x);

  if(plume->dc->move_slider){
    plume->dc->moveSlider(x);
  }

  if (plume->dc->change_height) 
    {
      // pan view around gaze center...
      // since y is up, move eye in z only to take it into and out of the screen
      //float change = y - last_y;
      //plume->dc->moveForwardorBack(change_x);
      plume->dc->setElevation(change_y,0.1);
    }

    if (plume->dc->change_look) 
    {
	// since y is up, move eye in z only to take it into and out of the screen
	//float change = x - last_x;
	//float rate = 0.1;

	//change = x - last_x;
	//rate = 0.1;
	//plume->dc->setAzimuth(change,rate);

	//change = y - last_y;
      //	rate = 0.1;
      //plume->dc->setElevation(change,rate);
       plume->dc->lookUporDown(change_y);

    }

    if (plume->dc->rotate_around)
    {
      //float change = x - last_x;
      //float rate = 0.1;
     
      plume->dc->setRotateAround(change_x);
      

    }
    last_x = x;
    last_y = y;
    

    glutPostRedisplay();
}
