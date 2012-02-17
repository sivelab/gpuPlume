#include <iostream>
#include <signal.h>

#include <GL/glew.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <cstdlib>
#include "GL/gl.h"      // System and OpenGL Stuff
#include "GL/glut.h"
#include "GL/glx.h"
#endif

#ifdef WIN32
#include <windows.h>
#include <stdio.h>
#include <conio.h>
#include <tchar.h>
#endif

#include "util/handlePlumeArgs.h"

#include "plumeControl.h"
// #include "nonGaussianModel.h"
// #include "GaussianModel.h"
// #include "Gaussian_2shaders_Model.h"
// #include "ReflectionModel.h"
#include "MultipleBuildingsModel.h"
// #include "GeomTest.h"

#include "Timer.h"
#include "glErrorUtil.h"

#include "CmdOptionInterpreter.h"

#define PBUFFER_WIDTH  1024
#define PBUFFER_HEIGHT 1024

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

void signalHandler(int sig);

void idle();
void reshape(int w, int h);
void display(void);
void keyboard_cb(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);

int winwidth = 1000, winheight = 1000;
// int winwidth = 1000, winheight = (int)(1000 * 9/16.0);

//int last_x, last_y;

#ifdef __linux__
// Store all system info in one place
typedef struct RenderContextRec
{
  GLXContext ctx;
  Display *dpy;
  GLXContext g_pbufferContext;
  GLXPbuffer g_pbuffer;
    Window win;
    int nWinWidth;
    int nWinHeight;
    int nMousePosX;
    int nMousePosY;

} RenderContext;

void init(RenderContext *rcx);
void Cleanup(RenderContext *rcx);

void CreateWindow(RenderContext *rcx)
{
    XSetWindowAttributes winAttribs;
    GLint winmask;
    GLint nMajorVer = 0;
    GLint nMinorVer = 0;
    XVisualInfo *visualInfo;
    static int attributes[] = { GLX_RGBA,
			        GLX_DOUBLEBUFFER,
			        GLX_RED_SIZE, 8,
			        GLX_BLUE_SIZE, 8,
			        GLX_GREEN_SIZE, 8,
			        0 };

    // Tell X we are going to use the display
    rcx->dpy = XOpenDisplay(NULL);

    // Get Version info
    glXQueryVersion(rcx->dpy, &nMajorVer, &nMinorVer);
    std::cout << "Supported GLX version - " << nMajorVer << "." << nMinorVer << std::endl;

    if(nMajorVer == 1 && nMinorVer < 2)
    {
      std::cerr << "ERROR: GLX 1.2 or greater is necessary\n";
        XCloseDisplay(rcx->dpy);
        exit(0);
    }

    // Get a new visual that meets our attrib requirements
    visualInfo = glXChooseVisual(rcx->dpy, DefaultScreen(rcx->dpy), attributes);

    // Now create an X window
    winAttribs.event_mask = ExposureMask | VisibilityChangeMask | 
                            KeyPressMask | PointerMotionMask    |
                            StructureNotifyMask ;

    winAttribs.border_pixel = 0;
    winAttribs.bit_gravity = StaticGravity;
    winAttribs.colormap = XCreateColormap(rcx->dpy, 
                                          RootWindow(rcx->dpy, visualInfo->screen), 
                                          visualInfo->visual, AllocNone);
    winmask = CWBorderPixel | CWBitGravity | CWEventMask| CWColormap;

    rcx->win = XCreateWindow(rcx->dpy, DefaultRootWindow(rcx->dpy), 20, 20,
			     rcx->nWinWidth, rcx->nWinHeight, 0, 
                             visualInfo->depth, InputOutput,
			     visualInfo->visual, winmask, &winAttribs);

    XMapWindow(rcx->dpy, rcx->win);

    // hopefully set the window behind all others
    XLowerWindow(rcx->dpy, rcx->win);

    // Also create a new GL context for rendering
    rcx->ctx = glXCreateContext(rcx->dpy, visualInfo, 0, True);
    glXMakeCurrent(rcx->dpy, rcx->win, rcx->ctx);
}
#endif // __linux__

int main(int argc, char** argv)
{
  // Setup a signal handler to catch the ^C when program exits
  signal(SIGINT, signalHandler);

#ifdef WIN32
  if (argc == 1)
  {	
    // Fix! -Pete
    // only do this if we have not been supplied a command line argument ... in other words, if we are likely running from Visual Studio!
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
  }
#endif
 
  // 
  // Setup command line argument parsing
  // 
  sivelab::PlumeArgs argParser;

  // args to set the problem id so we can tailor files names appropriately
  argParser.reg("problemID", "specifies the problem ID to use for the optimization run", sivelab::ArgumentParsing::INT, 'b');
  argParser.reg("probInstID", "specifies the problem instance ID to use for optimization", sivelab::ArgumentParsing::INT, 'z');

  // allocate memory for the timing values
  // keep 1000 values
  timing_N = 1000;
  timing_array = new double[timing_N];

  // Create a timer to use with timing calculations
  // & query the start time so we have a value in it
  // 
  // Attempt to create a high resolution timer.
  // plume_clock = new Timer(false);
  plume_clock = new Timer(true);
  plume_timer[0] = plume_clock->tic();

  util = new Util();

  // Parse any command line options specifed.
  argParser.process(argc, argv);  

  // Must supply an argument containing the .prof file to be read.
  std::string quicInputFile;
  if (argParser.isSet("quicproj", quicInputFile))
    {
      std::cout << "Reading input from file: \"" << quicInputFile << "\"" << std::endl;
      
      if (util->readInput( quicInputFile ) == false)
	{
	  std::cerr << "Could not open or parse input file: \"" << argv[1] << "\"\nExiting." << std::endl;
	  argParser.printUsage();
	  exit(EXIT_FAILURE);
	}
    }
  else 
    {
      std::cerr << "Need to provide QUIC .proj file for opening. Exiting." << std::endl;
      argParser.printUsage();
      exit(EXIT_FAILURE);
    }
  
  CmdOptionInterpreter cmdOI(&argParser, util);
  cmdOI.parse();

  // no Advect choice
#if 0
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
#endif

  long int seedVal = (long)time(0) % (long)getpid();
  Random random_gen( seedVal );

  plume = new MultipleBuildingsModel(util);

#ifdef __linux__
  RenderContext rcx;
#endif // __linux__

  if (util->offscreenRender == false)
    {
      glutInitDisplayMode( GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE );
      glutInit(&argc, argv);

      // A variable to allow us (eventually) to control full screen, hopefully, SLI rendering.
      bool use_game_mode = util->fullscreen;
      if (use_game_mode) 
	{
	  //glutGameModeString("1280x1024");
	  glutGameModeString("1024x768");
	  std::cout << "Game Mode Width: " << glutGameModeGet(GLUT_GAME_MODE_WIDTH) << std::endl;
	  std::cout << "Game Mode Height: " << glutGameModeGet(GLUT_GAME_MODE_HEIGHT) << std::endl;
	  util->winWidth = glutGameModeGet(GLUT_GAME_MODE_WIDTH);
	  util->winHeight = glutGameModeGet(GLUT_GAME_MODE_HEIGHT);
	  glutEnterGameMode();
	}
      else 
	{
	  float domainAspectRatio = (float)util->ny / (float)util->nx;
	  std::cout << "setting winsize = " << winwidth << " X " << winwidth/domainAspectRatio << std::endl;

	  util->winWidth = winwidth;
	  util->winHeight = winwidth / domainAspectRatio;

	  glutInitWindowSize(winwidth, winwidth/domainAspectRatio);
	  plume->winid = glutCreateWindow("gpuplume");
	}

      glutDisplayFunc(display);
      glutReshapeFunc(reshape);
      glutIdleFunc(idle);
      glutKeyboardFunc(keyboard_cb);
      glutMotionFunc(motion);
      glutMouseFunc(mouse);
    }
  else
    {
#ifdef __linux__
      // ////////////////////////////////////////////////////////////
      // IMPORTANT!!!
      // ////////////////////////////////////////////////////////////

      // This is X-based Windowing!!! it is only used for the
      // off-screen rendering (aka headless rendering) of the gpuplume
      // code.  This most certainly will not work on Windows and may
      // be problematic on OS X.  Therefore, by default, glut-based
      // windowing will be used.

      // Bool bWinMapped = False;

      float domainAspectRatio = (float)util->ny / (float)util->nx;

      // Set initial window size... don't need big window as we're
      // rendering to textures.
      rcx.nWinWidth = winwidth;
      rcx.nWinHeight = rcx.nWinWidth / domainAspectRatio;

      util->winWidth = rcx.nWinWidth;
      util->winHeight = rcx.nWinHeight;

      // Setup X window and GLX context
      CreateWindow(&rcx);
#else
      std::cerr << "Offscreen Rendering not supported on this architecture/OS!" << std::endl;
      exit(EXIT_FAILURE);
#endif // __linux__
    }

  GLenum err = glewInit();
  if (GLEW_OK != err) 
    {
      std::cout << "Error: " << glewGetErrorString(err) << std::endl;
    }
  
  // We really should place a function call here to do a scan and
  // check over the various OpenGL extensions and states that we
  // require for this code to function.
  //  if (!GL_ARB_vertex_buffer_object) 
  //    {
  //      std::cout << "GL_ARB_vertex_buffer_object is NOT available!  Exiting!" << std::endl;
  //      exit(-1);
  //    }

  curr = 0;
  glEnable(GL_DEPTH_TEST);
  
  plume->init(false); 
  plume->paused = false;
  plume->inPauseMode = util->pauseMode;

  plume->dc->increaseVisualLayer();

  // record the start time
  total_timer[0] = plume_clock->tic();  

  if (util->offscreenRender == false)
    {
      glutMainLoop();
    }
  else 
    {
#ifdef __linux__
      init(&rcx);

      // Execute loop the whole time the app runs
      for(;;)
	{
	  // XEvent newEvent;
	  // XWindowAttributes winData;

#if 0
	  // Watch for new X events
	  XNextEvent(rcx.dpy, &newEvent);

	  switch(newEvent.type)
	    {
	      case UnmapNotify:
		bWinMapped = False;
		break;
	      case MapNotify :
		bWinMapped = True;
	      case ConfigureNotify:
		XGetWindowAttributes(rcx.dpy, rcx.win, &winData);
		rcx.nWinHeight = winData.height;
		rcx.nWinWidth = winData.width;
		// SetupGLState(&rcx);
		break;
	      case MotionNotify:
		// display();
		break;
	      case KeyPress:
	      case DestroyNotify:
		Cleanup(&rcx);
		exit(0);
		break;
	    }
#endif
	  // if(bWinMapped)
	  // {
	  // Draw(&rcx);
	  display();
	  // }

	  // Display rendering
	  glXSwapBuffers(rcx.dpy, rcx.win);    
	}

      Cleanup(&rcx);
#endif // __linux__
    }

  return 0;
}

#ifdef __linux__
void Cleanup(RenderContext *rcx)
{
    // Unbind the context before deleting
    glXMakeCurrent(rcx->dpy, None, NULL);

    glXDestroyContext(rcx->dpy, rcx->ctx);
    rcx->ctx = NULL;

    XDestroyWindow(rcx->dpy, rcx->win);
    rcx->win = (Window)NULL;

    XCloseDisplay(rcx->dpy);
    rcx->dpy = 0;
}

void init(RenderContext *rcx)
{
	//
	// Create the p-buffer...
	//
    
    int scrnum;
    GLXFBConfig *fbconfig;
    XVisualInfo *visinfo;
    int nitems;

    int attrib[] = 
    {
        GLX_DOUBLEBUFFER,  False,
        GLX_RED_SIZE,      1,
        GLX_GREEN_SIZE,    1,
        GLX_BLUE_SIZE,     1,
        GLX_DEPTH_SIZE,    1,
        GLX_RENDER_TYPE,   GLX_RGBA_BIT,
        GLX_DRAWABLE_TYPE, GLX_PBUFFER_BIT | GLX_WINDOW_BIT,
        None
    };

    int pbufAttrib[] = 
    {
        GLX_PBUFFER_WIDTH,   PBUFFER_WIDTH,
        GLX_PBUFFER_HEIGHT,  PBUFFER_HEIGHT,
        GLX_LARGEST_PBUFFER, False,
        None
    };

    scrnum = DefaultScreen( rcx->dpy );

    fbconfig = glXChooseFBConfig( rcx->dpy,
                                  scrnum,
                                  attrib,
                                  &nitems );
    if( NULL == fbconfig )
    {
        cout << "Error: init - Couldn't get fbconfig" << endl;
        exit( 1 );
    }

      rcx->g_pbuffer = glXCreatePbuffer( rcx->dpy, fbconfig[0], pbufAttrib );

    visinfo = glXGetVisualFromFBConfig( rcx->dpy, fbconfig[0] );

    if( !visinfo )
    {
        cout << "Error: init - Couldn't get an RGBA, double-buffered visual" << endl;
        exit( 1 );
    }

	//
	// In this sample, we're going to speed things up by sharing display lists 
	// and textures between rendering contexts.
	//
	// With sharing turned on, we can simple create our dynamic texture by 
	// binding it to the p-buffer and calling glCopyTexSubImage2D, but this
	// will only work if the window's rendering context is sharing its display 
	// lists and texture space with the p-buffer's context.
	//
	
      rcx->g_pbufferContext = glXCreateContext( rcx->dpy,
                                         visinfo,
                                         rcx->ctx, // Share display lists and textures with the regular window
                                         GL_TRUE );

    if( !rcx->g_pbufferContext )
    {
        cout << "ERROR: init - Call to glXCreateContext failed!" << endl;
        exit( 1 );
     }

    XFree( fbconfig );
    XFree( visinfo );
    
    //
	// We were successful in creating a p-buffer. We can now make its context 
	// current and set it up just like we would a regular context 
	// attached to a window.
	//
	
    glXMakeCurrent( rcx->dpy, rcx->g_pbuffer, rcx->g_pbufferContext );



}
#endif // __linux__


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

void display(void)
{
  // For now, turn on visuals for displaying information
  plume->dc->tau_visual = draw_layers;

  if (compute_timings)
    {
      plume_timer[0] = plume_clock->tic();    
      //std::cout << timing_count << std::endl;
    }

  // if quitting the simulation, 0 is returned
  int quitSimulation = 1;

  Timer_t displayStart = plume_clock->tic();    

  plume->pe[curr]->emit = true;
  plume->pe[curr]->releaseType = onePerKeyPress;

  quitSimulation = plume->display();

  Timer_t displayEnd = plume_clock->tic();    
  // std::cout << "Display Time: " << plume_clock->deltam(displayStart, displayEnd) << " ms." << std::endl;  
  
  if (quitSimulation == 0)
    {
      total_timer[1] = plume_clock->tic();  
      std::cout << "Total Simulation Time: " << plume_clock->deltas(total_timer[0], total_timer[1]) << " seconds." << std::endl;

#ifdef WIN32
	// We do not want to do this when we run this through optimization experiments...  
	//  system("pause");
#endif

      // glutDestroyWindow(plume->winid);
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
	  //std::cout << timing_count << std::endl;
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
      plume->planeVisual->decreasePlaneLayer();
    }
  else if (key == 'K')
    {
      plume->dc->increaseVisualLayer();
      plume->contours->increaseContourLayer();
      plume->planeVisual->increasePlaneLayer();

    }
  else if (key == 'X')
    {
      plume->planeVisual->increasePitch();
    }
  else if (key == 'x')
    {
      plume->planeVisual->decreasePitch(); 
    }
  else if (key == 'Y')
    {
      plume->planeVisual->increaseYaw();
    }
  else if (key == 'y')
    {
      plume->planeVisual->decreaseYaw();
    }
  else if (key == 'R')
    {
      plume->planeVisual->increaseRoll();
    }
  else if (key == 'r')
    {
      plume->planeVisual->decreaseRoll();
    }

  else if (key == 'l')
    {
      plume->planeVisual->visual_field++;
      if(plume->planeVisual->visual_field > 4)
	plume->planeVisual->visual_field = 0;

      //value of 1 displays Tau11
      //value of 2 displays Tau22
      //value of 3 displays Tau33
      //value of 4 displays Tau13
      plume->contours->tauValue = plume->planeVisual->visual_field -1;
      

    }
  else if (key == '<')
    {
      plume->planeVisual->moveSliderDown();
    }
  else if (key == '>')
    {
      plume->planeVisual->moveSliderUp();
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
	// Note: We do not want to do this with optimization experiments as it will cause the system to 
	  // wait for the user to press a key... bad for thousands of simulation runs!
	 // system("pause");
#endif

      // Next, destroy the glut window and exit
      glutDestroyWindow(plume->winid);
      exit(0);
    }
  else if (key == 'f')
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
  else if (key == '+')
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
      plume->swapPauseMode();
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
      plume->dc->tau_visual = draw_contours;
    }
  else if(key == '2')
    {
      plume->dc->tau_visual = draw_layers;
    }
  else if(key == '3')
    {
      plume->planeVisual->rotationPlane = !plume->planeVisual->rotationPlane;
    }
  else if(key == '-')
    {
      plume->planeVisual->switchPlane();
      plume->contours->switchPlane();
    }
  else if (key == '.')
    {
      plume->dc->perform_cpu_sort = !plume->dc->perform_cpu_sort;
      std::cout << "Sorting: " << plume->dc->perform_cpu_sort << std::endl;
    }
  else if (key == 'i')
    {
      plume->drawIsoSurface = !plume->drawIsoSurface;
    }
  else if (key == 'u')
    {
      plume->isoSurface->solid = !plume->isoSurface->solid;
    }
  else if (key == ')')
    {
      plume->isoSurface->increaseMesh();
      plume->oneTime = 0;
    }
  else if (key == '(')
    {
      plume->isoSurface->decreaseMesh();
      plume->oneTime = 0;
    }
  else if (key == 'v')
    {
      plume->color_by_advect_terms = !plume->color_by_advect_terms;
    }
  else if (key == '/')
    {
      if(plume->dc->drawISD == false) {
        plume->dc->drawISD = true;
      } else {
        plume->dc->drawISD = false;
      }
    }
  else if (key == '?')
    {
      plume->writeShadowMapToFile();
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
    
    if(util->advectChoice == 4){
      plume->planeVisual->clickedRangeButton(x,y);

      if(plume->planeVisual->clickedSlider(x,y)){
	plume->planeVisual->move_slider = true;
	plume->dc->rotate_around = false;
      }
      else{
	plume->dc->rotate_around = true;
	plume->planeVisual->move_slider = false;
      }  
    }
    else{
      plume->dc->rotate_around = true;
    }
   
  }
  else{
    plume->dc->rotate_around = false;
    if(util->advectChoice == 4)
      plume->planeVisual->move_slider = false;
  }

  glutPostRedisplay();
}


void motion(int x, int y)
{
  
  float change_y = y - last_y;
  float change_x = x - last_x;
  //std::cout << "last: " << last_y << "   current: " << y << std::endl;
  //std::cout << "CHANGE: " << change_y << std::endl;
  
  //plume->dc->lookUporDown(change_y);
  //plume->dc->setRotateAround(change_x);

  if(util->advectChoice == 4){
    if(plume->planeVisual->move_slider){
      plume->planeVisual->moveSlider(x);
    }
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


void signalHandler(int sig)
{
  // this is set this way so that when we execute the gpuplume code
  // from within the population sampling code, we can actually return
  // from the gpuplume code and have the population sample code
  // cleanup well (fixzing database issues and state).  Unfortunately,
  // this means that we can't ^C the gpuplume code once this gets
  // going, but repeated ^C's will eventually get the pop sample code
  // to exit. 
  
  // eventually, a chained signal handler system might be better but
  // for now, this is simple and functional.

  // in any event, this only ignores the signal IFF the ignoreSignal
  // argument is provided to gpuPlume; otherwise it will exit
  std::cout << "GPUPlume Signal Caught! Cleaning up and exiting." << std::endl;

  if (util)
    {
      if (util->ignoreSignal) 
	{
	  std::cout << "\tignoring signal..." << std::endl;
	  return;
	}
    }

  exit(EXIT_SUCCESS);
}
