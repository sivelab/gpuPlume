#include <iostream>
#include <math.h>

#include "plumeControl.h"
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

// //////////////////////////////////////
// BEGIN -----> QUIC PLUME FORTRAN REFERENCES
// //////////////////////////////////////

#ifdef USE_PLUME_DATA

extern "C"
{
  void readfiles_();
}

// Domain size stored in nx, ny, and nz
extern "C" int __datamodule__nx;
extern "C" int __datamodule__ny;
extern "C" int __datamodule__nz;

extern "C" double __datamodule__dx;
extern "C" double __datamodule__dy;
extern "C" double __datamodule__dz;

// UVW contains the wind field
extern "C" double* __datamodule__u;
extern "C" double* __datamodule__v;
extern "C" double* __datamodule__w;

extern "C" int __datamodule__inumbuild;   // integer number of buildings
extern "C" double* __datamodule__xfo;
extern "C" double* __datamodule__yfo;
extern "C" double* __datamodule__zfo; 
extern "C" double* __datamodule__ht;
extern "C" double* __datamodule__wti;
extern "C" double* __datamodule__lti; 

#endif
// //////////////////////////////////////
// END ----> QUIC PLUME FORTRAN REFERENCES
// //////////////////////////////////////


PlumeControl::PlumeControl(int width, int height, int t){
 
#ifdef USE_PLUME_DATA
  // Call the PLUME code to read in the data files.
  std::cout << "Reading data using PLUME code..." << std::endl;
  readfiles_();

  //QuicPlume data for the domain
  nx = __datamodule__ny; //domain in the x direction
  ny = __datamodule__nz; //domain in the y direction(our orientation is y for up)
  nz = __datamodule__nx; //domain in the z direction

  //nx = (__datamodule__nx - 1) * __datamodule__dx; //domain in the x direction
  //ny = (__datamodule__nz - 1) * __datamodule__dz; //domain in the y direction(our orientation is y for up)
  //nz = (__datamodule__ny - 1) * __datamodule__dy; //domain in the z direction

  std::cout << "QUIC PLUME domain size: " << nx << " (in X) by " 
	    << ny << " (in Y) by " << nz << " (in Z)" << std::endl;

  //QuicPlume data for the windfield
  u = __datamodule__u;
  v = __datamodule__v;
  w = __datamodule__w;

  //QuicPlume data for the buildings
  numBuild = __datamodule__inumbuild;
  xfo = __datamodule__xfo;
  yfo = __datamodule__yfo;
  zfo = __datamodule__zfo;
  ht = __datamodule__ht;
  wti = __datamodule__wti;
  lti = __datamodule__lti;
  
#else
  nx = 60;
  ny = 20;
  nz = 60;
#endif

  //These valuse determine the number of particles
  twidth = width;  theight = height;

  texType = GL_TEXTURE_RECTANGLE_ARB;
  int_format = GL_RGBA32F_ARB;
  int_format_init = GL_RGBA;

  //testcase determines which data set to use for the windfield.
  //The value t is currently passed in from gpuPlume.  When it
  //equals 3, it runs the quicplume data set.  When it equals
  //4 it runs the uniform u-direction windfield.  
  testcase = t;
  
  //This is the time step of the simulation.
  //Set useRealTime to true for real-time simulation
  //or false to use the value given to time_step.
  time_step = 0.0012;
  useRealTime = true;

  odd = true; 
  dump_contents = false;
  emit = true;
  show_particle_visuals = true;

  for(int i = twidth*theight-1; i >= 0; i--)
    indices.push_back(i);

}
void PlumeControl::init(){
  
  pc = new ParticleControl(texType, twidth,theight,nx,ny,nz,u,v,w);

  dc = new DisplayControl(nx,ny,nz, texType);
  dc->initVars(numBuild,xfo,yfo,zfo,ht,wti,lti);

  if(testcase == 3){   
    dc->draw_buildings = true;
    //Creates a point emitter with position(10,10,10) , emitting 10 particles per second
    pe = new PointEmitter(10.0,10.0,10.0, 10.0, &twidth, &theight, &indices, &emit_shader);
  }
  else{
    dc->draw_buildings = false;
    //Creates a sphere emitter with position(30,10,10), emitting 10 pps, with a radius of 4
    pe = new SphereEmitter(30.0, 10.0, 30.0, 30.0, 4.0, &twidth, &theight, &indices, &emit_shader);
  }

  glEnable(texType);
  glGenTextures(8, texid);
  /////////////////////////////
  //Textures used:
  //texid[0] and texid[1] are the double buffered position textures
  //texid[2] can be used to initialize the positions
  //texid[3] is the wind field texture
  //texid[4] is ...
  /////////////////////////////
  setupTextures();

  //
  // set up vertex buffer
  // 
  glGenBuffersARB(1, &vertex_buffer);
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, vertex_buffer);
  glBufferDataARB(GL_ARRAY_BUFFER_ARB, twidth*theight*4*sizeof(GLfloat), 0, GL_STREAM_COPY);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  //Initialize FBO
  initFBO();

  //This shader is used to advect the particles using the windfield
  pc->setupAdvectShader(&time_step, &numInRow);
  //This shader is used to emmit particles
  emit_shader.addShader("Shaders/emitParticle_vp.glsl", GLSLObject::VERTEX_SHADER);
  emit_shader.addShader("Shaders/emitParticle_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  emit_shader.createProgram();

  //Initializes the particle positions in the domain
  //We don't need to do this anymore since we now can initialize data values
  //past 0-1 range in a RGBA32F texture without a shader.
  //pc->initParticlePositions(fbo, texid[2]); 
  CheckErrorsGL("END of init");

  if(useRealTime){
    display_clock = new Timer(true);
    //We need to initialize time 0;
    display_time[0] = display_clock->tic();
  }

} 

void PlumeControl::display(){
  //GLint draw_buffer;
  glGetIntegerv(GL_DRAW_BUFFER, &draw_buffer);
  
  // bind the framebuffer object so we can render to the 2nd texture
  fbo->Bind();

  ////////////////////////////////////////////////////////////
  // Emit Particles
  ///////////////////////////////////////////////////////////
  //record end time
  if(useRealTime)
    display_time[1] = display_clock->tic();
  
  if(emit){
    if(useRealTime){
      float runtime = display_clock->deltas(display_time[0],display_time[1]);
      time_step = runtime;
      //std::cout << runtime << std::endl;
    }
    if(pe->timeToEmit(time_step))
      pe->EmitParticle(fbo, odd);
  }
  //record start time
  if(useRealTime)
    display_time[0] = display_clock->tic();
 
  ////////////////////////////////////////////////////////////
  // Update Particle Positions 
  ////////////////////////////////////////////////////////////
  pc->advect(fbo, odd, texid[4], texid[3], texid[0], texid[1], time_step);

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
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      gluPerspective(60.0, glutGet(GLUT_WINDOW_WIDTH)/float(glutGet(GLUT_WINDOW_HEIGHT)), 1.0, 250.0);
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();		
      

      //plume->displayVisual(vertex_buffer);
      dc->drawVisuals(vertex_buffer, texid[3], numInRow, twidth, theight);
      pe->Draw();

      // If we've chose to display the 3D particle domain, we need to
      // set the projection and modelview matrices back to what is
      // needed for the particle advection step
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      gluOrtho2D(-1, 1, -1, 1);
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
     
      glDisable(texType);
      CheckErrorsGL("END : visualization");
      
      // Finally, swap the front and back buffers to display the
      // particle field to the monitor
      glutSwapBuffers();
    }

}

void PlumeControl::initFBO(void){
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

void PlumeControl::setupTextures()
{
	CheckErrorsGL("BEGIN : Creating textures");

	int sz = 4;
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

	CheckErrorsGL("\tcreated texid[2]...");
		
	// Creates wind field data texture
	pc->initWindTex(texid[3], &numInRow, testcase);
	CheckErrorsGL("\tcreated texid[3], the wind field texture...");

	for (int j=0; j<theight; j++)
		for (int i=0; i<twidth; i++)
		{
			int idx = j*twidth*sz + i*sz;
			data[idx] = data[idx] +  100;
			data[idx+1] = data[idx+1] + 100;
			data[idx+2] = data[idx+2] + 100;
		}
  
	// create the base texture with inital vertex positions
	pc->createTexture(texid[0], int_format, twidth, theight, data);
	CheckErrorsGL("\tcreated texid[0], the position texture...");

	// create a second texture to double buffer the vertex positions
	pc->createTexture(texid[1], int_format, twidth, theight, NULL);
	CheckErrorsGL("\tcreated texid[1], the position texture (double buffer)...");

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

			// normalize
			float mag = sqrt(data[idx]*data[idx] + data[idx+1]*data[idx+1] + data[idx+2]*data[idx+2]);
			data[idx] /= mag;
			data[idx+1] /= mag;
			data[idx+2] /= mag;
		}
	pc->createTexture(texid[4], int_format, twidth, theight, data);
	CheckErrorsGL("\tcreated texid[4], the random number texture...");

  delete [] data;

  CheckErrorsGL("END : Creating textures");
}

