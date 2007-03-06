
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


PlumeControl::PlumeControl(int w, int h, int t){
  //These valuse determine the number of particles
  twidth = w;  theight = h;

  testcase = t;
  time_step = 0.0012;
  odd = true;
  dump_contents = false;
  emit = false;
  show_particle_visuals = true;

  texType = GL_TEXTURE_RECTANGLE_ARB;
  int_format = GL_RGBA32F_ARB;
  int_format_init = GL_RGBA;

  for(int i = twidth*theight-1; i >= 0; i--)
    indices.push_back(i);

}
void PlumeControl::init(){
 
  pc = new ParticleControl(texType, twidth,theight);
  int nx, ny, nz;
  pc->getDomain(&nx,&ny,&nz);
  dc = new DisplayControl(nx,ny,nz, texType);
  //pe = new ParticleEmitter(10.0,10.0, 10.0, 10.0, &twidth, &theight, &indices, &emit_shader);
  if(testcase == 3){
    dc->draw_buildings = true;
    pe = new ParticleEmitter(10.0,10.0, 10.0, 10.0, &twidth, &theight, &indices, &emit_shader);
  }
  else{
    dc->draw_buildings = false;
    pe = new ParticleEmitter(30.0, 10.0, 30.0, 30.0, &twidth, &theight, &indices, &emit_shader);
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

  display_clock = new Timer(true);

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
  display_time[1] = display_clock->tic();
  
  if(emit){
    float runtime = display_clock->deltas(display_time[0],display_time[1]);
    time_step = runtime;
    //std::cout << runtime << std::endl;
    if(pe->timeToEmit(time_step))
      pe->EmitParticle(fbo, odd);
  }
  //record start time
  display_time[0] = display_clock->tic();
  //emit = false;
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

      //plume->displayVisual(vertex_buffer);
      dc->drawVisuals(vertex_buffer, texid[3], numInRow, twidth, theight);

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

void PlumeControl::setupTextures(){
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

  // Creates wind field data texture
  pc->initWindTex(texid[3], &numInRow, testcase);
  
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
}

