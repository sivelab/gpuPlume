
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
    pe = new ParticleEmitter(30.0, 10.0, 30.0, 10.0, &twidth, &theight, &indices, &emit_shader);
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

  //This shader is used to advect the particles using the windfield
  pc->setupAdvectShader(&time_step, &numInRow);
  //This shader is used to emmit particles
  emit_shader.addShader("Shaders/emitParticle_vp.glsl", GLSLObject::VERTEX_SHADER);
  emit_shader.addShader("Shaders/emitParticle_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  emit_shader.createProgram();
}
void PlumeControl::injectParticles(FramebufferObject* fbo, bool odd){
  if(pe->timeToEmit(time_step))
    pe->EmitParticle(fbo, odd);
}

void PlumeControl::advectParticles(FramebufferObject* fbo, bool odd){
  pc->advect(fbo, odd, texid[4], texid[3], texid[0], texid[1]);
}

void PlumeControl::displayVisual(GLuint vertex_buffer){
  dc->drawVisuals(vertex_buffer, texid[3], numInRow, twidth, theight);
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

