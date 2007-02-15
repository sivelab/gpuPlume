#include <iostream>
#include <math.h>
#include "particleControl.h"
using namespace std;

#ifdef WIN32
float randomVal() { return (float)(rand()/(float)RAND_MAX); }
#else
float randomVal(){ return drand48();}
#endif

// //////////////////////////////////////
// BEGIN -----> QUIC PLUME FORTRAN REFERENCES
// //////////////////////////////////////

#define USE_PLUME_DATA
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

#endif
// //////////////////////////////////////
// END ----> QUIC PLUME FORTRAN REFERENCES
// //////////////////////////////////////


ParticleControl::ParticleControl(GLenum type){

  texType = type;

#ifdef USE_PLUME_DATA
  // Call the PLUME code to read in the data files.
  std::cout << "Reading data using PLUME code..." << std::endl;
  readfiles_();

  nx = __datamodule__nx; //domain in the x direction
  ny = __datamodule__nz; //domain in the y direction(our orientation is y for up)
  nz = __datamodule__ny; //domain in the z direction

  // nx = (__datamodule__nx - 1) * __datamodule__dx; //domain in the x direction
  // ny = (__datamodule__nz - 1) * __datamodule__dz; //domain in the y direction(our orientation is y for up)
  // nz = (__datamodule__ny - 1) * __datamodule__dy; //domain in the z direction

  std::cout << "QUIC PLUME domain size: " << nx << " (in X) by " << ny << " (in Y) by " << nz << " (in Z)" << std::endl;

#else
  nx = 60;
  ny = 20;
  nz = 60;
#endif

}

void ParticleControl::getDomain(int* x, int* y, int* z){
  *x = nx;
  *y = ny;
  *z = nz;
}
void ParticleControl::dumpContents(int w, int h){
   buffer_mem = new GLfloat[ w * h * 4 ];
   glReadPixels(0, 0, w, h, GL_RGBA, GL_FLOAT, buffer_mem);
      for (int j=0; j<w; j++)
	for (int i=0; i<h; i++){
	  
	    int idx = j*w*4 + i*4;
	    std::cout << buffer_mem[idx] << " ";
	    std::cout << buffer_mem[idx+1] << " ";
	    std::cout << buffer_mem[idx+2] << " ";
	    std::cout << buffer_mem[idx+3] << std::endl;
	}
      delete [] buffer_mem;
}
void ParticleControl::initParticlePositions(FramebufferObject* fbo, int twidth, int theight, GLSLObject init_shader, GLuint texId){
  GLint vp[4];
  glGetIntegerv(GL_VIEWPORT, vp);
  GLint draw_buffer;
  glGetIntegerv(GL_DRAW_BUFFER, &draw_buffer);
  fbo->Bind();
  glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
  glViewport(0, 0, twidth, theight);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  GLuint texture = init_shader.createUniform("texture");
  GLuint u_nx = init_shader.createUniform("nx");
  GLuint u_ny = init_shader.createUniform("ny");
  GLuint u_nz = init_shader.createUniform("nz");
  glBindTexture(texType, texId);
  init_shader.activate();
  glUniform1iARB(texture, 0);
  glUniform1fARB(u_nx, nx);
  glUniform1fARB(u_ny, ny);
  glUniform1fARB(u_nz, nz);

  glBegin(GL_QUADS);
  {
    glTexCoord2f(0, 0);            glVertex3f(-1, -1, -0.5f);
    glTexCoord2f(twidth, 0);       glVertex3f( 1, -1, -0.5f);
    glTexCoord2f(twidth, theight); glVertex3f( 1,  1, -0.5f);
    glTexCoord2f(0, theight);      glVertex3f(-1,  1, -0.5f);
  }
  glEnd();

  init_shader.deactivate();
  glDisable(texType);
  FramebufferObject::Disable();
  glViewport(vp[0], vp[1], vp[2], vp[3]);
  glDrawBuffer(draw_buffer);

}

void ParticleControl::createTexture(GLuint texId, GLenum format, int w, int h, GLfloat* data){
  glBindTexture(texType, texId);
  glTexParameteri(texType, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(texType, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(texType, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(texType, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  glTexImage2D(texType, 0, format, w, h, 0, GL_RGBA, GL_FLOAT, data);

}
void ParticleControl::initWindTex(GLuint texId, int* numInRow){
  // Create wind data texture
  data3d = new wind[nx*ny*nz];
  //test1();

#ifdef USE_PLUME_DATA
  test3();
#else
  test2();
#endif
  /////////////////////////////////////////////////////////
  //Calculate width and height for wind texture
  //
  //This tries to minimize the width and height values
  //to try and fit the wind field into a 2D texture without
  //wasting too much space.  
  int total = nx*ny*nz;
  int width = (int)sqrt((float)total);
  
  int scaler;
  if(nx > nz) scaler = nx;
  else scaler = nz;

  width = width - (width%scaler);
  
  bool done = false;
  while(!done){ 
    int num = width/scaler;
    if((num*num) >= ny){
      done = true;
    }
    else{
      width = width+scaler;
    }
  }
  if(width%2 != 0) width++;
  int height = width;
  /////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////
  //Convert this to data array for a texture
  //
  //This will directly put the 3D data into an array
  //that is used to make the 2D texture.
  (*numInRow) = (width - (width % nz))/nz;
  int qi, qj, qk;
  int p2idx = 0, texidx = 0;
  int row = 0;
  
  GLfloat *dataTwo = new GLfloat[ width * height * 4 ];
  for (qk=0; qk<ny; qk++) 
    for (qi=0; qi<nx; qi++)
      for (qj=0; qj<nz; qj++)
	{
	  p2idx = qk*nz*nx + qj*nz + qi;
	    
	  row = qk / (*numInRow);
	  texidx = row * width * nx * 4 +
	  qi * width * 4 +
	  qk % (*numInRow) * nz * 4 +
	  qj * 4;
	  
	  dataTwo[texidx] = data3d[p2idx].u;
	  dataTwo[texidx+1] = data3d[p2idx].v;
	  dataTwo[texidx+2] = data3d[p2idx].w;
	  dataTwo[texidx+3] = 1.0;	      			      
        }

  createTexture(texId, GL_RGBA, width, height, dataTwo);
  
  delete [] dataTwo;
  
}
void ParticleControl::test1(){
  for(int k = 0; k < ny; k++){   
    for(int i = 0; i < nx; i++){
      for(int j = 0; j < nz; j++){
	int p2idx = k*nx*nz + i*nz + j;
  	if(i%2 == 0){
	  data3d[p2idx].u = 0;
	  data3d[p2idx].v = 1.0;
	  data3d[p2idx].w = 0;
	 }
	else {
	  data3d[p2idx].u = 0;
	  data3d[p2idx].v = 0;
	  data3d[p2idx].w = 0;
	  }
      }
    }
  }
}
//Creates a random value wind field.
void ParticleControl::test2(){
  for(int k = 0; k < ny; k++){   
    for(int i = 0; i < nx; i++){
      for(int j = 0; j < nz; j++){
	int p2idx = k*nx*nz + i*nz + j;
	//
	// Randomly generate a vector within the unit sphere
	// 	
	data3d[p2idx].u = randomVal();
	data3d[p2idx].v = randomVal();
	data3d[p2idx].w = randomVal();
      }
    }
  }
}
//Uses the QUIC_PLUME data for the wind field.
#ifdef USE_PLUME_DATA
void ParticleControl::test3(){
  for(int k = 0; k < ny; k++){   
    for(int i = 0; i < nx; i++){
      for(int j = 0; j < nz; j++){
	int p2idx = k*nx*nz + i*nz + j;
	data3d[p2idx].u = __datamodule__u[p2idx];
	data3d[p2idx].v = __datamodule__w[p2idx];
	data3d[p2idx].w = __datamodule__v[p2idx];
      }
    }
  }
}
#endif
