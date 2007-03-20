#include <iostream>
#include <math.h>
#include "particleControl.h"

#ifdef WIN32
float randomVal() { return (float)(rand()/(float)RAND_MAX); }
#else
float randomVal(){ return drand48();}
#endif


ParticleControl::ParticleControl(GLenum type,int width,int height,
				 int x, int y, int z,
				 double* u, double* v, double* w){

  texType = type;
  twidth = width;
  theight = height;
  nx = x;
  ny = y;
  nz = z;
  u_quicPlumeData = v;
  v_quicPlumeData = w;
  w_quicPlumeData = u;
}
void ParticleControl::setupAdvectShader(float* time_step, int* numInRow){

  //This shader is used to move the particles
  pass1_shader.addShader("Shaders/plumeAdvect_vp.glsl", GLSLObject::VERTEX_SHADER);
  pass1_shader.addShader("Shaders/plumeAdvect_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  pass1_shader.createProgram();

  // Get location of the sampler uniform
  uniform_postex = pass1_shader.createUniform("pos_texunit");
  uniform_wind = pass1_shader.createUniform("wind_texunit");
  uniform_randomTexture = pass1_shader.createUniform("random_texunit");
  uniform_timeStep = pass1_shader.createUniform("time_step");
  GLint unx = pass1_shader.createUniform("nx");
  GLint uny = pass1_shader.createUniform("ny");
  GLint unz = pass1_shader.createUniform("nz");
  GLint uNumInRow = pass1_shader.createUniform("numInRow");

  pass1_shader.activate();

  glUniform1fARB(unx, nx);
  glUniform1fARB(uny, ny);
  glUniform1fARB(unz, nz);
  glUniform1fARB(uniform_timeStep, *time_step);
  float numR= *numInRow;
  glUniform1fARB(uNumInRow, numR);

  pass1_shader.deactivate();

}
void ParticleControl::advect(FramebufferObject* fbo, bool odd, GLuint texid4, GLuint texid3,
			     GLuint texid0, GLuint texid1, float time_step)
{
  if (odd)
    glDrawBuffer(GL_COLOR_ATTACHMENT1_EXT);
  else 
    glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glViewport(0, 0, twidth, theight);
  
  
  glEnable(texType);
  pass1_shader.activate();
  glUniform1fARB(uniform_timeStep, time_step);

  // Bind the random data field to TEXTURE UNIT 2
  glActiveTexture(GL_TEXTURE2);
  glBindTexture(texType, texid4);
  glUniform1iARB(uniform_randomTexture, 2);

  // wind field can be stored here
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(texType, texid3);
  glUniform1iARB(uniform_wind, 1);
    
  glActiveTexture(GL_TEXTURE0);
  glUniform1iARB(uniform_postex, 0);

  if (odd)
    glBindTexture(texType, texid0);  // read from texture 0
  else 
    glBindTexture(texType, texid1);  // read from texture 1

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

}

void ParticleControl::getDomain(int* x, int* y, int* z){
  *x = nx;
  *y = ny;
  *z = nz;
}

void ParticleControl::dumpContents(){
   buffer_mem = new GLfloat[ twidth * theight * 4 ];
   glReadPixels(0, 0, twidth, theight, GL_RGBA, GL_FLOAT, buffer_mem);
      for (int j=0; j<theight; j++)
	for (int i=0; i<twidth; i++){
	  
	    int idx = j*twidth*4 + i*4;
	    std::cout << buffer_mem[idx] << " ";
	    std::cout << buffer_mem[idx+1] << " ";
	    std::cout << buffer_mem[idx+2] << " ";
	    std::cout << buffer_mem[idx+3] << std::endl;
	}
      delete [] buffer_mem;
}

void ParticleControl::createTexture(GLuint texId, GLenum format, int w, int h, GLfloat* data){
  glBindTexture(texType, texId);
  glTexParameteri(texType, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(texType, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(texType, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(texType, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  glTexImage2D(texType, 0, format, w, h, 0, GL_RGBA, GL_FLOAT, data);

}

void ParticleControl::createWrappedTexture(GLuint texId, GLenum format, int w, int h, GLfloat* data)
{
  glBindTexture(texType, texId);
  glTexParameteri(texType, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(texType, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(texType, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(texType, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexImage2D(texType, 0, format, w, h, 0, GL_RGBA, GL_FLOAT, data);
}

void ParticleControl::initWindTex(GLuint texId, int* numInRow, int dataSet){
  // Create wind data texture
  data3d = new wind[nx*ny*nz];
  switch(dataSet){

      case 1:
	 test1();
	 break;
      case 2:
	randomWindField();
	break;
      case 3:
	quicPlumeWindField();
	break;
      case 4:
	uniformUWindField();
	break;

  }

  /////////////////////////////////////////////////////////
  //Calculate width and height for wind texture
  //
  //This tries to minimize the width and height values
  //to try and fit the wind field into a 2D texture without
  //wasting too much space.  
  /////////////////////////////////////////////////////////
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

  ////////////////////////////////////////////////////////
  //Convert this to data array for a texture
  //
  //This will directly put the 3D data into an array
  //that is used to make the 2D texture.
  ///////////////////////////////////////////////////////
  (*numInRow) = (width - (width % nz))/nz;
  //std::cout << width << " " << *numInRow << std::endl;

  int qi, qj, qk;
  int p2idx = 0, texidx = 0;
  int row = 0;
  
  GLfloat *dataTwo = new GLfloat[ width * height * 4 ];
  for (qk=0; qk<ny; qk++) 
    for (qi=0; qi<nx; qi++)
      for (qj=0; qj<nz; qj++)
	{
	  p2idx = qk*nz*nx + qi*nz + qj;
	    
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
void ParticleControl::randomWindField(){
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

void ParticleControl::quicPlumeWindField()
{
  //#ifdef USE_PLUME_DATA
  for(int k = 0; k < ny; k++){   
    for(int i = 0; i < nx; i++){
      for(int j = 0; j < nz; j++){
	int p2idx = k*(nx)*(nz) + i*(nz) + j;
	int idx = k*nx*nz + j*nz + i;
	data3d[p2idx].u = u_quicPlumeData[idx+1];
	data3d[p2idx].v = v_quicPlumeData[idx+1];
	data3d[p2idx].w = w_quicPlumeData[idx+1];
      }
    }
  }
  //#endif
}

void ParticleControl::uniformUWindField(){
  for(int k = 0; k < ny; k++){   
    for(int i = 0; i < nx; i++){
      for(int j = 0; j < nz; j++){
	int p2idx = k*nx*nz + i*nz + j;
	data3d[p2idx].u = 0.0;
	data3d[p2idx].v = 0.0;
	data3d[p2idx].w = 1.0;
      }
    }
  }
}
void ParticleControl::initParticlePositions(FramebufferObject* fbo, GLuint texId){
  
   //This shader is used to initialize the particle positions
  init_shader.addShader("Shaders/initialize_vp.glsl", GLSLObject::VERTEX_SHADER);
  init_shader.addShader("Shaders/initialize_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  init_shader.createProgram();

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
