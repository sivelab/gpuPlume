#include <iostream>
#include <math.h>
#include "particleControl.h"
#include "Random.h"

ParticleControl::ParticleControl(GLenum type,int width,int height,
				 int x, int y, int z,
				 double* u, double* v, double* w){

  texType = type;
  twidth = width;
  theight = height;
  nx = x;
  ny = y;
  nz = z;
  u_quicPlumeData = u;
  v_quicPlumeData = v;
  w_quicPlumeData = w;

  outputPrime = false;

}
void ParticleControl::setUstarAndSigmas(float u){
  ustar = u;
  sigU = 2.0*ustar;
  sigV = 2.0*ustar;
  sigW = 1.3*ustar;

}

void ParticleControl::setupPrimeShader( int* numInRow){  //Included argument -- Balli(04/12/07)
  //This shader is used to update prime values
  prime_shader.addShader("Shaders/updatePrime_vp.glsl", GLSLObject::VERTEX_SHADER);
  prime_shader.addShader("Shaders/updatePrime_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  prime_shader.createProgram();

  uniform_pos = prime_shader.createUniform("pos"); //Balli (04/12/07)- position texture is need for epsilon calculations  
  uniform_prime = prime_shader.createUniform("primePrev");
  uniform_random = prime_shader.createUniform("random");
  uniform_windTex = prime_shader.createUniform("wind");
  uniform_lambda = prime_shader.createUniform("lambda");
  uniform_dt = prime_shader.createUniform("time_step");
  uniform_randomTexCoordOffset = prime_shader.createUniform("random_texCoordOffset");
  uniform_randomTexWidth = prime_shader.createUniform("random_texWidth");
  uniform_randomTexHeight = prime_shader.createUniform("random_texHeight");
  
  // Following is copy-paste from setupadvect shader
  //we need all of the following in prime shader too. --Balli(04/12/07)
  GLint unx = prime_shader.createUniform("nx");
  GLint uny = prime_shader.createUniform("ny");
  GLint unz = prime_shader.createUniform("nz");
  GLint uNumInRow = prime_shader.createUniform("numInRow");

  prime_shader.activate();

  glUniform1iARB(unx, nx);
  glUniform1iARB(uny, ny);
  glUniform1iARB(unz, nz);
  float numR= *numInRow;
  glUniform1fARB(uNumInRow, numR);

  prime_shader.deactivate();
  //End --Balli(04/12/07)
}
void ParticleControl::updatePrime(FramebufferObject* fbo, bool odd, GLuint positions0,GLuint positions1,GLuint prime0, 
				  GLuint prime1, GLuint windField, GLuint randomValues, 
				  GLuint lambda, float time_step){  // included two more argument for position--Balli(04/12/07)


  //Prints out the previous prime values
  if(outputPrime)
    printPrime(odd, true);

  if (odd)
    glDrawBuffer(GL_COLOR_ATTACHMENT3_EXT);
  else 
    glDrawBuffer(GL_COLOR_ATTACHMENT2_EXT);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glViewport(0, 0, twidth, theight);
  
  
  glEnable(texType);
  prime_shader.activate();
  glUniform1fARB(uniform_dt, time_step);

  // generate and set the random texture coordinate offset
  // 
  if (texType == GL_TEXTURE_RECTANGLE_ARB)
    {
      // texture coordinates will range from 0 to W in width and 0 to
      // H in height, so generate random value in this range
      float f1 = Random::uniform() * twidth;
      float f2 = Random::uniform() * theight;
      glUniform2fARB(uniform_randomTexCoordOffset, f1, f2);
    }
  else 
    {
      // texture coordinates will range from 0 to 1, so generate random value in this range
      glUniform2fARB(uniform_randomTexCoordOffset, Random::uniform(), Random::uniform());
    }

  // set the size of the texture width and height for the shader to use
  glUniform1iARB(uniform_randomTexWidth, twidth);
  glUniform1iARB(uniform_randomTexHeight, theight);

  //Bind the position texture to TEXTURE UNIT 4-- Balli (04/12/07)
  glActiveTexture(GL_TEXTURE4);
  glUniform1iARB(uniform_pos, 4);
  if (odd)
    glBindTexture(texType, positions0);  // read from texture 0
  else 
    glBindTexture(texType, positions1);  // read from texture 1

  //Bind the lambda texture to TEXTURE UNIT 3
  glActiveTexture(GL_TEXTURE3);
  glBindTexture(texType, lambda);
  glUniform1iARB(uniform_lambda, 3);

  // Bind the random data field to TEXTURE UNIT 2
  glActiveTexture(GL_TEXTURE2);
  glBindTexture(texType, randomValues);
  glUniform1iARB(uniform_random, 2);

  // wind field can be stored here in TEXTURE UNIT 1
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(texType, windField);
  glUniform1iARB(uniform_windTex, 1);
    
  glActiveTexture(GL_TEXTURE0);
  glUniform1iARB(uniform_prime, 0);

  if (odd)
    glBindTexture(texType, prime0);  // read from prime texture 0
  else 
    glBindTexture(texType, prime1);  // read from prime texture 1


  glBegin(GL_QUADS);
  {
    glTexCoord2f(0, 0);            glVertex3f(-1, -1, -0.5f);
    glTexCoord2f(twidth, 0);       glVertex3f( 1, -1, -0.5f);
    glTexCoord2f(twidth, theight); glVertex3f( 1,  1, -0.5f);
    glTexCoord2f(0, theight);      glVertex3f(-1,  1, -0.5f);
  }
  glEnd();
  
  prime_shader.deactivate();
  
  //Prints out the updated prime values
  if(outputPrime)
    printPrime(odd,false);  
 
  glBindTexture(texType, 0);

}
void ParticleControl::printPrime(bool odd, bool prev){
  glGetIntegerv(GL_READ_BUFFER, &currentbuffer);
  
  if(prev){
    std::cout << "Previous Prime Values" << std::endl;
    if(odd)
      glReadBuffer(GL_COLOR_ATTACHMENT2_EXT);
    else
      glReadBuffer(GL_COLOR_ATTACHMENT3_EXT);
  }
  else{

    std::cout << "Updated Prime Values" << std::endl;
    if(odd)
      glReadBuffer(GL_COLOR_ATTACHMENT3_EXT);
    else
      glReadBuffer(GL_COLOR_ATTACHMENT2_EXT);

    outputPrime = false;

  }

  dumpContents();
  
 
  glReadBuffer(currentbuffer);


}
void ParticleControl::setupAdvectShader(float* time_step, int* numInRow, float life_time){

  //This shader is used to move the particles
  pass1_shader.addShader("Shaders/plumeAdvect_vp.glsl", GLSLObject::VERTEX_SHADER);
  pass1_shader.addShader("Shaders/plumeAdvect_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  pass1_shader.createProgram();

  // Get location of the sampler uniform
  uniform_postex = pass1_shader.createUniform("pos_texunit");
  uniform_wind = pass1_shader.createUniform("wind_texunit");
  uniform_timeStep = pass1_shader.createUniform("time_step");
  uniform_primePrev = pass1_shader.createUniform("primePrev");
  uniform_primeCurr = pass1_shader.createUniform("primeCurr");

  GLint ulifeTime = pass1_shader.createUniform("life_time");
  GLint unx = pass1_shader.createUniform("nx");
  GLint uny = pass1_shader.createUniform("ny");
  GLint unz = pass1_shader.createUniform("nz");
  GLint uNumInRow = pass1_shader.createUniform("numInRow");

  pass1_shader.activate();

  glUniform1fARB(ulifeTime, life_time);
  glUniform1iARB(unx, nx);
  glUniform1iARB(uny, ny);
  glUniform1iARB(unz, nz);
  glUniform1fARB(uniform_timeStep, *time_step);
  float numR= *numInRow;
  glUniform1fARB(uNumInRow, numR);

  pass1_shader.deactivate();

}
void ParticleControl::advect(FramebufferObject* fbo, bool odd, 
			     GLuint windField, GLuint positions0, GLuint positions1, 
			     GLuint prime0, GLuint prime1, float time_step)
{

  //fbo->Bind();

  if (odd)
    glDrawBuffer(GL_COLOR_ATTACHMENT1_EXT);
  else 
    glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glViewport(0, 0, twidth, theight);
  
  
  glEnable(texType);
  pass1_shader.activate();
  glUniform1fARB(uniform_timeStep, time_step);

  //Bind current prime values to TEXTURE UNIT 4
  glActiveTexture(GL_TEXTURE3);
  glUniform1iARB(uniform_primeCurr, 3);
  if(odd)
    glBindTexture(texType, prime1);
  else
    glBindTexture(texType, prime0);

  //Bind previous prime values to TEXTURE UNIT 3
  glActiveTexture(GL_TEXTURE2);
  glUniform1iARB(uniform_primePrev, 2);
  if(odd)
    glBindTexture(texType, prime0);
  else
    glBindTexture(texType, prime1);


  // wind field can be stored here
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(texType, windField);
  glUniform1iARB(uniform_wind, 1);
    
  glActiveTexture(GL_TEXTURE0);
  glUniform1iARB(uniform_postex, 0);

  if (odd)
    glBindTexture(texType, positions0);  // read from texture 0
  else 
    glBindTexture(texType, positions1);  // read from texture 1
 
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
   std::cout << "IDX  X     Y     Z" << std::endl;
   int pn =0;
      for (int j=0; j<theight; j++)
	for (int i=0; i<twidth; i++){
	  
	    int idx = j*twidth*4 + i*4;
	    std::cout << pn << " ";
	    std::cout << buffer_mem[idx] << " ";
	    std::cout << buffer_mem[idx+1] << " ";
	    std::cout << buffer_mem[idx+2] << " ";
	    std::cout << buffer_mem[idx+3] << std::endl;
	    pn++;
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
  //
  // NOTE: GL_REPEAT IS NOT SUPPORTED FOR GL_TEXTURE_RECTANGLE_ARB
  // Thus, we must use the normal clamp to edge texture, but do the work in the shader to correct the random
  // offset.
  //
  glTexParameteri(texType, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(texType, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexImage2D(texType, 0, format, w, h, 0, GL_RGBA, GL_FLOAT, data);
}

void ParticleControl::initWindTex(GLuint windField, GLuint lambda, int* numInRow,
				  int dataSet){
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
      case 5:
        gravity();
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
  if(ny > nx) scaler = ny;
  else scaler = nx;

  width = width - (width%scaler);
  
  bool done = false;
  while(!done){ 
    int num = width/scaler;
    if((num*num) >= nz){
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
  (*numInRow) = (width - (width % nx))/nx;
  //std::cout << width << " " << *numInRow << std::endl;

  int qi, qj, qk;
  int p2idx = 0, texidx = 0;
  int row = 0;
  
  GLfloat *dataTwo = new GLfloat[ width * height * 4 ];
  for (qk=0; qk<nz; qk++) 
    for (qi=0; qi<ny; qi++)
      for (qj=0; qj<nx; qj++)
	{
	  p2idx = qk*ny*nx + qi*nx + qj;
	    
	  row = qk / (*numInRow);
	  texidx = row * width * ny * 4 +
	  qi * width * 4 +
	  qk % (*numInRow) * nx * 4 +
	  qj * 4;
	  
	  dataTwo[texidx] = data3d[p2idx].u;
	  dataTwo[texidx+1] = data3d[p2idx].v;
	  dataTwo[texidx+2] = data3d[p2idx].w;	  
	  dataTwo[texidx+3] = (0.5*5.7)*(ustar*ustar*ustar)/(0.4*(qk+1));//This value is the '0.5*CoEps' value	

        }

  createTexture(windField, GL_RGBA32F_ARB, width, height, dataTwo);

  //Create lambda texture. --Balli (04/12/07)
  float tau11=sigU*sigU;
  float tau22=sigV*sigV;
  float tau33=sigW*sigW;
  float tau13=ustar*ustar;
  float tauDetInv=1/((tau11*tau22*tau33)-(tau13*tau13*tau22));

   for (qk=0; qk<nz; qk++) 
    for (qi=0; qi<ny; qi++)
      for (qj=0; qj<nx; qj++)
	{
	  p2idx = qk*ny*nx + qi*nx + qj;
	    
	  row = qk / (*numInRow);
	  texidx = row * width * ny * 4 +
	  qi * width * 4 +
	  qk % (*numInRow) * nx * 4 +
	  qj * 4;
	  
	  dataTwo[texidx]   = tauDetInv*(tau22*tau33);              //Lam11
	  dataTwo[texidx+1] = tauDetInv*(tau11*tau33-tau13*tau13);//Lam22
	  dataTwo[texidx+2] = tauDetInv*(tau11*tau22);	          //Lam33
	  dataTwo[texidx+3] = tauDetInv*(-tau13*tau22);           //Lam13
        }

  createTexture(lambda, GL_RGBA32F_ARB, width,height, dataTwo);
  //Lambda Texture Ends-- Balli (04/12/07)
  delete [] dataTwo;
  
}

void ParticleControl::test1(){
  for(int k = 0; k < nz; k++){   
    for(int i = 0; i < ny; i++){
      for(int j = 0; j < nx; j++){
	int p2idx = k*nx*ny + i*nx + j;
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
  for(int k = 0; k < nz; k++){   
    for(int i = 0; i < ny; i++){
      for(int j = 0; j < nx; j++){
	int p2idx = k*nx*ny + i*nx + j;
	//
	// Not currently random
	// 	
	data3d[p2idx].u = 0.0;//randVal();
	data3d[p2idx].v = 1.0;//randVal();
	data3d[p2idx].w = 0.0;//randVal();
      }
    }
  }
}
//Uses the QUIC_PLUME data for the wind field.

void ParticleControl::quicPlumeWindField()
{
  //#ifdef USE_PLUME_DATA
  for(int k = 0; k < nz; k++){   
    for(int i = 0; i < ny; i++){
      for(int j = 0; j < nx; j++){
	int p2idx = k*(nx)*(ny) + i*(nx) + j;
	int idx = k*nx*ny + i*nx + j;
	data3d[p2idx].u = u_quicPlumeData[idx+1];
	data3d[p2idx].v = v_quicPlumeData[idx+1];
	data3d[p2idx].w = w_quicPlumeData[idx+1];
      }
    }
  }
  //#endif
}

void ParticleControl::uniformUWindField(){
  for(int k = 0; k < nz; k++){   
    for(int i = 0; i < ny; i++){
      for(int j = 0; j < nx; j++){
	int p2idx = k*nx*ny + i*nx + j;
	data3d[p2idx].u = 1.0;
	data3d[p2idx].v = 0.0;
	data3d[p2idx].w = 0.0;
      }
    }
  }
}
void ParticleControl::gravity(){
   for(int k = 0; k < nz; k++){   
    for(int i = 0; i < ny; i++){
      for(int j = 0; j < nx; j++){
	int p2idx = k*nx*ny + i*nx + j;
	data3d[p2idx].u = 0.0;
	data3d[p2idx].v = 0.0;
	data3d[p2idx].w = -9.81;
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
