#include <stdio.h>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "particleControl.h"
#include "Random.h"
#include "glErrorUtil.h"

ParticleControl::ParticleControl(GLenum type,int width,int height,
				 int x, int y, int z,
				 float c_dx, float c_dy, float c_dz){

  texType = type;
  twidth = width;
  theight = height;
  nx = x;
  ny = y;
  nz = z;

  outputPrime = false;
  alreadyOpen = false;
  for(int i=0; i <4; i++){
    tauMax[i] = 0;
    tauMin[i] = 10;
    windMax[i] = 0;
    windMin[i] = 10;
  }

  //cell_dx = 0.5;
  //cell_dy = 0.5;
  //cell_dz = 0.5;

  nzdz = (int)(nz*(1.0/c_dz));
  nydy = (int)(ny*(1.0/c_dy));
  nxdx = (int)(nx*(1.0/c_dx));

 
  cell_dx = c_dx;
  cell_dy = c_dy;
  cell_dz = c_dz;


}
void ParticleControl::setBuildingParameters(int nB,int* ns,float* x,float* y,float* z,
					   float* h,float* w,float* l){
  numBuild = nB;
  numSides = ns;
  xfo = x;
  yfo = y;
  zfo = z;
  ht = h;
  wti = w;
  lti = l; 
}
void ParticleControl::setQuicFilesPath(std::string path){
  quicFilesPath = path;
}
void ParticleControl::setUstarAndSigmas(float u){
  ustar = u;
  sigU = float(2.0*ustar);
  sigV = float(2.0*ustar);
  sigW = float(1.3*ustar);
}
void ParticleControl::setRandomTexCoords(){
  t1 = Random::uniform() * twidth;
  t2 = Random::uniform() * theight;
}
void ParticleControl::setupMultipleBuildingsShader(float life_time, int shader){
  multipleBuildings_shader.addShader("Shaders/multipleBuildingsAdvect_vp.glsl", GLSLObject::VERTEX_SHADER);
  if(shader == 0)
    multipleBuildings_shader.addShader("Shaders/multipleBuildingsAdvect_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  else
    multipleBuildings_shader.addShader("Shaders/multipleBuildingsAdvectVariedinU_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  multipleBuildings_shader.createProgram();

  // Get location of the sampler uniform
  uniform_postex = multipleBuildings_shader.createUniform("pos_texunit");
  uniform_wind = multipleBuildings_shader.createUniform("wind_texunit");
  uniform_timeStep = multipleBuildings_shader.createUniform("time_step");
  uniform_primePrev = multipleBuildings_shader.createUniform("primePrev");
  uniform_random = multipleBuildings_shader.createUniform("random");
  uniform_lambda = multipleBuildings_shader.createUniform("lambda");
  uniform_tau_dz = multipleBuildings_shader.createUniform("tau_dz");
  uniform_duvw_dz = multipleBuildings_shader.createUniform("duvw_dz");
  uniform_randomTexCoordOffset = multipleBuildings_shader.createUniform("random_texCoordOffset");
  uniform_randomTexWidth = multipleBuildings_shader.createUniform("random_texWidth");
  uniform_randomTexHeight = multipleBuildings_shader.createUniform("random_texHeight");
  uniform_buildings = multipleBuildings_shader.createUniform("buildings");
  uniform_cellType = multipleBuildings_shader.createUniform("cellType");
  uniform_colorAdvectTerms = multipleBuildings_shader.createUniform("color_advect_terms");
  
  GLint ulifeTime = multipleBuildings_shader.createUniform("life_time");
  GLint unx = multipleBuildings_shader.createUniform("nx");
  GLint uny = multipleBuildings_shader.createUniform("ny");
  GLint unz = multipleBuildings_shader.createUniform("nz");
  GLint unxdx = multipleBuildings_shader.createUniform("nxdx");
  GLint unydy = multipleBuildings_shader.createUniform("nydy");
  GLint unzdz = multipleBuildings_shader.createUniform("nzdz");
  GLint uNumInRow = multipleBuildings_shader.createUniform("numInRow");
  GLint udx = multipleBuildings_shader.createUniform("dx");
  GLint udy = multipleBuildings_shader.createUniform("dy");
  GLint udz = multipleBuildings_shader.createUniform("dz");
    
  multipleBuildings_shader.activate();

  glUniform1fARB(ulifeTime, life_time);
  glUniform1iARB(unx, nx);
  glUniform1iARB(uny, ny);
  glUniform1iARB(unz, nz);
  glUniform1fARB(udx, cell_dx);
  glUniform1fARB(udy, cell_dy);
  glUniform1fARB(udz, cell_dz);
  glUniform1iARB(unxdx, nxdx);
  glUniform1iARB(unydy, nydy);
  glUniform1iARB(unzdz, nzdz);
  glUniform1iARB(uNumInRow, numInRow);

  multipleBuildings_shader.deactivate();
}
void ParticleControl::multipleBuildingsAdvect(bool odd, GLuint windField, GLuint positions0, 
			     GLuint positions1, GLuint prime0, GLuint prime1, 
			     GLuint randomValues,GLuint lambda, GLuint tau_dz, 
			     GLuint duvw_dz, float time_step, GLuint buildings,
					      GLuint cellType, GLuint advect_terms)
{
  //Prints out the previous prime values
  if(outputPrime)
    printPrime(odd, true);


  //If coloring particles using advect_terms, advect_terms will not be NULL here.
  //Will need to render to a third texture.
  if(advect_terms != 0){

    if(odd){
      GLenum buffers[] = {GL_COLOR_ATTACHMENT1_EXT,GL_COLOR_ATTACHMENT3_EXT,GL_COLOR_ATTACHMENT7_EXT};
      glDrawBuffers(3,buffers);
    }
    else{ 
      GLenum buffers[] = {GL_COLOR_ATTACHMENT0_EXT,GL_COLOR_ATTACHMENT2_EXT,GL_COLOR_ATTACHMENT7_EXT};
      glDrawBuffers(3,buffers);
    }

  }
  else{

    if(odd){
      GLenum buffers[] = {GL_COLOR_ATTACHMENT1_EXT,GL_COLOR_ATTACHMENT3_EXT};
      glDrawBuffers(2,buffers);
    }
    else{ 
      GLenum buffers[] = {GL_COLOR_ATTACHMENT0_EXT,GL_COLOR_ATTACHMENT2_EXT};
      glDrawBuffers(2,buffers);
    }

  }

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glViewport(0, 0, twidth, theight);
  
  
  glEnable(texType);
  multipleBuildings_shader.activate();

  if(advect_terms != 0){
    glUniform1iARB(uniform_colorAdvectTerms, 1);
  }
  else{
    glUniform1iARB(uniform_colorAdvectTerms, 0);
  }


  glUniform1fARB(uniform_timeStep, time_step);
  // generate and set the random texture coordinate offset
  // 
  if (texType == GL_TEXTURE_RECTANGLE_ARB)
    {
      // texture coordinates will range from 0 to W in width and 0 to
      // H in height, so generate random value in this range
      //float f1 = Random::uniform() * twidth;
      //float f2 = Random::uniform() * theight;
      //if(!osgPlume)
      //setRandomTexCoords();
      float f1 = Random::uniform() * twidth;
      float f2 = Random::uniform() * theight;
      //output << t1 << "\n";
      //output << t2 << "\n";

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

  //Bind the cell type to TEXTURE UNIT 8
  glActiveTexture(GL_TEXTURE8);
  glBindTexture(texType, cellType);
  glUniform1iARB(uniform_cellType, 8);

  //Bind the building to TEXTURE UNIT 7
  glActiveTexture(GL_TEXTURE7);
  glBindTexture(texType, buildings);
  glUniform1iARB(uniform_buildings, 7);

  //Bind the du/dz texture to TEXTURE UNIT 6
  glActiveTexture(GL_TEXTURE6);
  glBindTexture(texType, duvw_dz);
  glUniform1iARB(uniform_duvw_dz, 6);

  //Bind the tau texture to TEXTURE UNIT 5
  glActiveTexture(GL_TEXTURE5);
  glBindTexture(texType, tau_dz);
  glUniform1iARB(uniform_tau_dz, 5);

  //Bind the lambda texture to TEXTURE UNIT 4
  glActiveTexture(GL_TEXTURE4);
  glBindTexture(texType, lambda);
  glUniform1iARB(uniform_lambda, 4);

  // Bind the random data field to TEXTURE UNIT 3
  glActiveTexture(GL_TEXTURE3);
  glBindTexture(texType, randomValues);
  glUniform1iARB(uniform_random, 3);

  
  //Bind previous prime values to TEXTURE UNIT 2
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
  if (odd)
    glBindTexture(texType, positions0);  // read from texture 0
  else 
    glBindTexture(texType, positions1);  // read from texture 1
 
  glUniform1iARB(uniform_postex, 0);

  glBegin(GL_QUADS);
  {
    glTexCoord2f(0, 0);					glVertex3f(-1, -1, -0.5f);
    glTexCoord2f(float(twidth), 0);			glVertex3f( 1, -1, -0.5f);
    glTexCoord2f(float(twidth), float(theight));	glVertex3f( 1,  1, -0.5f);
    glTexCoord2f(0, float(theight));			glVertex3f(-1,  1, -0.5f);
  }
  glEnd();
  
  multipleBuildings_shader.deactivate();
 
  glBindTexture(texType, 0);
  
  //Prints out the updated prime values
  if(outputPrime)
    printPrime(odd,false);

}

void ParticleControl::setupReflectionShader(float life_time){
  reflection_shader.addShader("Shaders/reflectionAdvect_vp.glsl", GLSLObject::VERTEX_SHADER);
  reflection_shader.addShader("Shaders/reflectionAdvect_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  reflection_shader.createProgram();

  // Get location of the sampler uniform
  uniform_postex = reflection_shader.createUniform("pos_texunit");
  uniform_wind = reflection_shader.createUniform("wind_texunit");
  uniform_timeStep = reflection_shader.createUniform("time_step");
  uniform_primePrev = reflection_shader.createUniform("primePrev");
  uniform_random = reflection_shader.createUniform("random");
  uniform_lambda = reflection_shader.createUniform("lambda");
  uniform_tau_dz = reflection_shader.createUniform("tau_dz");
  uniform_duvw_dz = reflection_shader.createUniform("duvw_dz");
  uniform_randomTexCoordOffset = reflection_shader.createUniform("random_texCoordOffset");
  uniform_randomTexWidth = reflection_shader.createUniform("random_texWidth");
  uniform_randomTexHeight = reflection_shader.createUniform("random_texHeight");
  
  //Building uniform variables
  uniform_xfo = reflection_shader.createUniform("xfo");
  uniform_yfo = reflection_shader.createUniform("yfo");
  uniform_zfo = reflection_shader.createUniform("zfo");
  uniform_ht = reflection_shader.createUniform("ht");
  uniform_wti = reflection_shader.createUniform("wti");
  uniform_lti = reflection_shader.createUniform("lti");

  GLint ulifeTime = reflection_shader.createUniform("life_time");
  GLint unx = reflection_shader.createUniform("nx");
  GLint uny = reflection_shader.createUniform("ny");
  GLint unz = reflection_shader.createUniform("nz");
  GLint unxdx = reflection_shader.createUniform("nxdx");
  GLint unydy = reflection_shader.createUniform("nydy");
  GLint unzdz = reflection_shader.createUniform("nzdz");
  GLint uNumInRow = reflection_shader.createUniform("numInRow");

  reflection_shader.activate();

  glUniform1fARB(ulifeTime, life_time);
  glUniform1iARB(unx, nx);
  glUniform1iARB(uny, ny);
  glUniform1iARB(unz, nz);
  glUniform1iARB(unxdx, nxdx);
  glUniform1iARB(unydy, nydy);
  glUniform1iARB(unzdz, nzdz);
  glUniform1iARB(uNumInRow, numInRow);

  reflection_shader.deactivate();
}
void ParticleControl::reflectionAdvect(bool odd, GLuint windField, GLuint positions0, 
			     GLuint positions1, GLuint prime0, GLuint prime1, 
			     GLuint randomValues,GLuint lambda, GLuint tau_dz, 
			     GLuint duvw_dz, float time_step, float* buildParam)
{

  /*std::ofstream output;
  if(!alreadyOpen){
    output.open(randomFile.c_str());
    alreadyOpen = true;
  }
  else{
    output.open(randomFile.c_str(),std::ios::app);
    }*/
  

  //Prints out the previous prime values
  if(outputPrime)
    printPrime(odd, true);

  if(odd){
    GLenum buffers[] = {GL_COLOR_ATTACHMENT1_EXT,GL_COLOR_ATTACHMENT3_EXT};
    glDrawBuffers(2,buffers);
  }
  else{ 
    GLenum buffers[] = {GL_COLOR_ATTACHMENT0_EXT,GL_COLOR_ATTACHMENT2_EXT};
    glDrawBuffers(2,buffers);
  }

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glViewport(0, 0, twidth, theight);
  
  
  glEnable(texType);
  reflection_shader.activate();

  //Building varibles
  if(buildParam[0] > 0){
    glUniform1fARB(uniform_xfo, buildParam[1]);
    glUniform1fARB(uniform_yfo, buildParam[2]);
    glUniform1fARB(uniform_zfo, buildParam[3]);
    glUniform1fARB(uniform_ht, buildParam[4]);
    glUniform1fARB(uniform_wti, buildParam[5]);
    glUniform1fARB(uniform_lti, buildParam[6]);
  }

  glUniform1fARB(uniform_timeStep, time_step);
  // generate and set the random texture coordinate offset
  // 
  if (texType == GL_TEXTURE_RECTANGLE_ARB)
    {
      // texture coordinates will range from 0 to W in width and 0 to
      // H in height, so generate random value in this range
      //float f1 = Random::uniform() * twidth;
      //float f2 = Random::uniform() * theight;
      //if(!osgPlume)
      //setRandomTexCoords();
      float f1 = Random::uniform() * twidth;
      float f2 = Random::uniform() * theight;
      //output << t1 << "\n";
      //output << t2 << "\n";

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

  //Bind the lambda texture to TEXTURE UNIT 6
  glActiveTexture(GL_TEXTURE6);
  glBindTexture(texType, duvw_dz);
  glUniform1iARB(uniform_duvw_dz, 6);

  //Bind the lambda texture to TEXTURE UNIT 5
  glActiveTexture(GL_TEXTURE5);
  glBindTexture(texType, tau_dz);
  glUniform1iARB(uniform_tau_dz, 5);

  //Bind the lambda texture to TEXTURE UNIT 4
  glActiveTexture(GL_TEXTURE4);
  glBindTexture(texType, lambda);
  glUniform1iARB(uniform_lambda, 4);

  // Bind the random data field to TEXTURE UNIT 3
  glActiveTexture(GL_TEXTURE3);
  glBindTexture(texType, randomValues);
  glUniform1iARB(uniform_random, 3);

  
  //Bind previous prime values to TEXTURE UNIT 2
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
  if (odd)
    glBindTexture(texType, positions0);  // read from texture 0
  else 
    glBindTexture(texType, positions1);  // read from texture 1
 
  glUniform1iARB(uniform_postex, 0);

  glBegin(GL_QUADS);
  {
    glTexCoord2f(0, 0);								glVertex3f(-1, -1, -0.5f);
    glTexCoord2f(float(twidth), 0);					glVertex3f( 1, -1, -0.5f);
    glTexCoord2f(float(twidth), float(theight));	glVertex3f( 1,  1, -0.5f);
    glTexCoord2f(0, float(theight));				glVertex3f(-1,  1, -0.5f);
  }
  glEnd();
  
  reflection_shader.deactivate();
 
  glBindTexture(texType, 0);
  
  //Prints out the updated prime values
  if(outputPrime)
    printPrime(odd,false);

}
void ParticleControl::setupNonGaussianShader(float life_time){
  nonGaussian_shader.addShader("Shaders/nonGaussianAdvect_vp.glsl", GLSLObject::VERTEX_SHADER);
  nonGaussian_shader.addShader("Shaders/nonGaussianAdvect_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  nonGaussian_shader.createProgram();

  // Get location of the sampler uniform
  uniform_postex = nonGaussian_shader.createUniform("pos_texunit");
  uniform_wind = nonGaussian_shader.createUniform("wind_texunit");
  uniform_timeStep = nonGaussian_shader.createUniform("time_step");
  uniform_primePrev = nonGaussian_shader.createUniform("primePrev");
  uniform_random = nonGaussian_shader.createUniform("random");
  uniform_lambda = nonGaussian_shader.createUniform("lambda");
  uniform_tau_dz = nonGaussian_shader.createUniform("tau_dz");
  uniform_duvw_dz = nonGaussian_shader.createUniform("duvw_dz");
  uniform_randomTexCoordOffset = nonGaussian_shader.createUniform("random_texCoordOffset");
  uniform_randomTexWidth = nonGaussian_shader.createUniform("random_texWidth");
  uniform_randomTexHeight = nonGaussian_shader.createUniform("random_texHeight");

  GLint ulifeTime = nonGaussian_shader.createUniform("life_time");
  GLint unx = nonGaussian_shader.createUniform("nx");
  GLint uny = nonGaussian_shader.createUniform("ny");
  GLint unz = nonGaussian_shader.createUniform("nz");
  GLint unxdx = nonGaussian_shader.createUniform("nxdx");
  GLint unydy = nonGaussian_shader.createUniform("nydy");
  GLint unzdz = nonGaussian_shader.createUniform("nzdz");
  GLint uNumInRow = nonGaussian_shader.createUniform("numInRow");

  nonGaussian_shader.activate();

  glUniform1fARB(ulifeTime, life_time);
  glUniform1iARB(unx, nx);
  glUniform1iARB(uny, ny);
  glUniform1iARB(unz, nz);
  glUniform1iARB(unxdx, nxdx);
  glUniform1iARB(unydy, nydy);
  glUniform1iARB(unzdz, nzdz);
  glUniform1iARB(uNumInRow, numInRow);

  nonGaussian_shader.deactivate();
}
void ParticleControl::nonGaussianAdvect(bool odd, GLuint windField, GLuint positions0, 
			     GLuint positions1, GLuint prime0, GLuint prime1, 
			     GLuint randomValues,GLuint lambda, GLuint tau_dz, 
			     GLuint duvw_dz, float time_step)
{

  //Prints out the previous prime values
  if(outputPrime)
    printPrime(odd, true);

  if (odd){
    GLenum buffers[] = {GL_COLOR_ATTACHMENT1_EXT,GL_COLOR_ATTACHMENT3_EXT};
    glDrawBuffers(2,buffers);
  }
  else{ 
    GLenum buffers[] = {GL_COLOR_ATTACHMENT0_EXT,GL_COLOR_ATTACHMENT2_EXT};
    glDrawBuffers(2,buffers);
  }

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glViewport(0, 0, twidth, theight);
  
  
  glEnable(texType);
  nonGaussian_shader.activate();

  glUniform1fARB(uniform_timeStep, time_step);
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

  //Bind the lambda texture to TEXTURE UNIT 6
  glActiveTexture(GL_TEXTURE6);
  glBindTexture(texType, duvw_dz);
  glUniform1iARB(uniform_duvw_dz, 6);

  //Bind the lambda texture to TEXTURE UNIT 5
  glActiveTexture(GL_TEXTURE5);
  glBindTexture(texType, tau_dz);
  glUniform1iARB(uniform_tau_dz, 5);

  //Bind the lambda texture to TEXTURE UNIT 4
  glActiveTexture(GL_TEXTURE4);
  glBindTexture(texType, lambda);
  glUniform1iARB(uniform_lambda, 4);

  // Bind the random data field to TEXTURE UNIT 3
  glActiveTexture(GL_TEXTURE3);
  glBindTexture(texType, randomValues);
  glUniform1iARB(uniform_random, 3);

  
  //Bind previous prime values to TEXTURE UNIT 2
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
    glTexCoord2f(0, 0);								glVertex3f(-1, -1, -0.5f);
    glTexCoord2f(float(twidth), 0);					glVertex3f( 1, -1, -0.5f);
    glTexCoord2f(float(twidth), float(theight));	glVertex3f( 1,  1, -0.5f);
    glTexCoord2f(0, float(theight));				glVertex3f(-1,  1, -0.5f);
  }
  glEnd();
  
  nonGaussian_shader.deactivate();
 
  glBindTexture(texType, 0);
  
  //Prints out the updated prime values
  if(outputPrime)
    printPrime(odd,false);

}

void ParticleControl::setupPrime_and_AdvectShader(float life_time){
  //This shader is used to move the particles
  
  mrt_shader.addShader("Shaders/updatePrime_and_Advect_vp.glsl", GLSLObject::VERTEX_SHADER);
  mrt_shader.addShader("Shaders/updatePrime_and_Advect_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  mrt_shader.createProgram();
   
  // Get location of the sampler uniform
  uniform_postex = mrt_shader.createUniform("pos_texunit");
  uniform_wind = mrt_shader.createUniform("wind_texunit");
  uniform_timeStep = mrt_shader.createUniform("time_step");
  uniform_primePrev = mrt_shader.createUniform("primePrev");
  uniform_random = mrt_shader.createUniform("random");
  uniform_lambda = mrt_shader.createUniform("lambda");
  uniform_randomTexCoordOffset = mrt_shader.createUniform("random_texCoordOffset");
  uniform_randomTexWidth = mrt_shader.createUniform("random_texWidth");
  uniform_randomTexHeight = mrt_shader.createUniform("random_texHeight");


  GLint ulifeTime = mrt_shader.createUniform("life_time");
  GLint unx = mrt_shader.createUniform("nx");
  GLint uny = mrt_shader.createUniform("ny");
  GLint unz = mrt_shader.createUniform("nz");
  GLint unxdx = mrt_shader.createUniform("nxdx");
  GLint unydy = mrt_shader.createUniform("nydy");
  GLint unzdz = mrt_shader.createUniform("nzdz");
  GLint uNumInRow = mrt_shader.createUniform("numInRow");
  mrt_shader.activate();

  glUniform1fARB(ulifeTime, life_time);
  glUniform1iARB(unx, nx);
  glUniform1iARB(uny, ny);
  glUniform1iARB(unz, nz);
  glUniform1iARB(unxdx, nxdx);
  glUniform1iARB(unydy, nydy);
  glUniform1iARB(unzdz, nzdz);
  glUniform1iARB(uNumInRow, numInRow);

  mrt_shader.deactivate();


}
void ParticleControl::updatePrimeAndAdvect(bool odd, 
			     GLuint windField, GLuint positions0, GLuint positions1, 
			     GLuint prime0, GLuint prime1, GLuint randomValues, 
				  GLuint lambda, float time_step)
{

  //Prints out the previous prime values
  if(outputPrime)
    printPrime(odd, true);

  if (odd){
    GLenum buffers[] = {GL_COLOR_ATTACHMENT1_EXT,GL_COLOR_ATTACHMENT3_EXT};
    glDrawBuffers(2,buffers);
  }
  else{ 
    GLenum buffers[] = {GL_COLOR_ATTACHMENT0_EXT,GL_COLOR_ATTACHMENT2_EXT};
    glDrawBuffers(2,buffers);
  }

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glViewport(0, 0, twidth, theight);
  
  
  glEnable(texType);
  mrt_shader.activate();

  glUniform1fARB(uniform_timeStep, time_step);
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


  //Bind the lambda texture to TEXTURE UNIT 4
  glActiveTexture(GL_TEXTURE4);
  glBindTexture(texType, lambda);
  glUniform1iARB(uniform_lambda, 4);

  // Bind the random data field to TEXTURE UNIT 3
  glActiveTexture(GL_TEXTURE3);
  glBindTexture(texType, randomValues);
  glUniform1iARB(uniform_random, 3);

  
  //Bind previous prime values to TEXTURE UNIT 2
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
    glTexCoord2f(0, 0);					glVertex3f(-1, -1, -0.5f);
    glTexCoord2f(float(twidth), 0);			glVertex3f( 1, -1, -0.5f);
    glTexCoord2f(float(twidth), float(theight));	glVertex3f( 1,  1, -0.5f);
    glTexCoord2f(0, float(theight));			glVertex3f(-1,  1, -0.5f);
  }
  glEnd();
  
  mrt_shader.deactivate();
 
  glBindTexture(texType, 0);
  
  //Prints out the updated prime values
  if(outputPrime)
    printPrime(odd,false);

}

void ParticleControl::setupPrimeShader(){  //Included argument -- Balli(04/12/07)
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
  GLint unxdx = prime_shader.createUniform("nxdx");
  GLint unydy = prime_shader.createUniform("nydy");
  GLint unzdz = prime_shader.createUniform("nzdz");
  GLint uNumInRow = prime_shader.createUniform("numInRow");

  prime_shader.activate();

  glUniform1iARB(unx, nx);
  glUniform1iARB(uny, ny);
  glUniform1iARB(unz, nz);
  glUniform1iARB(unxdx, nxdx);
  glUniform1iARB(unydy, nydy);
  glUniform1iARB(unzdz, nzdz);
  glUniform1iARB(uNumInRow, numInRow);

  prime_shader.deactivate();
  //End --Balli(04/12/07)
}
void ParticleControl::updatePrime(bool odd, GLuint positions0,GLuint positions1,GLuint prime0, 
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
    glTexCoord2f(0, 0);								glVertex3f(-1, -1, -0.5f);
    glTexCoord2f(float(twidth), 0);					glVertex3f( 1, -1, -0.5f);
    glTexCoord2f(float(twidth), float(theight));	glVertex3f( 1,  1, -0.5f);
    glTexCoord2f(0, float(theight));				glVertex3f(-1,  1, -0.5f);
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
void ParticleControl::setupAdvectShader(float life_time){

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
  GLint unxdx = pass1_shader.createUniform("nxdx");
  GLint unydy = pass1_shader.createUniform("nydy");
  GLint unzdz = pass1_shader.createUniform("nzdz");
  GLint uNumInRow = pass1_shader.createUniform("numInRow");

  pass1_shader.activate();

  glUniform1fARB(ulifeTime, life_time);
  glUniform1iARB(unx, nx);
  glUniform1iARB(uny, ny);
  glUniform1iARB(unz, nz);
  glUniform1iARB(unxdx, nxdx);
  glUniform1iARB(unydy, nydy);
  glUniform1iARB(unzdz, nzdz);
  glUniform1iARB(uNumInRow, numInRow);

  pass1_shader.deactivate();

}
void ParticleControl::advect(bool odd, 
			     GLuint windField, GLuint positions0, GLuint positions1, 
			     GLuint prime0, GLuint prime1, float time_step)
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
    glTexCoord2f(0, 0);								glVertex3f(-1, -1, -0.5f);
    glTexCoord2f(float(twidth), 0);					glVertex3f( 1, -1, -0.5f);
    glTexCoord2f(float(twidth), float(theight));	glVertex3f( 1,  1, -0.5f);
    glTexCoord2f(0, float(theight));				glVertex3f(-1,  1, -0.5f);
  }
  glEnd();
  
  pass1_shader.deactivate();
 
  glBindTexture(texType, 0);

}
void ParticleControl::setupParticleColor_shader(){
  currDir_shader.addShader("Shaders/direction_to_color_vp.glsl",GLSLObject::VERTEX_SHADER);
  currDir_shader.addShader("Shaders/direction_to_color_fp.glsl",GLSLObject::FRAGMENT_SHADER);
  currDir_shader.createProgram();

  uniform_currentPrime = currDir_shader.createUniform("currPrime");
  uniform_windVelocity = currDir_shader.createUniform("windVel");
  uniform_partPos = currDir_shader.createUniform("position");
  uniform_prevPartPos = currDir_shader.createUniform("position_prev");
  
  GLint unir = currDir_shader.createUniform("numInRow");

  GLint unx = currDir_shader.createUniform("nx");
  GLint uny = currDir_shader.createUniform("ny");
  GLint unz = currDir_shader.createUniform("nz");

  currDir_shader.activate();

  glUniform1iARB(unx, nxdx);
  glUniform1iARB(uny, nydy);
  glUniform1iARB(unz, nzdz);
  glUniform1iARB(unir,numInRow);

  currDir_shader.deactivate();

}
void ParticleControl::updateParticleColors(bool odd, GLuint prime0,GLuint prime1,GLuint windField,
				    GLuint positions0,GLuint positions1){

  glDrawBuffer(particleColorBuffer);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glViewport(0, 0, twidth, theight);
    
  glEnable(texType);
  currDir_shader.activate();

  glActiveTexture(GL_TEXTURE3);
  if(odd)
    glBindTexture(texType,positions0);
  else
    glBindTexture(texType,positions1);
  glUniform1iARB(uniform_prevPartPos,3);

  glActiveTexture(GL_TEXTURE2);
  if(odd)
    glBindTexture(texType,positions1);
  else
    glBindTexture(texType,positions0);
  glUniform1iARB(uniform_partPos,2);

  glActiveTexture(GL_TEXTURE1);
  glUniform1iARB(uniform_windVelocity, 1);
  glBindTexture(texType, windField);

  glActiveTexture(GL_TEXTURE0);
  glUniform1iARB(uniform_currentPrime, 0);
  if(odd)
    glBindTexture(texType, prime1);
  else
    glBindTexture(texType, prime0);


  glBegin(GL_QUADS);
  {
    glTexCoord2f(0, 0);								glVertex3f(-1, -1, -0.5f);
    glTexCoord2f(float(twidth), 0);					glVertex3f( 1, -1, -0.5f);
    glTexCoord2f(float(twidth), float(theight));	glVertex3f( 1,  1, -0.5f);
    glTexCoord2f(0, float(theight));				glVertex3f(-1,  1, -0.5f);
  }
  glEnd();

  currDir_shader.deactivate();

  glBindTexture(texType,0);
}
void ParticleControl::setupMeanVel_shader(){
  meanVel_shader.addShader("Shaders/meanVel_vp.glsl", GLSLObject::VERTEX_SHADER);
  meanVel_shader.addShader("Shaders/meanVel_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  meanVel_shader.createProgram();

  // Get location of the sampler uniform
  uniform_prevMean = meanVel_shader.createUniform("prevMean");
  uniform_currVel = meanVel_shader.createUniform("currVel");
  uniform_position = meanVel_shader.createUniform("position");
  uniform_windVel = meanVel_shader.createUniform("windVel");

  GLint unir = meanVel_shader.createUniform("numInRow");

  GLint unx = meanVel_shader.createUniform("nx");
  GLint uny = meanVel_shader.createUniform("ny");
  GLint unz = meanVel_shader.createUniform("nz");

  meanVel_shader.activate();

  glUniform1iARB(unx, nxdx);
  glUniform1iARB(uny, nydy);
  glUniform1iARB(unz, nzdz);
  glUniform1iARB(unir,numInRow);

  meanVel_shader.deactivate();


}

void ParticleControl::findMeanVel(bool odd,GLuint prime0,GLuint prime1,
				  GLuint meanVel0,GLuint meanVel1,
				  GLuint positions0,GLuint positions1,
				  GLuint windField){

  if(odd){
    glDrawBuffer(meanVelBuffer1);
  }
  else 
    glDrawBuffer(meanVelBuffer0);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glViewport(0, 0, twidth, theight);
    
  glEnable(texType);
  meanVel_shader.activate();

  glActiveTexture(GL_TEXTURE3);
  glBindTexture(texType, windField);
  glUniform1iARB(uniform_windVel,3);

  glActiveTexture(GL_TEXTURE2);
  if(odd)
    glBindTexture(texType,positions1);
  else
    glBindTexture(texType,positions0);
  glUniform1iARB(uniform_position,2);

  glActiveTexture(GL_TEXTURE1);
  if(odd)
    glBindTexture(texType,prime1);
  else
    glBindTexture(texType,prime0);
  glUniform1iARB(uniform_currVel,1);

  glActiveTexture(GL_TEXTURE0);
  glUniform1iARB(uniform_prevMean, 0);
  if(odd)
    glBindTexture(texType, meanVel0);
  else
    glBindTexture(texType, meanVel1);


  glBegin(GL_QUADS);
  {
    glTexCoord2f(0, 0);								glVertex3f(-1, -1, -0.5f);
    glTexCoord2f(float(twidth), 0);					glVertex3f( 1, -1, -0.5f);
    glTexCoord2f(float(twidth), float(theight));	glVertex3f( 1,  1, -0.5f);
    glTexCoord2f(0, float(theight));				glVertex3f(-1,  1, -0.5f);
  }
  glEnd();

  meanVel_shader.deactivate();

  glBindTexture(texType,0);

}

void ParticleControl::printPositions(bool odd){
  
  glGetIntegerv(GL_READ_BUFFER, &currentbuffer);

  if(odd)
    glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
  else
    glReadBuffer(GL_COLOR_ATTACHMENT1_EXT);

  std::cout << "Previous Positions" << std::endl;
  dumpContents();

  if(odd)
    glReadBuffer(GL_COLOR_ATTACHMENT1_EXT);
  else
    glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
       
  std::cout << "Updated Positions" << std::endl;
  dumpContents();

  glReadBuffer(currentbuffer);

}
void ParticleControl::printMeanVelocities(bool odd){
  glGetIntegerv(GL_READ_BUFFER, &currentbuffer);
  double tts;
  //std::ofstream output;

  //output.open("MeanVel.dat");

  std::cout << "Mean Velocities" << std::endl;

  if(odd)
    glReadBuffer(meanVelBuffer1);
  else
    glReadBuffer(meanVelBuffer0);

  buffer_mem = new GLfloat[ twidth * theight * 4 ];  
  glReadPixels(0, 0, twidth, theight, GL_RGBA, GL_FLOAT, buffer_mem);

  std::cout << "IDX  X     Y     Z" << std::endl;
  //output << "IDX  X     Y     Z  \n" ;

  int pn =0;
  for (int j=0; j<theight; j++)
    for (int i=0; i<twidth; i++){
	  
       int idx = j*twidth*4 + i*4;
       
       tts = buffer_mem[idx+3];

       std::cout << pn << " ";
       std::cout << buffer_mem[idx]/tts << " ";
       std::cout << buffer_mem[idx+1]/tts << " ";
       std::cout << buffer_mem[idx+2]/tts << " ";
       std::cout << "#ts: " << buffer_mem[idx+3] << std::endl;
       /*output << pn << " ";
       output << buffer_mem[idx]/tts << " ";
       output << buffer_mem[idx+1]/tts << " ";
       output << buffer_mem[idx+2]/tts << "\n";
       //output << "#ts: " << buffer_mem[idx+3] << "\n";*/


       pn++;
    }
  delete [] buffer_mem;

  //output.close();

  glReadBuffer(currentbuffer);
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

void ParticleControl::createPrimeImages(bool odd){
  glGetIntegerv(GL_READ_BUFFER, &currentbuffer);

  if(odd)
    glReadBuffer(GL_COLOR_ATTACHMENT2_EXT);
  else
    glReadBuffer(GL_COLOR_ATTACHMENT3_EXT);

  buffer_mem = new GLfloat[ twidth * theight * 4 ];
  glReadPixels(0, 0, twidth, theight, GL_RGBA, GL_FLOAT, buffer_mem);
  
  int alpha = 3;

  //Find min and max of previous position values
  max = buffer_mem[0];
  min = buffer_mem[0];

  for(int j=1; j < theight*twidth*4; j++){
    if(j != alpha){
      if(buffer_mem[j] > max)
	max = buffer_mem[j];
      if(buffer_mem[j] < min)
	min = buffer_mem[j];
    }
    else
      alpha +=4;
  }
  //Now update min and max with updated position values
  if(odd)
    glReadBuffer(GL_COLOR_ATTACHMENT3_EXT);
  else
    glReadBuffer(GL_COLOR_ATTACHMENT2_EXT);

  buffer_mem = new GLfloat[ twidth * theight * 4 ];
  glReadPixels(0, 0, twidth, theight, GL_RGBA, GL_FLOAT, buffer_mem);
  
  alpha = 3;
  for(int j=0; j < theight*twidth*4; j++){
    if(j != alpha){
      if(buffer_mem[j] > max)
	max = buffer_mem[j];
      if(buffer_mem[j] < min)
	min = buffer_mem[j];
    }
    else
      alpha +=4;
  }

  //Now create the ppm image files
  if(odd)
    glReadBuffer(GL_COLOR_ATTACHMENT2_EXT);
  else
    glReadBuffer(GL_COLOR_ATTACHMENT3_EXT);

  writePPM("prevPrime.ppm");

  if(odd)
    glReadBuffer(GL_COLOR_ATTACHMENT3_EXT);
  else
    glReadBuffer(GL_COLOR_ATTACHMENT2_EXT);

  writePPM("updatePrime.ppm");

  glReadBuffer(currentbuffer);

}
void ParticleControl::createPosImages(bool odd){

  glGetIntegerv(GL_READ_BUFFER, &currentbuffer);

  if(odd)
    glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
  else
    glReadBuffer(GL_COLOR_ATTACHMENT1_EXT);

  buffer_mem = new GLfloat[ twidth * theight * 4 ];
  glReadPixels(0, 0, twidth, theight, GL_RGBA, GL_FLOAT, buffer_mem);
  
  int alpha = 3;

  //Find min and max of previous position values
  max = buffer_mem[0];
  min = buffer_mem[0];

  for(int j=1; j < theight*twidth*4; j++){
    if(j != alpha){
      if(buffer_mem[j] > max)
	max = buffer_mem[j];
      if(buffer_mem[j] < min)
	min = buffer_mem[j];
    }
    else
      alpha +=4;
  }
  //Now update min and max with updated position values
  if(odd)
    glReadBuffer(GL_COLOR_ATTACHMENT1_EXT);
  else
    glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);

  buffer_mem = new GLfloat[ twidth * theight * 4 ];
  glReadPixels(0, 0, twidth, theight, GL_RGBA, GL_FLOAT, buffer_mem);
  
  alpha = 3;
  for(int j=0; j < theight*twidth*4; j++){
    if(j != alpha){
      if(buffer_mem[j] > max)
	max = buffer_mem[j];
      if(buffer_mem[j] < min)
	min = buffer_mem[j];
    }
    else
      alpha +=4;
  }

  //Now create the ppm image files
  if(odd)
    glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
  else
    glReadBuffer(GL_COLOR_ATTACHMENT1_EXT);

  writePPM("prevPos.ppm");

  if(odd)
    glReadBuffer(GL_COLOR_ATTACHMENT1_EXT);
  else
    glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);

  writePPM("updatePos.ppm");

  glReadBuffer(currentbuffer);
}
void ParticleControl::writePPM(const std::string& filename)
{

  buffer_mem = new GLfloat[ twidth * theight * 4 ];
  glReadPixels(0, 0, twidth, theight, GL_RGBA, GL_FLOAT, buffer_mem);

  std::ofstream outfile;
  outfile.open(filename.c_str(), std::ios::out | std::ios::trunc);
  if (!(outfile.is_open()))
    {
      std::cerr << "ERROR: Could not open " << filename << " for writing output image." << std::endl;
      exit(-1);
    }

  outfile << "P3\n" << twidth << " " << theight << "\n255\n";
 
  int idx;
  int count = 1;
  int N = twidth;
  
  //std::cout << filename << " max is " << max << " min is " << min << std::endl;
  std::ostringstream strstream;
  //for (int j=theight-1; j>=0; j--)
  for(int j=0; j<theight; j++)
    for (int i=0; i<twidth; i++)
      {
    idx = j*twidth*4 + i*4;
     
    // output a RGB tuple to the string
      strstream << c2Short(buffer_mem[idx]) << ' '
          << c2Short(buffer_mem[idx+1]) << ' '
          << c2Short(buffer_mem[idx+2]) <<  ' ';

    /*std::cout << c2Short(buffer_mem[idx]) << ' '
          << c2Short(buffer_mem[idx+1]) << ' '
	      << c2Short(buffer_mem[idx+2]) <<  std::endl;

    std::cout << buffer_mem[idx] << ' '
          << buffer_mem[idx+1] << ' '
	  << buffer_mem[idx+2] <<  std::endl;*/
     
    // put at least N tuples on a line
    if (count%N == 0)
      {
        outfile << strstream.str();
        outfile << std::endl;
        strstream.str("");
      }
    count++;
      }
  outfile << strstream.str();
  outfile << std::endl;
 
  outfile.close();
}
short ParticleControl::c2Short(float num){
  
  float range = max - min;
  num = num - min;
  num = num/range;
  num = num*255;
  
  return (short)num;

}
void ParticleControl::createTexture(GLuint texId, GLenum format, int w, int h, GLfloat* data){
  glBindTexture(texType, texId);
  glTexParameteri(texType, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(texType, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(texType, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(texType, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  glTexImage2D(texType, 0, format, w, h, 0, GL_RGBA, GL_FLOAT, data);

}
void ParticleControl::createIntTexture(GLuint texId, GLenum format, int w, int h, GLint* data){
  glBindTexture(texType, texId);
  glTexParameteri(texType, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(texType, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(texType, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(texType, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  glTexImage2D(texType, 0, format, w, h, 0, GL_RGBA, GL_INT, data);

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

void ParticleControl::initWindTex(GLuint windField, int* numberInRow, int dataSet){
  // Create wind velocity data texture
  int arrSize = nxdx*nydy*nzdz;

  wind_vel = new wind[arrSize];
  cellQuic = new cellType[arrSize];

  //Matrix tau11,tau22,tau33,and tau13
  tau = new Matrix[arrSize];

  //sigU,sigV,and sigW
  sig = new wind[arrSize];

  switch(dataSet){

  case 1:
    test1();
    break;
  case 2:
    randomWindField();
    break;
  case 3:
    uniformUWindField();
    break;
  case 4:
    variedUWindField();
    break;
  case 5:
    QUICWindField();
    break;
  case 6:
    QUICWindField();
    break;

  }

  /////////////////////////////////////////////////////////
  //Calculate width and height for wind texture
  //
  //This tries to minimize the width and height values
  //to try and fit the wind field into a 2D texture without
  //wasting too much space.  
  /////////////////////////////////////////////////////////
  int total = arrSize;
  width = (int)sqrt((float)total);
  
  int scaler;
  if(nydy > nxdx) scaler = nydy;
  else scaler = nxdx;

  width = width - (width%scaler);
  
  bool done = false;
  while(!done){ 
    int num = width/scaler;
    if((num*num) >= nzdz){
      done = true;
    }
    else{
      width = width+scaler;
    }
  }
  if(width%2 != 0) width++;
  height = width;
 
  ////////////////////////////////////////////////////////
  //Convert this to data array for a texture
  //
  //This will directly put the 3D data into an array
  //that is used to make the 2D texture.
  ///////////////////////////////////////////////////////
  
  (*numberInRow) = (width - (width % nxdx))/nxdx;
  numInRow = *numberInRow;

  if(dataSet < 5){

    int qi, qj, qk;
    int p2idx = 0, texidx = 0;
    int row = 0;
  
    GLfloat *data = new GLfloat[ width * height * 4 ];
   
    for (qk=0; qk<nzdz; qk++) 
      for (qi=0; qi<nydy; qi++)
	for (qj=0; qj<nxdx; qj++)
	  {
	    p2idx = qk*nydy*nxdx + qi*nxdx + qj;
	    
	    row = qk / (numInRow);
	    texidx = row * width * nydy * 4 +
	      qi * width * 4 +
	      qk % (numInRow) * nxdx * 4 +
	      qj * 4;
	  
	    data[texidx] = wind_vel[p2idx].u;
	    data[texidx+1] = wind_vel[p2idx].v;
	    data[texidx+2] = wind_vel[p2idx].w;	 
	    data[texidx+3] = float((0.5*5.7)*(ustar*ustar*ustar)/(0.4*(qk+1)));//This value is the '0.5*CoEps' value	
	  }

    createTexture(windField, GL_RGBA32F_ARB, width, height, data);

    delete [] data;

  }
}
void ParticleControl::initLambda_and_TauTex(GLuint lambda, GLuint tau_dz, GLuint duvw_dz){
  GLfloat *data = new GLfloat[ width * height * 4 ];
  GLfloat *dataTwo = new GLfloat[ width * height * 4 ];
  GLfloat *data3 = new GLfloat[width*height*4];

  int qi, qj, qk;
  int p2idx = 0, texidx = 0, texidxBelow = 0;
  int row = 0, rowBelow = 0;
  
  //cell below, cell above, cell two below
  int idxBelow, idxAbove, idx2Below;

  float du_dz;

  //Cell height
  //This is now a variable set in settings file
  //float dz = 1.0;

  //Surface roughness
  float znaut = 0.01f;
  
  for (qk=0; qk<nzdz; qk++) 
    for (qi=0; qi<nydy; qi++)
      for (qj=0; qj<nxdx; qj++)
	{
	  p2idx = qk*nydy*nxdx + qi*nxdx + qj;
	    
	  row = qk / (numInRow);
	  texidx = row * width * nydy * 4 +
	  qi * width * 4 +
	  qk % (numInRow) * nxdx * 4 +
	  qj * 4;

	  rowBelow = (qk-1) / (numInRow);
	  texidxBelow = rowBelow * width * nydy * 4 +
	    qi * width * 4 +
	    (qk-1) % (numInRow) * nxdx * 4 +
	    qj * 4;
	  
	  //The cells right above and below the current cell
	  idxBelow = (qk-1)*nydy*nxdx + qi*nxdx + qj;
	  idxAbove = (qk+1)*nydy*nxdx + qi*nxdx + qj;
	  //The cell two below the current cell
	  idx2Below = (qk-2)*nydy*nxdx + qi*nxdx + qj;

	  //Calculate du/dz

	  //Cell just above the ground
	  //du/dz = Uat 1st cell/0.5dz(log(dz/znaut))
	  if(qk == 0){
	    du_dz = wind_vel[p2idx].u/(cell_dz*log(cell_dz/znaut));
	  }
	  //Cell at the top of the domain
	  else if(qk == (nzdz-1)){
	    du_dz = data3[texidxBelow];  // du_dz at k-1th cell
	  }
	  /*
	  //Cell at the top of the domain
	  //du/dz = du/dz at k-1th cell
	  //which is the same as du/dz = Uk-Uk-2/2dz if domain is higher than 2
	  else if((qk == (nz-1)) && (qk > 1)){
	    du_dz = (wind_vel[p2idx].u - wind_vel[idx2Below].u)/(2.0*dz);
	  }
	  //Cell at the top of the domain
	  //which is the same as cell just above ground if domain is only 2 high
	  else if((qk == (nz-1)) && (qk == 1)){
	    du_dz = wind_vel[idxBelow].u/(0.5*dz*log(0.5*dz/znaut));
	  }
	  */

	  //All other cells. i.e not boundary cells
	  //du/dz = (Uk+1 - Uk-1)/2dz
	  else{
	    du_dz = (wind_vel[idxAbove].u - wind_vel[idxBelow].u)/(2.0f*cell_dz);
	  }
	  	  
	  data3[texidx] = du_dz;    //du_dz
	  data3[texidx+1] = 0.0f;    //eventually dv_dz
	  data3[texidx+2] = 0.0f;    //eventually dw_dz
	  data3[texidx+3] = 0.0f;
	  
	  ustar = 0.4f*(qk+1)*du_dz;

	  sigU = 2.0f*ustar;
	  sigV = 2.0f*ustar;
	  sigW = 1.3f*ustar;

	  sig[p2idx].u = sigU;   //sigU
	  sig[p2idx].v = sigV;   //sigV
	  sig[p2idx].w = sigW;   //sigW

	  float tau11=sigU*sigU;
	  float tau22=sigV*sigV;
	  float tau33=sigW*sigW;
	  float tau13=ustar*ustar;
	  float tauDetInv= float(1.0/((tau11*tau22*tau33)-(tau13*tau13*tau22)));

	  tau[p2idx].t11   = tau11;             //Tau11
	  tau[p2idx].t22 = tau22;             //Tau22
	  tau[p2idx].t33 = tau33;             //Tau33
	  tau[p2idx].t13 = tau13;             //Tau13

	  dataTwo[texidx]   = tauDetInv*(tau22*tau33);            //Lam11
	  dataTwo[texidx+1] = tauDetInv*(tau11*tau33-tau13*tau13);//Lam22
	  dataTwo[texidx+2] = tauDetInv*(tau11*tau22);	          //Lam33
	  dataTwo[texidx+3] = tauDetInv*(-tau13*tau22);           //Lam13
        }

  createTexture(lambda, GL_RGBA32F_ARB, width,height, dataTwo);
  createTexture(duvw_dz, GL_RGBA32F_ARB, width,height, data3);
  
  delete [] data3;
  delete [] dataTwo;

  //Create the Tau/dz texture
  float tau11_dz, tau22_dz, tau33_dz, tau13_dz;
  
  for (qk=0; qk<nzdz; qk++) 
    for (qi=0; qi<nydy; qi++)
      for (qj=0; qj<nxdx; qj++)
	{
	  p2idx = qk*nydy*nxdx + qi*nxdx + qj;
	    
	  row = qk / (numInRow);
	  texidx = row * width * nydy * 4 +
	  qi * width * 4 +
	  qk % (numInRow) * nxdx * 4 +
	  qj * 4;

	  rowBelow = (qk-1) / (numInRow);
	  texidxBelow = rowBelow * width * nydy * 4 +
	    qi * width * 4 +
	    (qk-1) % (numInRow) * nxdx * 4 +
	    qj * 4;
	  
	  //The cells right above and below the current cell
	  idxBelow = (qk-1)*nydy*nxdx + qi*nxdx + qj;
	  idxAbove = (qk+1)*nydy*nxdx + qi*nxdx + qj;

	  //The cell two below the current cell
	  idx2Below = (qk-2)*nydy*nxdx + qi*nxdx + qj;

	  //Calculate Tau/dz
	  
	  //Cell just above the ground
	  //dT/dz = -2(T/(0.5dz)log((0.5dz/znaut)))
	  if(qk == 0){
	    tau11_dz = -2.0f*(tau[p2idx].t11/(cell_dz*log(cell_dz/znaut)));
	    tau22_dz = -2.0f*(tau[p2idx].t22/(cell_dz*log(cell_dz/znaut)));
	    tau33_dz = -2.0f*(tau[p2idx].t33/(cell_dz*log(cell_dz/znaut)));
	    tau13_dz = -2.0f*(tau[p2idx].t13/(cell_dz*log(cell_dz/znaut)));
	  }
	  //Cell at the top of the domain
	  else if(qk == (nzdz-1)){
	    tau11_dz = data[texidxBelow];
	    tau22_dz = data[texidxBelow+1];
	    tau33_dz = data[texidxBelow+2];
	    tau13_dz = data[texidxBelow+3];
	  }
	  /*
	  //Cell at the top of the domain
	  //dT/dz = dT/dz at k-1th cell
	  //which is the same as dT/dz = Tk-Tk-2/2dz if domain is higher than 2		     
	  else if((qk == (nzdz-1)) && qk > 1){
	    tau11_dz = (tau[p2idx].t11 - tau[idx2Below].t11)/(2.0*dz);
	    tau22_dz = (tau[p2idx].t22 - tau[idx2Below].t22)/(2.0*dz);
	    tau33_dz = (tau[p2idx].t33 - tau[idx2Below].t33)/(2.0*dz);
	    tau13_dz = (tau[p2idx].t13 - tau[idx2Below].t13)/(2.0*dz);
	  }
	  //Cell at the top of the domain
	  //which is the same as cell just above ground if domain is only 2 high
	  else if((qk == (nzdz-1)) && qk == 1){
	    tau11_dz = -2.0*(tau[idxBelow].t11/(0.5*dz*log((0.5*dz)/znaut)));
	    tau22_dz = -2.0*(tau[idxBelow].t22/(0.5*dz*log((0.5*dz)/znaut)));
	    tau33_dz = -2.0*(tau[idxBelow].t33/(0.5*dz*log((0.5*dz)/znaut)));
	    tau13_dz = -2.0*(tau[idxBelow].t13/(0.5*dz*log((0.5*dz)/znaut)));
	    }*/


	  //All other cells
	  //dT/dz = (Tk+1 - Tk-1)/2dz
	  else{
	    tau11_dz = (tau[idxAbove].t11 - tau[idxBelow].t11)/(2.0f*cell_dz);
	    tau22_dz = (tau[idxAbove].t22 - tau[idxBelow].t22)/(2.0f*cell_dz);
	    tau33_dz = (tau[idxAbove].t33 - tau[idxBelow].t33)/(2.0f*cell_dz);
	    tau13_dz = (tau[idxAbove].t13 - tau[idxBelow].t13)/(2.0f*cell_dz);
	  }

	  data[texidx] = tau11_dz;
	  data[texidx+1] = tau22_dz;
	  data[texidx+2] = tau33_dz;
	  data[texidx+3] = tau13_dz;

	}

  createTexture(tau_dz, GL_RGBA32F_ARB, width,height, data);
  
  delete [] data;

}
void ParticleControl::initLambda_and_TauTex_fromQUICFILES(GLuint windField,GLuint lambda, GLuint tau_dz, GLuint duvw_dz, 
							  GLuint tauTex){
  GLfloat *data = new GLfloat[ width * height * 4 ];
  GLfloat *dataWind = new GLfloat[ width * height * 4 ];
  GLfloat *dataTwo = new GLfloat[ width * height * 4 ];
  GLfloat *dataTau = new GLfloat[width*height*4];

  GLfloat *data3 = new GLfloat[width*height*4];
  //GLfloat *data4 = new GLfloat[width*height*4];
  //GLfloat *data5 = new GLfloat[width*height*4];

  int qi, qj, qk;
  int p2idx = 0, texidx = 0, texidxBelow = 0;
  int row = 0, rowBelow = 0;
  
  //cell below, cell above, cell two below
  int idxBelow, idxAbove, idx2Below, idxFront, idxBehind;
  //int idxRight, idxLeft;
  //float du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,dw_dx,dw_dy,dw_dz;
  //double S11,S22,S33,S21,S12,S31,S13,S23,S32;
  
  //Cell height
  //These variables are now set in the setting file
  //float dz = 1.0,dx=1.0,dy=1.0;

  //Surface roughness
  float znaut = 0.1f;

  //initializes the array cellQuic[]
  initCellType();
 
  //Reading turbulence data from the QUIC generated file
  std::ifstream turbulence;
  std::string path;
  if(quicFilesPath.c_str() != NULL){
    path = quicFilesPath + "QP_turbfield.dat";
  }
  else
    path = "Settings/QP_turbfield.dat";

  //turbulence.open("Settings/QP_turbfield.dat");
  turbulence.open(path.c_str());
  
  std::string header;
  
  turbulence>>header>>header>>header>>header>>header>>header>>header>>header>>header>>header;
  turbulence>>header>>header>>header>>header>>header>>header>>header>>header>>header>>header;
  turbulence>>header>>header>>header>>header>>header>>header>>header>>header>>header>>header;
  turbulence>>header>>header>>header>>header>>header>>header>>header>>header>>header>>header;

  turbulence>>header>>header>>header;
  //End reading the header, Now the values are read in the loop

  float indexVal,extraVal,elz,eps; // indexVal is for x,y and z locations , eps and elz are required for ustar calculations.

#if 0
    //
    // Create a VBO to hold the data for drawing arrows
    //
  
  GLfloat *windVBOdata = new GLfloat[ width * height * 3];  // Enough to hold the wind field direction vectors AND the center position of each cell
#endif

  for (qk=0; qk<nzdz; qk++) 
    for (qi=0; qi<nydy; qi++)
      for (qj=0; qj<nxdx; qj++)
	{	
	  sigU=0;
	  sigV=0;
	  sigW=0;

	  turbulence>>indexVal; // reading the X value from the file
	  turbulence>>indexVal; // reading the Y value from the file
	  turbulence>>indexVal; // reading the Z value from the file
             
	  turbulence>>sigU;// = 2.5*ustar;
	  turbulence>>sigV;// = 2.0*ustar;
	  turbulence>>sigW;// = 1.3*ustar;

	  turbulence>>extraVal;// this value in not required

	  turbulence>>elz;
	  turbulence>>eps;
	  turbulence>>extraVal;  //this value in not required
	  turbulence>>extraVal;  //this value in not required
	  turbulence>>extraVal;  //this value in not required

	  p2idx = qk*nydy*nxdx + qi*nxdx + qj;
	  row = qk / (numInRow);
	  texidx = row * width * nydy * 4 +
	    qi * width * 4 +
	    qk % (numInRow) * nxdx * 4 +
	    qj * 4;

	  rowBelow = (qk-1) / (numInRow);
	  texidxBelow = rowBelow * width * nydy * 4 +
	    qi * width * 4 +
	    (qk-1) % (numInRow) * nxdx * 4 +
	    qj * 4;
	  
		
	  if(cellQuic[p2idx].c != 0){ //Calculate gradients ONLY if its not a building or ground cell    
	   
	         
	    //The cells right above and below the current cell
	    idxBelow = (qk-1)*nydy*nxdx + qi*nxdx + qj;
	    idxAbove = (qk+1)*nydy*nxdx + qi*nxdx + qj;         

	    //The cells at front and behind of the current cell   
	    idxBehind = qk*nydy*nxdx + qi*nxdx + (qj-1);
	    idxFront  = qk*nydy*nxdx + qi*nxdx + (qj+1);
	    //The cell two below the current cell
	    idx2Below = (qk-2)*nydy*nxdx + qi*nxdx + qj;
            
		  	  
	    data3[texidx] = (sigU/2.5f)/elz;  //du_dz;      //du_dz
	    data3[texidx+1] = 0.0;    //dv_dz
	    data3[texidx+2] = 0.0;    //dw_dz
	    data3[texidx+3] = 0.0;
		  
	    ustar=sigU/2.5f;
	    data3[texidx+3] = ustar;//temporarily used for the ustar.

	    sig[p2idx].u = sigU;   //sigU
	    sig[p2idx].v = sigV;   //sigV
	    sig[p2idx].w = sigW;   //sigW

	    float tau11=sigU*sigU;
	    float tau22=sigV*sigV;
	    float tau33=sigW*sigW;
	    float tau13=ustar*ustar;
	    float tauDetInv=1.0f/((tau11*tau22*tau33)-(tau13*tau13*tau22));

	  
	    updateMaxandMinTaus(tau11,tau22,tau33,tau13);
			    
	    tau[p2idx].t11 = tau11;             //Tau11
	    tau[p2idx].t22 = tau22;             //Tau22
	    tau[p2idx].t33 = tau33;             //Tau33
	    tau[p2idx].t13 = tau13;             //Tau13
	    //Make tau's a texture so that they can be visualized as horizontal layers in the domain
	    dataTau[texidx] = tau11;
	    dataTau[texidx+1] = tau22;
	    dataTau[texidx+2] = tau33;
	    dataTau[texidx+3] = tau13;

	    dataTwo[texidx]   =  1.0f/(tau11-tau13*tau13/tau33);// tauDetInv*(tau22*tau33);      //Lam11
	    dataTwo[texidx+1] =  1.0f/tau22;// tauDetInv*(tau11*tau33-tau13*tau13);//Lam22
	    dataTwo[texidx+2] =  1.0f/(tau33-tau13*tau13/tau11);//tauDetInv*(tau11*tau22);	 //Lam33
	    dataTwo[texidx+3] =  -tau13/(tau11*tau33-tau13*tau13);//tauDetInv*(-tau13*tau22);    //Lam13
        
	    dataWind[texidx] = wind_vel[p2idx].u;
	    dataWind[texidx+1] = wind_vel[p2idx].v;
	    dataWind[texidx+2] = wind_vel[p2idx].w;	  
	    dataWind[texidx+3] = (0.5f*5.7f)*eps;//(0.5*5.7)*(ustar*ustar*ustar)/(0.4*(minDistance));
	    //This value is the '0.5*CoEps' value	

	    // Each of these are going to be considered lines so the first will be the cell center..
	    // Store the
	    // std::cout << "(" << qj << ", " << qi << ", " << qk << ")" << std::endl;

#if 0
	      windVBOdata[texidx] = qj;
	    windVBOdata[texidx+1] = qi;
	    windVBOdata[texidx+2] = qk;


	    float windMag = sqrt(wind_vel[p2idx].u*wind_vel[p2idx].u + wind_vel[p2idx].v*wind_vel[p2idx].v + wind_vel[p2idx].w*wind_vel[p2idx].w);	 

	    windVBOdata[texidx+3] = qj + 0.5 + wind_vel[p2idx].u/windMag;
	    windVBOdata[texidx+4] = qi + 0.5 + wind_vel[p2idx].v/windMag;
	    windVBOdata[texidx+5] = qk + 0.5 + wind_vel[p2idx].w/windMag;	 
#endif

	    updateMaxandMinWindVel(dataWind[texidx],dataWind[texidx+1],dataWind[texidx+2],dataWind[texidx+3]);

	  }
	  
	  else{
	    
	    //std::cout << qk << " " << qi << " " << qj << std::endl;

	    data3[texidx] = 0.0;      //du_dz
	    data3[texidx+1] = 0.0;    //dv_dz
	    data3[texidx+2] = 0.0;    //dw_dz

	    data3[texidx+3] = 0.0;//ustar;
	    sig[p2idx].u = 0.0;   //sigU
	    sig[p2idx].v = 0.0;   //sigV
	    sig[p2idx].w = 0.0;   //sigW
	    //updateMaxandMinTaus(0.0,0.0,0.0,0.0); 
	    tau[p2idx].t11 = 0.0;
	    tau[p2idx].t22 = 0.0;
	    tau[p2idx].t33 = 0.0;
	    tau[p2idx].t13 = 0.0;
	    dataTau[texidx] = 0.0;
	    dataTau[texidx+1] = 0.0;
	    dataTau[texidx+2] = 0.0;
	    dataTau[texidx+3] = 0.0;

	    dataTwo[texidx]   = 0.0;
	    dataTwo[texidx+1] =  0.0;
	    dataTwo[texidx+2] =  0.0;
	    dataTwo[texidx+3] =  0.0;
        	   
	    dataWind[texidx] = wind_vel[p2idx].u;
	    dataWind[texidx+1] = wind_vel[p2idx].v;
	    dataWind[texidx+2] = wind_vel[p2idx].w;	  
	    dataWind[texidx+3] = (0.5f*5.7f)*eps;//This value is the '0.5*CoEps' value

#if 0
	    windVBOdata[texidx] = -1.0;
	    windVBOdata[texidx+1] = -1.0;
	    windVBOdata[texidx+2] = -1.0;
#endif

	  }
	}

  find_tauLocalMax();


#if 0
    windFieldVector_w = width;
    windFieldVector_h = height;
    glGenBuffersARB(1, &windFieldVector_vbo);
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, windFieldVector_vbo);
    glBufferDataARB(GL_ARRAY_BUFFER_ARB, width*height*4*sizeof(GLfloat), windVBOdata, GL_STATIC_DRAW_ARB);
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);

    delete [] windVBOdata;
#endif




  /*std::cout << "WindMax x = " << windMax[0] << std::endl;
  std::cout << "WindMax y = " << windMax[1] << std::endl;
  std::cout << "WindMax z = " << windMax[2] << std::endl;
  std::cout << "WindMax c = " << windMax[3] << std::endl;
  std::cout << "WindMin x = " << windMin[0] << std::endl;
  std::cout << "WindMin y = " << windMin[1] << std::endl;
  std::cout << "WindMin z = " << windMin[2] << std::endl;
  std::cout << "WindMin c = " << windMin[3] << std::endl;*/

  createTexture(lambda, GL_RGBA32F_ARB, width,height, dataTwo);
  createTexture(duvw_dz, GL_RGBA32F_ARB, width,height, data3);
  createTexture(windField, GL_RGBA32F_ARB, width, height, dataWind);
  createTexture(tauTex, GL_RGBA32F_ARB, width, height, dataTau);
  
  delete [] dataTau;
  delete [] dataWind;
  delete [] data3;
  delete [] dataTwo;

  //Create the Tau/dz texture
  float tau11_dz, tau22_dz, tau33_dz, tau13_dz;
  
  for (qk=0; qk<nzdz; qk++) 
    for (qi=0; qi<nydy; qi++)
      for (qj=0; qj<nxdx; qj++)
	{
	  p2idx = qk*nydy*nxdx + qi*nxdx + qj;
	    
	  row = qk / (numInRow);
	  texidx = row * width * nydy * 4 +
	    qi * width * 4 +
	    qk % (numInRow) * nxdx * 4 +
	    qj * 4;

	  rowBelow = (qk-1) / (numInRow);
	  texidxBelow = rowBelow * width * nydy * 4 +
	    qi * width * 4 +
	    (qk-1) % (numInRow) * nxdx * 4 +
	    qj * 4;
	  
	  //The cells right above and below the current cell
	  idxBelow = (qk-1)*nydy*nxdx + qi*nxdx + qj;
	  idxAbove = (qk+1)*nydy*nxdx + qi*nxdx + qj;

	  //The cell two below the current cell
	  idx2Below = (qk-2)*nydy*nxdx + qi*nxdx + qj;

	  //Calculate Tau/dz
	  
	  //Cell just above the ground
	  //dT/dz = -2(T/(0.5dz)log((0.5dz/znaut)))
	  if(qk == 0 || (cellQuic[idxBelow].c==0 && cellQuic[p2idx].c==1)){
	    tau11_dz = -2.0f*(tau[p2idx].t11/(cell_dz*log(cell_dz/znaut)));
	    tau22_dz = -2.0f*(tau[p2idx].t22/(cell_dz*log(cell_dz/znaut)));
	    tau33_dz = -2.0f*(tau[p2idx].t33/(cell_dz*log(cell_dz/znaut)));
	    tau13_dz = -2.0f*(tau[p2idx].t13/(cell_dz*log(cell_dz/znaut)));
	  }
	  //Cell at the top of the domain
	  else if(qk == (nzdz-1)){
	    tau11_dz = data[texidxBelow];
	    tau22_dz = data[texidxBelow+1];
	    tau33_dz = data[texidxBelow+2];
	    tau13_dz = data[texidxBelow+3];
	  }
	  //All other cells
	  //dT/dz = (Tk+1 - Tk-1)/2dz
	  else{
	    tau11_dz = (tau[idxAbove].t11 - tau[idxBelow].t11)/(2.0f*cell_dz);
	    tau22_dz = (tau[idxAbove].t22 - tau[idxBelow].t22)/(2.0f*cell_dz);
	    tau33_dz = (tau[idxAbove].t33 - tau[idxBelow].t33)/(2.0f*cell_dz);
	    tau13_dz = (tau[idxAbove].t13 - tau[idxBelow].t13)/(2.0f*cell_dz);
	  }

	  data[texidx] = tau11_dz;
	  data[texidx+1] = tau22_dz;
	  data[texidx+2] = tau33_dz;
	  data[texidx+3] = tau13_dz;

	}

  createTexture(tau_dz, GL_RGBA32F_ARB, width,height, data);
  
  delete [] data;

}
void ParticleControl::initLambda_and_Taus_withCalculations(GLuint windField,GLuint lambda, GLuint tau_dz, GLuint duvw_dz, 
							   GLuint tauTex){

  GLfloat *data = new GLfloat[ width * height * 4 ];
  GLfloat *dataWind = new GLfloat[ width * height * 4 ];
  GLfloat *dataTwo = new GLfloat[ width * height * 4 ];
  GLfloat *dataTau = new GLfloat[width*height*4];

  GLfloat *data3 = new GLfloat[width*height*4];
  //GLfloat *data4 = new GLfloat[width*height*4];
  //GLfloat *data5 = new GLfloat[width*height*4];

  int qi, qj, qk;
  int p2idx = 0, texidx = 0, texidxBelow = 0;
  int row = 0, rowBelow = 0;
  
  //cell below, cell above, cell two below
  int idxBelow, idxAbove, idx2Below, idxFront, idxBehind;
  float dz = 1.0,dx=1.0,dy=1.0; //Hardwired** should be in input file
  //Surface roughness
  float znaut = 0.1f;//Hardwired** should be in input file

  //initializes the array cellQuic[]
  initCellType();  

  for (qk=0; qk<nzdz; qk++) 
    for (qi=0; qi<nydy; qi++)
      for (qj=0; qj<nxdx; qj++){
        
	float minDistance=0.0f;
	sigU=0;
	sigV=0;
	sigW=0;
	float tau11=0;
	float tau22=0;
	float tau33=0;
	float tau13=0;

	p2idx = qk*nydy*nxdx + qi*nxdx + qj;
	

	row = qk / (numInRow);
	texidx = row * width * nydy * 4 +
	  qi * width * 4 +
	  qk % (numInRow) * nxdx * 4 +
	  qj * 4;
                
	rowBelow = (qk-1) / (numInRow);
	texidxBelow = rowBelow * width * nydy * 4 +
	  qi * width * 4 +
	  (qk-1) % (numInRow) * nxdx * 4 +
	  qj * 4;

	if(cellQuic[p2idx].c!=0){ //Calculate gradients ONLY if it is a fluid cell    
         
	  minDistance = getMinDistance(qj,qi,qk);
	  	         
	  //The cells right above and below the current cell
	  idxBelow = (qk-1)*nydy*nxdx + qi*nxdx + qj;
	  idxAbove = (qk+1)*nydy*nxdx + qi*nxdx + qj;         
          
	  //The cells at front and behind of the current cell   
	  idxBehind = qk*nydy*nxdx + qi*nxdx + (qj-1);
	  idxFront  = qk*nydy*nxdx + qi*nxdx + (qj+1);
	  //The cell two below the current cell
	  idx2Below = (qk-2)*nydy*nxdx + qi*nxdx + qj;
            
	  //Calculate du/dz
            
	    
	  float du_dz=0.0f;
	  //Cell just above the ground

	  if(qk == 0 || (cellQuic[idxBelow].c==0 && cellQuic[p2idx].c==1)){
	    du_dz = wind_vel[p2idx].u/(0.5*dz*log(0.5*dz/znaut));
	  }
	    
	  //Cell at the top of the domain
	    
	  else if(qk == (nzdz-1)){
	    du_dz = data3[texidxBelow];  // du_dz at k-1th cell
	  }          
	  //All other cells. i.e not boundary cells
	  //du/dz = (Uk+1 - Uk-1)/2dz
	  else{
	    du_dz = (wind_vel[idxAbove].u - wind_vel[idxBelow].u)/(2.0f*dz); 
	  }

	  float VertGradFactor=pow( (1.0-(minDistance/20.0)) ,3.0/4.0); 
	  ustar=0.4*minDistance*du_dz; //Note: ustar doesn't include the vertGradFactor; sigmas do have vertGradFactor

	  data3[texidx] = du_dz; //du_dz
	  data3[texidx+1] = 0.0; //dv_dz
	  data3[texidx+2] = VertGradFactor;    //dw_dz
	  data3[texidx+3] = ustar;
         
	  sigU = 2.5*ustar*VertGradFactor;
	  sigV = 2.0*ustar*VertGradFactor;
	  sigW = 1.3*ustar*VertGradFactor;
          
	  sig[p2idx].u = sigU;   //sigU
	  sig[p2idx].v = sigV;   //sigV
	  sig[p2idx].w = sigW;   //sigW
         
	  tau11=sigU*sigU;
	  tau22=sigV*sigV;
	  tau33=sigW*sigW;
	  tau13=ustar*ustar;
	  float tauDetInv=1.0f/((tau11*tau22*tau33)-(tau13*tau13*tau22));
		
	  updateMaxandMinTaus(tau11,tau22,tau33,tau13);
			    
	  tau[p2idx].t11   = tau11;             //Tau11
	  tau[p2idx+1].t22 = tau22;             //Tau22
	  tau[p2idx+2].t33 = tau33;             //Tau33
	  tau[p2idx+3].t13 = tau13;             //Tau13
	  //Make tau's a texture so that they can be visualized as horizontal layers in the domain
	  dataTau[texidx] = tau11;
	  dataTau[texidx+1] = tau22;
	  dataTau[texidx+2] = tau33;
	  dataTau[texidx+3] = tau13;
         
	  dataTwo[texidx]   =  1.0f/(tau11-tau13*tau13/tau33);// tauDetInv*(tau22*tau33);   //Lam11
	  dataTwo[texidx+1] =  1.0f/tau22;// tauDetInv*(tau11*tau33-tau13*tau13);           //Lam22
	  dataTwo[texidx+2] =  1.0f/(tau33-tau13*tau13/tau11);//tauDetInv*(tau11*tau22);    //Lam33
	  dataTwo[texidx+3] =  -tau13/(tau11*tau33-tau13*tau13);//tauDetInv*(-tau13*tau22); //Lam13
         
	  dataWind[texidx] = wind_vel[p2idx].u;
	  dataWind[texidx+1] = wind_vel[p2idx].v;
	  dataWind[texidx+2] = wind_vel[p2idx].w;	  
	  dataWind[texidx+3] = (0.5*5.7)*(ustar*ustar*ustar)*
	    pow((1.0f-0.85f*minDistance/20.0f),(1.5f))/(0.4*(minDistance));//This value is the '0.5*CoEps' value	

	  updateMaxandMinWindVel(dataWind[texidx],dataWind[texidx+1],dataWind[texidx+2],dataWind[texidx+3]);

	}
	else{
	    
	  data3[texidx] = 0.0;      //du_dz
	  data3[texidx+1] = 0.0;    //dv_dz
	  data3[texidx+2] = 0.0;    //dw_dz

	  data3[texidx+3] = 0.0;//ustar;
	  sig[p2idx].u = 0.0;   //sigU
	  sig[p2idx].v = 0.0;   //sigV
	  sig[p2idx].w = 0.0;   //sigW
	  //updateMaxandMinTaus(0.0,0.0,0.0,0.0); 
	  tau[p2idx].t11 = 0.0;
	  tau[p2idx].t22 = 0.0;
	  tau[p2idx].t33 = 0.0;
	  tau[p2idx].t13 = 0.0;
	  dataTau[texidx] = 0.0;
	  dataTau[texidx+1] = 0.0;
	  dataTau[texidx+2] = 0.0;
	  dataTau[texidx+3] = 0.0;

	  dataTwo[texidx]   = 0.0;
	  dataTwo[texidx+1] =  0.0;
	  dataTwo[texidx+2] =  0.0;
	  dataTwo[texidx+3] =  0.0;
        
	  dataWind[texidx] = 0.0;
	  dataWind[texidx+1] = 0.0;
	  dataWind[texidx+2] = 0.0;	 
	  dataWind[texidx+3] = 0.0;

	  //dataWind[texidx] = wind_vel[p2idx].u;
	  //dataWind[texidx+1] = wind_vel[p2idx].v;
	  //dataWind[texidx+2] = wind_vel[p2idx].w;	  
	  //dataWind[texidx+3] = (0.5f*5.7f)*eps;//(0.5*5.7)*(ustar*ustar*ustar)/(0.4*(minDistance));

	  //This value is the '0.5*CoEps' value

	}
      }
      
  //std::cout<< (static_cast<float>(std::clock()-mytimer))/CLOCKS_PER_SEC<< " seconds" << std::endl;
     
  find_tauLocalMax();
  createTexture(lambda, GL_RGBA32F_ARB, width,height, dataTwo);
  createTexture(duvw_dz, GL_RGBA32F_ARB, width,height, data3);
  createTexture(windField, GL_RGBA32F_ARB, width, height, dataWind);
  createTexture(tauTex, GL_RGBA32F_ARB, width, height, dataTau);
  
  delete [] dataTau;
  delete [] dataWind;
  delete [] data3;
  delete [] dataTwo;

  //Create the Tau/dz texture
  float tau11_dz, tau22_dz, tau33_dz, tau13_dz;
  for (qk=0; qk<nzdz; qk++) 
    for (qi=0; qi<nydy; qi++)
      for (qj=0; qj<nxdx; qj++){
	p2idx = qk*nydy*nxdx + qi*nxdx + qj;
	    
	row = qk / (numInRow);
	texidx = row * width * nydy * 4 +
	  qi * width * 4 +
	  qk % (numInRow) * nxdx * 4 +
	  qj * 4;

	rowBelow = (qk-1) / (numInRow);
	texidxBelow = rowBelow * width * nydy * 4 +
	  qi * width * 4 +
	  (qk-1) % (numInRow) * nxdx * 4 +
	  qj * 4;
	  
	//The cells right above and below the current cell
	idxBelow = (qk-1)*nydy*nxdx + qi*nxdx + qj;
	idxAbove = (qk+1)*nydy*nxdx + qi*nxdx + qj;

	//The cell two below the current cell
	idx2Below = (qk-2)*nydy*nxdx + qi*nxdx + qj;

	//Calculate Tau/dz
	  
	//Cell just above the ground
	//dT/dz = -2(T/(0.5dz)log((0.5dz/znaut)))
	if(qk == 0 || (cellQuic[idxBelow].c==0 && cellQuic[p2idx].c==1)){
	  tau11_dz = -2.0f*(tau[p2idx].t11/(dz*log(dz/znaut)));
	  tau22_dz = -2.0f*(tau[p2idx].t22/(dz*log(dz/znaut)));
	  tau33_dz = -2.0f*(tau[p2idx].t33/(dz*log(dz/znaut)));
	  tau13_dz = -2.0f*(tau[p2idx].t13/(dz*log(dz/znaut)));
	}
	//Cell at the top of the domain
	else if(qk == (nzdz-1)){
	  tau11_dz = data[texidxBelow];
	  tau22_dz = data[texidxBelow+1];
	  tau33_dz = data[texidxBelow+2];
	  tau13_dz = data[texidxBelow+3];
	}
	//All other cells
	//dT/dz = (Tk+1 - Tk-1)/2dz
	else{
	  tau11_dz = (tau[idxAbove].t11 - tau[idxBelow].t11)/(2.0f*dz);
	  tau22_dz = (tau[idxAbove].t22 - tau[idxBelow].t22)/(2.0f*dz);
	  tau33_dz = (tau[idxAbove].t33 - tau[idxBelow].t33)/(2.0f*dz);
	  tau13_dz = (tau[idxAbove].t13 - tau[idxBelow].t13)/(2.0f*dz);
	}
	data[texidx] = tau11_dz;
	data[texidx+1] = tau22_dz;
	data[texidx+2] = tau33_dz;
	data[texidx+3] = tau13_dz;

      }

  createTexture(tau_dz, GL_RGBA32F_ARB, width,height, data);
  
  delete [] data;



}

float ParticleControl::getMinDistance(int qj, int qi, int qk){
  
  //Following is hardwired, but we need to read the building info mentioned in the input file
  int totBuild=6;
  int xfo[]={20,20,35,35,50,50};//{25};
  int yfo[]={18,32,18,32,18,32};//{25};
  int zfo[]={0,0,0,0,0,0};//{0};
  int length[]={10,10,10,10,10,10};//{10};
  int width[]={6,6,6,6,6,6};//{10};
  int height[]={6,6,10,10,14,14};//{10};
  float distance[6]; //stores minimum distance to each building

  int dz=1;
  int dy=1;
  int dx=1;
  // Adding 0.5*gridResolution as the QUIC-URB grid is shifted by 0.5*gridResolution 
  float jCell=qj+0.5*dx; // converting units in meters (original position of the cell in meters)
  float iCell=qi+0.5*dy;
  float kCell=qk+0.5*dz;


  for(int build=0;build<totBuild;++build){
  
    int x=xfo[build]; //storing the building parameters 
    int y=yfo[build];
    int z=zfo[build];
    int l=length[build];
    int w=width[build];
    int h=height[build];

 

    float minDisFaces=0.0; //minimum distance to 4 faces(sides) from a cell
    float minDisTop=0.0;//minimum distance to the top(roof) of the building from a cell 
  
    if(kCell<h){//For this condition we have only 4 planes each building
          
      float actualDis[4];// absolute value of the perpendDis or the actual distance, we have 4 faces
    
      for(int i=0;i<4;++i){
	// i=0 is front face 
	// i=1, back face
	// i=2, right side face (face towards front face of building)
	// i=3, left side face:

	float iedge;  // edges of the suface, declared as floats as ...
	float jedge;  // one of the edge value for cells perpendicular...
	//float kedge;  // to faces can be float 

	if(i==0 ){//front face
	  int edge1=y-(w/2);//right edge of the front plane
	  int edge2=y+(w/2);//left edge of the front plane
	  jedge=x;// to get the edge in X-Direction
	  if( iCell<=edge1 || iCell>=edge2 ){//for cells (qj,qi,qk) off the plane			  
	    if(fabs(edge2-iCell)< fabs(edge1-iCell))//for cells which are closer to "edge2"	   	   		  
	      iedge=edge2;
	    else
	      iedge=edge1;
		  
	  }
	  else{ //for cells perpendicular to the faces
	    iedge=iCell;
	  }
	  actualDis[i]=pow( (pow((jCell-jedge),2.0f)) + (pow((iCell-iedge),2.0f)) , 0.5f );
		
		
	}// if condition for i==0 ends
            
	if(i==1){//back face
	  int edge1=y-(w/2);
	  int edge2=y+(w/2);
	  jedge=x+l; //back face
	  if(iCell<edge1 || iCell>edge2){ 
	    if(fabs(edge2-iCell)< fabs(edge1-iCell)) //for cells which are closer to "edge2"	
	      iedge=edge2;
	    else
	      iedge=edge1;           
	  }
	  else{
	    iedge=iCell;
	  }
	  actualDis[i]=pow( (pow((jCell-jedge),2.0f)) + (pow((iCell-iedge),2.0f)) , 0.5f );
	}//if condition for i==1 ends

	if(i==2){//right side face
	  int edge1=x;
	  int edge2=x+l;
	  iedge=y-(w/2);
	  if(jCell>edge2 || jCell<edge1){
	    if(fabs(edge1-jCell) < fabs(edge2-jCell))
	      jedge=edge1;
	    else
	      jedge=edge2;	   
	  }
	  else{
	    jedge=jCell;
	  }
	  actualDis[i]=pow( (pow((jCell-jedge),2.0f)) + (pow((iCell-iedge),2.0f)) , 0.5f );
	}//if condition for i==2 ends
	if(i==3){// left side face
	  int edge1=x;
	  int edge2=x+l;
	  iedge=y+(w/2);
	  if(jCell>edge2 || jCell<edge1){
	    if(fabs(edge1-jCell) < fabs(edge2-jCell))
	      jedge=edge1;
	    else
	      jedge=edge2;  
	  }
	  else{
	    jedge=jCell;
	  }
	  actualDis[i]=pow( (pow((jCell-jedge),2.0f)) + (pow((iCell-iedge),2.0f)) , 0.5f );
	}//if condition for i==3 ends
      }// For Loop for number of faces ends
      minDisFaces=actualDis[1];//assuming one is minimum
       
      for(int i=0;i<(int)(sizeof(actualDis)/sizeof(*actualDis));++i){  //sizeof() provide number of bytes
	if(minDisFaces>actualDis[i])
	  minDisFaces=actualDis[i]; 	   
      }

      if(minDisFaces>kCell) // checking if ground is closer than any of the faces
        minDisFaces=kCell;
     
      distance[build]=minDisFaces;
    }
    else{ //if qk>=h
       
      float iedge;
      float jedge = 0.0;
      float kedge;
      
      int edgeX1=x;
      int edgeX2=x+l;
      int edgeY1=y-(w/2);
      int edgeY2=y+(w/2);
         
      if((jCell<edgeX1 || jCell>edgeX2 || iCell<edgeY1 || iCell>edgeY2)  ) { // for all the off plane cells (areas B0 and B1 in the PPT)
        iedge=iCell;
	kedge=h;
	if(jCell<=edgeX1){ // cells in front of front face
	  jedge=edgeX1;
	  if(iCell<edgeY1)
	    iedge=edgeY1;
	      
	  if(iCell>edgeY2)
	    iedge=edgeY2;   
        }
        if(jCell>=edgeX2){//cells behind the back face
          jedge=edgeX2;
	  if(iCell<=edgeY1)
	    iedge=edgeY1;
	  if(iCell>edgeY2)
	    iedge=edgeY2;    
	}
	if(jCell>edgeX1 && jCell<edgeX2){ //cells  on either side of side faces

	  jedge=jCell;
	  kedge=h;
	  if(iCell<=edgeY1)
	    iedge=edgeY1;
	  if(iCell>edgeY2)
	    iedge=edgeY2;
	}	
	
      }
      else{//if the prependicular from the cell lies on the roof.
	iedge=iCell;
	jedge=jCell;
	kedge=h;
      }	  
	  
      minDisTop=pow( (pow((jCell-jedge),2.0f)) + (pow((iCell-iedge),2.0f)) + (pow((kCell-kedge),2.0f))  , 0.5f );	
      if(minDisTop>kCell) // checking if ground is closer than the distance to the roof.
	minDisTop=kCell;
    
      distance[build]=minDisTop;

    }//if else of qk>h or qk<h ends
  }//For loop for buildings

  float minDisAll=distance[1];//assuming one is minimum
       
  for(int i=0;i<(int)(sizeof(distance)/sizeof(*distance));++i){  //sizeof() provide number of bytes
    if(minDisAll>distance[i])
      minDisAll=distance[i]; 	   
  }

  return minDisAll;
 
}


void ParticleControl::initLambdaTex(GLuint lambda){

  GLfloat *data = new GLfloat[ width * height * 4 ];

  //Create lambda texture. --Balli (04/12/07)
  float tau11=sigU*sigU;
  float tau22=sigV*sigV;
  float tau33=sigW*sigW;
  float tau13=ustar*ustar;
  float tauDetInv=1/((tau11*tau22*tau33)-(tau13*tau13*tau22));

  int qi, qj, qk;
  int p2idx = 0, texidx = 0;
  int row = 0;
  
 
  for (qk=0; qk<nzdz; qk++) 
    for (qi=0; qi<nydy; qi++)
      for (qj=0; qj<nxdx; qj++)
	{
	  p2idx = qk*nydy*nxdx + qi*nxdx + qj;
	    
	  row = qk / (numInRow);
	  texidx = row * width * nydy * 4 +
	  qi * width * 4 +
	  qk % (numInRow) * nxdx * 4 +
	  qj * 4;

	  sig[p2idx].u = sigU;   //sigU
	  sig[p2idx].v = sigV;   //sigV
	  sig[p2idx].w = sigW;   //sigW
	  sig[p2idx].id = -1.0;
	  
	  data[texidx]   = tauDetInv*(tau22*tau33);            //Lam11
	  data[texidx+1] = tauDetInv*(tau11*tau33-tau13*tau13);//Lam22
	  data[texidx+2] = tauDetInv*(tau11*tau22);	       //Lam33
	  data[texidx+3] = tauDetInv*(-tau13*tau22);           //Lam13
	}
  createTexture(lambda, GL_RGBA32F_ARB, width,height, data);
  //Lambda Texture Ends-- Balli (04/12/07)
  delete [] data;
  
}

void ParticleControl::test1(){
  //int center = (nz/2)*nx*ny + (ny/2)*nx + (nx/2);

  for(int k = 0; k < nzdz; k++){   
    for(int i = 0; i < nydy; i++){
      for(int j = 0; j < nxdx; j++){
	int p2idx = k*nxdx*nydy + i*nxdx + j;
  	
	if((j > (nxdx-1)/2) && (i > (nydy-1)/2) && (k <= ((nzdz/2)+ 1)) && (k >= ((nzdz/2)- 1))){
	  wind_vel[p2idx].u = 1.0;
	  wind_vel[p2idx].v = 0.0;
	  wind_vel[p2idx].w = 0.0;
	  wind_vel[p2idx].id = -1.0;
	}
	else if ((j <= (nxdx-1)/2) && ( i <= (nydy-1)/2) && (k <= ((nzdz/2)+ 1)) && (k >= ((nzdz/2)- 1))){
	  wind_vel[p2idx].u = -1.0;
	  wind_vel[p2idx].v = 0.0;
	  wind_vel[p2idx].w = 0.0;
	  wind_vel[p2idx].id = -1.0;
	}
	else if((i > (nydy-1)/2) && (j <= (nxdx-1)/2) && (k <= ((nzdz/2)+ 1)) && (k >= ((nzdz/2)- 1))){
	  wind_vel[p2idx].u = 0.0;
	  wind_vel[p2idx].v = 1.0;
	  wind_vel[p2idx].w = 0.0;
	  wind_vel[p2idx].id = -1.0;
	}
	else if((j > (nxdx-1)/2) && (i <= (nydy-1)/2) && (k <= ((nzdz/2)+ 1)) && (k >= ((nzdz/2)- 1))){
	  wind_vel[p2idx].u = 0.0;
	  wind_vel[p2idx].v = -1.0;
	  wind_vel[p2idx].w = 0.0;
	  wind_vel[p2idx].id = -1.0;
	  
	}
	else if( k > ((nzdz/2) +1) ){
	  wind_vel[p2idx].u = 0.0;
	  wind_vel[p2idx].v = 0.0;
	  wind_vel[p2idx].w = 1.0;
	  wind_vel[p2idx].id = -1.0;
	}
	else{
	  wind_vel[p2idx].u = 0.0;
	  wind_vel[p2idx].v = 0.0;
	  wind_vel[p2idx].w = -1.0;
	  wind_vel[p2idx].id = -1.0;
	}
      }
    }
  }
}
//Creates a random value wind field.
void ParticleControl::randomWindField(){
  for(int k = 0; k < nzdz; k++){   
    for(int i = 0; i < nydy; i++){
      for(int j = 0; j < nxdx; j++){
	int p2idx = k*nxdx*nydy + i*nxdx + j;
	//
	// Not currently random
	// 	
	wind_vel[p2idx].u = 0.0;//randVal();
	wind_vel[p2idx].v = 1.0;//randVal();
	wind_vel[p2idx].w = 0.0;//randVal();
	wind_vel[p2idx].id = -1.0;
      }
    }
  }
}

void ParticleControl::uniformUWindField(){
 
  for(int k = 0; k < nzdz; k++){   
    for(int i = 0; i < nydy; i++){
      for(int j = 0; j < nxdx; j++){
	int p2idx = k*nxdx*nydy + i*nxdx + j;
	wind_vel[p2idx].u = 1.0;
	wind_vel[p2idx].v = 0.0;
	wind_vel[p2idx].w = 0.0;
	wind_vel[p2idx].id = -1.0;
      }
    }
  }
}
void ParticleControl::variedUWindField(){
   for(int k = 0; k < nzdz; k++){   
    for(int i = 0; i < nydy; i++){
      for(int j = 0; j < nxdx; j++){
	int p2idx = k*nxdx*nydy + i*nxdx + j;
	wind_vel[p2idx].u = 7.52f*pow(((k+1)/20.0f),0.15f);
	wind_vel[p2idx].v = 0.0f;
	wind_vel[p2idx].w = 0.0f;
	//wind_vel[p2idx].id = -1.0;
      }
    }
   }
   /*if(xfo != NULL)
     addBuildingsInWindField();
   else
   std::cout << "NO BUILDINGS ADDED TO WIND FIELD" << std::endl;*/
}
void ParticleControl::addBuildingsInWindField(GLuint cellType){

  GLfloat* data = new GLfloat[width*height*4];
  wind* cell_type = new wind[nxdx*nydy*nzdz];



  for(int k = 0; k < nzdz; k++){   
    for(int i = 0; i < nydy; i++){
      for(int j = 0; j < nxdx; j++){
	int p2idx = k*nxdx*nydy + i*nxdx + j;
	cell_type[p2idx].u = 1.0;
	cell_type[p2idx].v = 1.0;
	cell_type[p2idx].w = 1.0;
	cell_type[p2idx].id = -1.0;
      }
    }
   }
  //// READ A LINE
  // CONVERT POS TO INDEX
  // IF CELL TYPE == 0 THEN
  //     U,V,W = 1.0
  //       ID = -1.0
  // ELSE 
  //     U,V,W = 0.0
  //    ID = CELL_TYPE VAL
  std::string path;
  if(quicFilesPath.c_str() != NULL){
    path = quicFilesPath + "QU_celltype.dat";
  }
  else
    path = "Settings/QU_celltype.dat";
  
  std::ifstream in;
  in.open(path.c_str(),std::ios::in);
    	
  if(in == NULL) std::cout << "input file didn't open" << std::endl;
    
  char line[1024];

  while(  !in.eof() )
  {  
	 in.getline(line, 1024);
	 if( line[ strlen(line)] == '\n' ){
		   line[ strlen(line)] = '\0';
	  }
          
          std::istringstream ist(line);
	  float x1,y1,z1,id1;
	  ist >> x1 >> y1 >> z1 >> id1;
          float x=x1;
          float y=y1;
          float z=z1;
          int id=int(id1);
          if(id==0 && z>0)
          {
            
	     int p2idx = int(z)*nxdx*nydy + int(y)*nxdx + int(x);
             cell_type[p2idx].u = 0.0;
	     cell_type[p2idx].v = 0.0;
	     cell_type[p2idx].w = 0.0;
              for(int n=0; n < numBuild; n++)
              {
                  /*if(numSides[n]==4)
                  {
                     float lk = zfo[n];
                     float uk = zfo[n]+ht[n];
                     float li = yfo[n]-(wti[n]/2.0);
                     float ui = yfo[n]+(wti[n]/2.0);
                     float lj = xfo[n];
                     float uj = xfo[n]+lti[n];

                     if(x >= lj && x <= uj && y >= li && y <= ui && z >= lk && z <= uk) cell_type[p2idx].id = float(n);
                      
                  }*/
                  if(numSides[n]==5)
                  {
                      float lk = zfo[n];
                      float uk = zfo[n]+ht[n];
                      float j = yfo[n];
                      float i = xfo[n];
                     if(z >= lk && z<= uk)
                     {
                         float xi[5],yi[5];
                         float cp1[5];
                         int count=0;
                         xi[0]=xfo[n]-(lti[n]/3.0);
                         yi[0]=yfo[n]-(wti[n]/2.0);
                         xi[1]=xfo[n]+(lti[n]/3.0);
                         yi[1]=yfo[n]-(wti[n]/2.0);
                         xi[2]=xfo[n]+(lti[n]/2.0);
                         yi[2]=yfo[n]+(wti[n]/6.0);
                         xi[3]=xfo[n];
                         yi[3]=yfo[n]+(wti[n]/2.0);
                         xi[4]=xfo[n]-(lti[n]/2.0);
                         yi[4]=yfo[n]+(wti[n]/6.0);

                         for(int m=0;m<5;m++)
                         {
                              if(m<4)
                              {
                                cp1[m]=(((xi[m+1]-xi[m])*(j-yi[m]))*((xi[m+1]-xi[m])*(y-yi[m])))+(((yi[m+1]-yi[m])*(i-xi[m]))*((yi[m+1]-yi[m])*(x-xi[m])));
                              }
                              else
                              {
                                cp1[m]=(((xi[0]-xi[m])*(j-yi[m]))*((xi[0]-xi[m])*(y-yi[m])))+(((yi[0]-yi[m])*(i-xi[m]))*((yi[0]-yi[m])*(x-xi[m])));
                              
                              }
                              if(cp1[m]>=0)count++;
                         }
                         if(count==5){ cell_type[p2idx].id = float(n); //printf("%f %f %f %d %d\n",x,y,z,p2idx,n);}

                     }
                        
                  }
                  else if(numSides[n]==1)
                  {
                     float lk = zfo[n];
                     float uk = zfo[n]+ht[n];
                     float j = yfo[n];
                     float r = lti[n]/2.0;
                     float i = xfo[n]+(lti[n]/2.0);
                     if(z >= lk && z<= uk)
                     {
                          float dist=sqrt((pow(x-i,2)+pow(y-j,2)));
                          if(dist<=r) cell_type[p2idx].id = float(n);
                     }
                  }
                  
              }  
          }

  }
  }  
  in.close();


 /* for(int k = 0; k < nzdz; k++){   
    for(int i = 0; i < nydy; i++){
      for(int j = 0; j < nxdx; j++){
	int p2idx = k*nxdx*nydy + i*nxdx + j;
	cell_type[p2idx].u = 1.0;
	cell_type[p2idx].v = 1.0;
	cell_type[p2idx].w = 1.0;
	cell_type[p2idx].id = -1.0;
      }
    }
   }*/
  

 /*
  for(int nb=0;nb<numBuild;nb++)
  {
      /* if(numSides[nb]==4)
       {
          for(int hr=int(zfo[nb]);hr<=int(zfo[nb]+ht[nb]);hr++)
          {
           for(int wr=0;wr<=wti[nb];wr++)
           {
                int p2idx = int(zfo[nb]+hr)*nxdx*nydy + int((yfo[nb]- (wti[nb]/2.0))+wr)*nydy + int(xfo[nb]); 
                cell_type[p2idx].u = 0.0;
	        cell_type[p2idx].v = 0.0;
	        cell_type[p2idx].w = 0.0;
                cell_type[p2idx].id = float(nb);
           }
           for(int wr=0;wr<=wti[nb];wr++)
           {
                int p2idx = int(zfo[nb]+hr)*nxdx*nydy + int((yfo[nb]- (wti[nb]/2.0))+wr)*nydy + int(xfo[nb]+lti[nb]); 
                cell_type[p2idx].u = 0.0;
	        cell_type[p2idx].v = 0.0;
	        cell_type[p2idx].w = 0.0;
                cell_type[p2idx].id = float(nb);       
           }
           for(int lr=0;lr<=lti[nb];lr++)
           {
                int p2idx = int(zfo[nb]+hr)*nxdx*nydy + int(yfo[nb]- (wti[nb]/2.0))*nydy + int(xfo[nb]+lr); 
                cell_type[p2idx].u = 0.0;
	        cell_type[p2idx].v = 0.0;
	        cell_type[p2idx].w = 0.0;
                cell_type[p2idx].id = float(nb);
           }
           for(int lr=0;lr<=lti[nb];lr++)
           {
                int p2idx = int(zfo[nb]+hr)*nxdx*nydy + int((yfo[nb]+ (wti[nb]/2.0)))*nydy + int(xfo[nb]+lr); 
                cell_type[p2idx].u = 0.0;
	        cell_type[p2idx].v = 0.0;
	        cell_type[p2idx].w = 0.0;
                cell_type[p2idx].id = float(nb);
           }
          }
       }
       else*//*
       if(numSides[nb]==5)
       {

           float s1,s2,s3,s4,s5,c1,c2,c3,c4,c5;
           s1=(sin(-54.0*3.14/180)-sin(18.0*3.14/180))/(cos(-54.0*3.14/180)-cos(18.0*3.14/180));
           s2=(sin(-126.0*3.14/180)-sin(-54.0*3.14/180))/(cos(-126.0*3.14/180)-cos(-54.0*3.14/180));
           s3=(sin(162.0*3.14/180)-sin(-126.0*3.14/180))/(cos(162.0*3.14/180)-cos(-126.0*3.14/180));
           s4=(sin(90.0*3.14/180)-sin(162.0*3.14/180))/(cos(90.0*3.14/180)-cos(162.0*3.14/180));
           s5=(sin(18.0*3.14/180)-sin(90.0*3.14/180))/(cos(18.0*3.14/180)-cos(90.0*3.14/180));
           c1=(yfo[nb]+((wti[nb]/2.0)*sin((18.0*3.14)/180.0)))-s1*(xfo[nb]+((wti[nb]/2.0)*cos((18.0*3.14)/180.0)));
           c2=(yfo[nb]+((wti[nb]/2.0)*sin((-54.0*3.14)/180.0)))-s2*(xfo[nb]+((wti[nb]/2.0)*cos((-54.0*3.14)/180.0)));
           c3=(yfo[nb]+((wti[nb]/2.0)*sin((-126.0*3.14)/180.0)))-s3*(xfo[nb]+((wti[nb]/2.0)*cos((-126.0*3.14)/180.0)));
           c4=(yfo[nb]+((wti[nb]/2.0)*sin((162.0*3.14)/180.0)))-s4*(xfo[nb]+((wti[nb]/2.0)*cos((162.0*3.14)/180.0)));
           c5=(yfo[nb]+((wti[nb]/2.0)*sin((90.0*3.14)/180.0)))-s5*(xfo[nb]+((wti[nb]/2.0)*cos((90.0*3.14)/180.0)));
           
           for(int h=int(zfo[nb]);h<=int(zfo[nb]+ht[nb]);h++)
           {
              //printf("%f %f \n",xfo[nb]+((wti[nb]/2.0)*cos((-54.0*3.14)/180.0)),yfo[nb]+((wti[nb]/2.0)*sin((-54.0*3.14)/180.0)));
              for(float x1=xfo[nb]+(wti[nb]/2.0)*cos(18.0*3.14/180);x1>=xfo[nb]+(wti[nb]/2.0)*cos(-54.0*3.14/180);x1=x1-0.2)
              {
                 float y1=x1*s1+c1;
                 int p2idx = int(zfo[nb]+h)*nxdx*nydy + int(y1)*nydy + int(x1);
                 cell_type[p2idx].u = 0.0;
	         cell_type[p2idx].v = 0.0;
	         cell_type[p2idx].w = 0.0; 
                 cell_type[p2idx].id = float(nb);
                 //printf("%f %f \n",xfo[nb]+(wti[nb]/2.0)*cos((18.0*3.14)/180.0),xfo[nb]+(wti[nb]/2.0)*sin((-54.0*3.14)/180.0));
                 //printf("%f %f \n",x1,y1);
              }
              for(float x1=xfo[nb]+(wti[nb]/2.0)*cos(-54.0*3.14/180);x1>=xfo[nb]+(wti[nb]/2.0)*cos(-126.0*3.14/180);x1=x1-0.2)
              {
                 float y1=x1*s2+c2;
                  int p2idx = int(zfo[nb]+h)*nxdx*nydy + int(y1)*nydy + int(x1);
                 cell_type[p2idx].u = 0.0;
	         cell_type[p2idx].v = 0.0;
	         cell_type[p2idx].w = 0.0; 
                 cell_type[p2idx].id = float(nb);
                 //printf("%f %f \n",xfo[nb]+(wti[nb]/2.0)*cos((18.0*3.14)/180.0),xfo[nb]+(wti[nb]/2.0)*sin((-54.0*3.14)/180.0));
                 //printf("%f %f %d \n",x1,y1,h);
              }
              for(float x1=xfo[nb]+(wti[nb]/2.0)*cos(-126.0*3.14/180);x1>=xfo[nb]+(wti[nb]/2.0)*cos(162.0*3.14/180);x1=x1-0.2)
              {
                 float y1=x1*s3+c3;
                  int p2idx = int(zfo[nb]+h)*nxdx*nydy + int(y1)*nydy + int(x1);
                 cell_type[p2idx].u = 0.0;
	         cell_type[p2idx].v = 0.0;
	         cell_type[p2idx].w = 0.0; 
                 cell_type[p2idx].id = float(nb);
                 //printf("%f %f \n",xfo[nb]+(wti[nb]/2.0)*cos((18.0*3.14)/180.0),xfo[nb]+(wti[nb]/2.0)*sin((-54.0*3.14)/180.0));
                 //printf("%f %f %d \n",x1,y1,h);
              }
             for(float x1=xfo[nb]+(wti[nb]/2.0)*cos(162.0*3.14/180);x1<=xfo[nb]+(wti[nb]/2.0)*cos(90.0*3.14/180);x1=x1+0.2)
              {
                 float y1=x1*s4+c4;
                  int p2idx = int(zfo[nb]+h)*nxdx*nydy + int(y1)*nydy + int(x1);
                 cell_type[p2idx].u = 0.0;
	         cell_type[p2idx].v = 0.0;
	         cell_type[p2idx].w = 0.0; 
                 cell_type[p2idx].id = float(nb);
                 //printf("%f %f \n",xfo[nb]+(wti[nb]/2.0)*cos((18.0*3.14)/180.0),xfo[nb]+(wti[nb]/2.0)*sin((-54.0*3.14)/180.0));
                 //printf("%f %f %d \n",x1,y1,h);
              }
              for(float x1=xfo[nb]+(wti[nb]/2.0)*cos(90.0*3.14/180);x1<=xfo[nb]+(wti[nb]/2.0)*cos(18.0*3.14/180);x1=x1+0.2)
              {
                 float y1=x1*s5+c5;
                 int p2idx = int(zfo[nb]+h)*nxdx*nydy + int(y1)*nydy + int(x1);
                 cell_type[p2idx].u = 0.0;
	         cell_type[p2idx].v = 0.0;
	         cell_type[p2idx].w = 0.0; 
                 cell_type[p2idx].id = float(nb);
                 //printf("%f %f \n",xfo[nb]+(wti[nb]/2.0)*cos((18.0*3.14)/180.0),xfo[nb]+(wti[nb]/2.0)*sin((-54.0*3.14)/180.0));
                
              }
           }
       }
       else if(numSides[nb]==1)
       {
           float incr=1.0;
           for(int h=int(zfo[nb]);h<=int(zfo[nb]+ht[nb]);h++)
           {
              for(float Theta=0.0;Theta<=360.0;Theta+=incr)
              {
                
                    float x1=xfo[nb]-((lti[nb]/2.0)*cos(Theta*(3.14/180.0))-lti[nb]/2.0);
                    float y1=yfo[nb]-((lti[nb]/2.0)*sin(Theta*(3.14/180.0)));
                    int p2idx = int(zfo[nb]+h)*nxdx*nydy + int(y1)*nydy + int(x1);
                    cell_type[p2idx].u = 0.0;
	            cell_type[p2idx].v = 0.0;
	            cell_type[p2idx].w = 0.0; 
                    cell_type[p2idx].id = float(nb);
	      }         
           }          
       }
       
  }*/
 for(int n=0; n < numBuild; n++){
   if(numSides[n]==4)
   {
    int lk = int(zfo[n]);
    int uk = int(zfo[n]+ht[n]);
    int li = int(yfo[n]-(wti[n]/2.0));
    int ui = int(yfo[n]+(wti[n]/2.0));
    int lj = int(xfo[n]);
    int uj = int(xfo[n]+lti[n]);

     for(int k= lk; k < uk; k++){
       for(int i= li; i < ui; i++){
	 for(int j= lj; j < uj; j++){
	   int p2idx = k*nxdx*nydy + i*nxdx + j;
	   
	   cell_type[p2idx].u = 0.0;
	   cell_type[p2idx].v = 0.0;
	   cell_type[p2idx].w = 0.0;
	   cell_type[p2idx].id = float(n);

	 }
       }
     }
    }
   }
   int qk,qi,qj;
   int row,texidx,p2idx;

   for (qk=0; qk<nzdz; qk++) 
      for (qi=0; qi<nydy; qi++)
	for (qj=0; qj<nxdx; qj++)
	  {
	    p2idx = qk*nydy*nxdx + qi*nxdx + qj;
	    
	    row = qk / (numInRow);
	    texidx = row * width * nydy * 4 +
	      qi * width * 4 +
	      qk % (numInRow) * nxdx * 4 +
	      qj * 4;
	  
	    data[texidx] = cell_type[p2idx].u;
	    data[texidx+1] = cell_type[p2idx].v;
	    data[texidx+2] = cell_type[p2idx].w;	    
	    data[texidx+3] = cell_type[p2idx].id;
	    
	  }

    createTexture(cellType, GL_RGBA16F_ARB, width, height, data);

    delete [] data;
    delete [] cell_type;

}
void ParticleControl::QUICWindField(){
  std::ifstream QUICWindField;//,QUICCellType;
	
  std::string path;
  if(quicFilesPath.c_str() != NULL){
    path = quicFilesPath + "QU_velocity.dat";
  }
  else
    path = "Settings/QU_velocity.dat";

  //QUICWindField.open("Settings/QU_velocity.dat"); //opening the wind file  to read
  QUICWindField.open(path.c_str()); //opening the wind file  to read
  if(!QUICWindField){
    std::cerr<<"Unable to open QUIC Windfield file : QU_velocity.dat ";
    exit(1);
  }

  std::string header;  //I am just using a very crude method to read the header of the wind file
  QUICWindField>>header>>header>>header>>header>>header>>header>>header>>header>>header>>header>>header>>header>>header>>header>>header>>header;
  QUICWindField>>header>>header>>header>>header>>header;
    
  double groundVal; // ignoring the ground values, so that k=0 have the first cell having non-zero velocity 
  //Balli had ++k?

  for(int k=0;k<(6*nxdx*nydy);++k){ // there are 6 columns in the wind file 
    QUICWindField>>groundVal;
  }

  double quicIndex;

  for(int k = 0; k < nzdz; k++){   
    for(int i = 0; i < nydy; i++){
      for(int j = 0; j < nxdx; j++){
	int p2idx = k*nxdx*nydy + i*nxdx + j;
	QUICWindField>>quicIndex; // ignoring the X,Y and Z values
	QUICWindField>>quicIndex;
	QUICWindField>>quicIndex;

	QUICWindField>>wind_vel[p2idx].u ;//storing the velocity values in the wind structure
	QUICWindField>>wind_vel[p2idx].v ;
	QUICWindField>>wind_vel[p2idx].w ;

	/*QUICCellType>>quicIndex;// ignoring the X,Y and Z values
	QUICCellType>>quicIndex;
	QUICCellType>>quicIndex;

	QUICCellType>>cellQuic[p2idx].c ;//storing the Celltype values in the Cell structure*/
		
      }
    }
  }
  
  QUICWindField.close();
}
void ParticleControl::initCellType(){
  std::ifstream QUICCellType;
  std::string path;

  if(quicFilesPath.c_str() != NULL){
    path = quicFilesPath + "QU_celltype.dat";
  }
  else
    path = "Settings/QU_celltype.dat";
  
  QUICCellType.open(path.c_str());//opening the Celltype file  to read
  //QUICCellType.open("Settings/QU_celltype.dat");//opening the Celltype file  to read
  if(!QUICCellType.is_open()){
    std::cerr<<"Unable to open QUIC Celltype file : QU_celltype.dat ";
    exit(1);
  }
  double groundVal;
  // ignoring the ground values for the cellype also
  //Balli had ++k ?

  for(int k=0;k<(4*nxdx*nydy);++k){// there are 4 columns in the wind file
    QUICCellType>>groundVal;
  }
  double quicIndex;
  for(int k = 0; k < nzdz; k++){   
    for(int i = 0; i < nydy; i++){
      for(int j = 0; j < nxdx; j++){
	int p2idx = k*nxdx*nydy + i*nxdx + j;
	
	QUICCellType>>quicIndex;// ignoring the X,Y and Z values
	QUICCellType>>quicIndex;
	QUICCellType>>quicIndex;

	QUICCellType>>cellQuic[p2idx].c ;//storing the Celltype values in the Cell structure

      }
    }
  }
  
  QUICCellType.close();
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
  glUniform1iARB(u_nx, nxdx);
  glUniform1iARB(u_ny, nydy);
  glUniform1iARB(u_nz, nzdz);

  glBegin(GL_QUADS);
  {
    glTexCoord2f(0, 0);								glVertex3f(-1, -1, -0.5f);
    glTexCoord2f(float(twidth), 0);					glVertex3f( 1, -1, -0.5f);
    glTexCoord2f(float(twidth), float(theight));	glVertex3f( 1,  1, -0.5f);
    glTexCoord2f(0, float(theight));		glVertex3f(-1,  1, -0.5f);
  }
  glEnd();

  init_shader.deactivate();
  glDisable(texType);
  FramebufferObject::Disable();
  glViewport(vp[0], vp[1], vp[2], vp[3]);
  glDrawBuffer(draw_buffer);

}
void ParticleControl::updateMaxandMinWindVel(float x, float y, float z, float c){
  if(x > windMax[0])
    windMax[0] = x;
  if(y > windMax[1])
    windMax[1] = y;
  if(z > windMax[2])
    windMax[2] = z;
  if(c > windMax[3])
    windMax[3] = c;

  if(x < windMin[0])
    windMin[0] = x;
  if(y < windMin[1])
    windMin[1] = y;
  if(z < windMin[2])
    windMin[2] = z;
  if(c < windMin[3])
    windMin[3] = c;

}

void ParticleControl::updateMaxandMinTaus(float tau11,float tau22,float tau33,float tau13){

  //update Max
  if(tau11 > tauMax[0])
    tauMax[0] = tau11;
  if(tau22 > tauMax[1])
    tauMax[1] = tau22;
  if(tau33 > tauMax[2])
    tauMax[2] = tau33;
  if(tau13 > tauMax[3])
    tauMax[3] = tau13;

  //update Min
  if(tau11 < tauMin[0])
    tauMin[0] = tau11;
  if(tau22 < tauMin[1])
    tauMin[1] = tau22;
  if(tau33 < tauMin[2])
    tauMin[2] = tau33;
  if(tau13 < tauMin[3])
    tauMin[3] = tau13;

}
void ParticleControl::find_tauLocalMax(){
  tauLocalMax = new float[4*nzdz];
  tauLocalMin = new float[4*nzdz];
 
  //Initialize max and min 
  for(int k=0; k<nzdz; k++){
   
    int idx = k*nydy*nxdx;
    int tidx = k*4;

    tauLocalMax[tidx] = tau[idx].t11;
    tauLocalMax[tidx+1] = tau[idx].t22;
    tauLocalMax[tidx+2] = tau[idx].t33;
    tauLocalMax[tidx+3] = tau[idx].t13;

    
    tauLocalMin[tidx] = tau[idx].t11;
    tauLocalMin[tidx+1] = tau[idx].t22;
    tauLocalMin[tidx+2] = tau[idx].t33;
    tauLocalMin[tidx+3] = tau[idx].t13;
    
  }

  //Find max and min tau values for each height value

  for(int k=0; k<nzdz; k++){
    for(int i=0; i<nydy; i++){
      for(int j=0; j<nxdx; j++){
	int idx = k*nydy*nxdx + i*nxdx + j;
	int tidx = k*4;

	if(tau[idx].t11 > tauLocalMax[tidx])
	  tauLocalMax[tidx] = tau[idx].t11;
	if(tau[idx].t22 > tauLocalMax[tidx+1])
	  tauLocalMax[tidx+1] = tau[idx].t22;
	if(tau[idx].t33 > tauLocalMax[tidx+2])
	  tauLocalMax[tidx+2] = tau[idx].t33;
	if(tau[idx].t13 > tauLocalMax[tidx+3])
	  tauLocalMax[tidx+3] = tau[idx].t13;

	if(tau[idx].t11 < tauLocalMin[tidx])
	  tauLocalMin[tidx] = tau[idx].t11;
	if(tau[idx].t22 < tauLocalMin[tidx+1])
	  tauLocalMin[tidx+1] = tau[idx].t22;
	if(tau[idx].t33 < tauLocalMin[tidx+2])
	  tauLocalMin[tidx+2] = tau[idx].t33;
	if(tau[idx].t13 < tauLocalMin[tidx+3])
	  tauLocalMin[tidx+3] = tau[idx].t13;



      }
    }
  }

  std::cout << "Local Max t11 values" << std::endl;
  for(int i=0; i < nzdz; i++){
    //std::cout << tauLocalMax[(i*4)+3] << std::endl;
    
  }

}

