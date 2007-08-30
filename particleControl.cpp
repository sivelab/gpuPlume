#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "particleControl.h"
#include "Random.h"
#include "glErrorUtil.h"

ParticleControl::ParticleControl(GLenum type,int width,int height,
				 int x, int y, int z){

  texType = type;
  twidth = width;
  theight = height;
  nx = x;
  ny = y;
  nz = z;

  outputPrime = false;
  alreadyOpen = false;

}
void ParticleControl::setBuildingParameters(int nB,float* x,float* y,float* z,
					   float* h,float* w,float* l){
  numBuild = nB;
  xfo = x;
  yfo = y;
  zfo = z;
  ht = h;
  wti = w;
  lti = l; 
}
void ParticleControl::setUstarAndSigmas(float u){
  ustar = u;
  sigU = 2.0*ustar;
  sigV = 2.0*ustar;
  sigW = 1.3*ustar;
}
void ParticleControl::setRandomTexCoords(){
  t1 = Random::uniform() * twidth;
  t2 = Random::uniform() * theight;
}
void ParticleControl::setupMultipleBuildingsShader(int numInRow, float life_time){
  multipleBuildings_shader.addShader("Shaders/multipleBuildingsAdvect_vp.glsl", GLSLObject::VERTEX_SHADER);
  multipleBuildings_shader.addShader("Shaders/multipleBuildingsAdvect_fp.glsl", GLSLObject::FRAGMENT_SHADER);
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

  GLint ulifeTime = multipleBuildings_shader.createUniform("life_time");
  GLint unx = multipleBuildings_shader.createUniform("nx");
  GLint uny = multipleBuildings_shader.createUniform("ny");
  GLint unz = multipleBuildings_shader.createUniform("nz");
  GLint uNumInRow = multipleBuildings_shader.createUniform("numInRow");

  multipleBuildings_shader.activate();

  glUniform1fARB(ulifeTime, life_time);
  glUniform1iARB(unx, nx);
  glUniform1iARB(uny, ny);
  glUniform1iARB(unz, nz);
  glUniform1fARB(uNumInRow, numInRow);

  multipleBuildings_shader.deactivate();
}
void ParticleControl::multipleBuildingsAdvect(bool odd, GLuint windField, GLuint positions0, 
			     GLuint positions1, GLuint prime0, GLuint prime1, 
			     GLuint randomValues,GLuint lambda, GLuint tau_dz, 
			     GLuint duvw_dz, float time_step, GLuint buildings,
			     GLuint cellType)
{
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
  multipleBuildings_shader.activate();

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
    glTexCoord2f(0, 0);            glVertex3f(-1, -1, -0.5f);
    glTexCoord2f(twidth, 0);       glVertex3f( 1, -1, -0.5f);
    glTexCoord2f(twidth, theight); glVertex3f( 1,  1, -0.5f);
    glTexCoord2f(0, theight);      glVertex3f(-1,  1, -0.5f);
  }
  glEnd();
  
  multipleBuildings_shader.deactivate();
 
  glBindTexture(texType, 0);
  
  //Prints out the updated prime values
  if(outputPrime)
    printPrime(odd,false);

}
void ParticleControl::setupReflectionShader(int numInRow, float life_time){
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
  uniform_numBuild = reflection_shader.createUniform("numBuild");

  uniform_xfo = reflection_shader.createUniform("xfo");
  uniform_yfo = reflection_shader.createUniform("yfo");
  uniform_zfo = reflection_shader.createUniform("zfo");
  uniform_ht = reflection_shader.createUniform("ht");
  uniform_wti = reflection_shader.createUniform("wti");
  uniform_lti = reflection_shader.createUniform("lti");

  uniform_xfo2 = reflection_shader.createUniform("xfo2");
  uniform_yfo2 = reflection_shader.createUniform("yfo2");
  uniform_zfo2 = reflection_shader.createUniform("zfo2");
  uniform_ht2 = reflection_shader.createUniform("ht2");
  uniform_wti2 = reflection_shader.createUniform("wti2");
  uniform_lti2 = reflection_shader.createUniform("lti2");

  GLint ulifeTime = reflection_shader.createUniform("life_time");
  GLint unx = reflection_shader.createUniform("nx");
  GLint uny = reflection_shader.createUniform("ny");
  GLint unz = reflection_shader.createUniform("nz");
  GLint uNumInRow = reflection_shader.createUniform("numInRow");

  reflection_shader.activate();

  glUniform1fARB(ulifeTime, life_time);
  glUniform1iARB(unx, nx);
  glUniform1iARB(uny, ny);
  glUniform1iARB(unz, nz);
  glUniform1fARB(uNumInRow, numInRow);

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

  glUniform1fARB(uniform_numBuild, buildParam[0]);
  //Building varibles
  if(buildParam[0] > 0){
    glUniform1fARB(uniform_xfo, buildParam[1]);
    glUniform1fARB(uniform_yfo, buildParam[2]);
    glUniform1fARB(uniform_zfo, buildParam[3]);
    glUniform1fARB(uniform_ht, buildParam[4]);
    glUniform1fARB(uniform_wti, buildParam[5]);
    glUniform1fARB(uniform_lti, buildParam[6]);
  }
  if(buildParam[0] > 1){
    glUniform1fARB(uniform_xfo2, buildParam[7]);
    glUniform1fARB(uniform_yfo2, buildParam[8]);
    glUniform1fARB(uniform_zfo2, buildParam[9]);
    glUniform1fARB(uniform_ht2, buildParam[10]);
    glUniform1fARB(uniform_wti2, buildParam[11]);
    glUniform1fARB(uniform_lti2, buildParam[12]);
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
    glTexCoord2f(0, 0);            glVertex3f(-1, -1, -0.5f);
    glTexCoord2f(twidth, 0);       glVertex3f( 1, -1, -0.5f);
    glTexCoord2f(twidth, theight); glVertex3f( 1,  1, -0.5f);
    glTexCoord2f(0, theight);      glVertex3f(-1,  1, -0.5f);
  }
  glEnd();
  
  reflection_shader.deactivate();
 
  glBindTexture(texType, 0);
  
  //Prints out the updated prime values
  if(outputPrime)
    printPrime(odd,false);

}
void ParticleControl::setupNonGaussianShader(int numInRow, float life_time){
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
  GLint uNumInRow = nonGaussian_shader.createUniform("numInRow");

  nonGaussian_shader.activate();

  glUniform1fARB(ulifeTime, life_time);
  glUniform1iARB(unx, nx);
  glUniform1iARB(uny, ny);
  glUniform1iARB(unz, nz);
  glUniform1fARB(uNumInRow, numInRow);

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
    glTexCoord2f(0, 0);            glVertex3f(-1, -1, -0.5f);
    glTexCoord2f(twidth, 0);       glVertex3f( 1, -1, -0.5f);
    glTexCoord2f(twidth, theight); glVertex3f( 1,  1, -0.5f);
    glTexCoord2f(0, theight);      glVertex3f(-1,  1, -0.5f);
  }
  glEnd();
  
  nonGaussian_shader.deactivate();
 
  glBindTexture(texType, 0);
  
  //Prints out the updated prime values
  if(outputPrime)
    printPrime(odd,false);

}

void ParticleControl::setupPrime_and_AdvectShader(int numInRow,float life_time){
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
  GLint uNumInRow = mrt_shader.createUniform("numInRow");

  mrt_shader.activate();

  glUniform1fARB(ulifeTime, life_time);
  glUniform1iARB(unx, nx);
  glUniform1iARB(uny, ny);
  glUniform1iARB(unz, nz);
  glUniform1fARB(uNumInRow, numInRow);

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
    glTexCoord2f(0, 0);            glVertex3f(-1, -1, -0.5f);
    glTexCoord2f(twidth, 0);       glVertex3f( 1, -1, -0.5f);
    glTexCoord2f(twidth, theight); glVertex3f( 1,  1, -0.5f);
    glTexCoord2f(0, theight);      glVertex3f(-1,  1, -0.5f);
  }
  glEnd();
  
  mrt_shader.deactivate();
 
  glBindTexture(texType, 0);
  
  //Prints out the updated prime values
  if(outputPrime)
    printPrime(odd,false);

}

void ParticleControl::setupPrimeShader(int numInRow){  //Included argument -- Balli(04/12/07)
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
  glUniform1fARB(uNumInRow, numInRow);

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
void ParticleControl::setupAdvectShader(int numInRow, float life_time){

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
  //glUniform1fARB(uniform_timeStep, *time_step);
  glUniform1fARB(uNumInRow, numInRow);

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
    glTexCoord2f(0, 0);            glVertex3f(-1, -1, -0.5f);
    glTexCoord2f(twidth, 0);       glVertex3f( 1, -1, -0.5f);
    glTexCoord2f(twidth, theight); glVertex3f( 1,  1, -0.5f);
    glTexCoord2f(0, theight);      glVertex3f(-1,  1, -0.5f);
  }
  glEnd();
  
  pass1_shader.deactivate();
 
  glBindTexture(texType, 0);

}
void ParticleControl::setupCurrVel_shader(int numInRow){
  currVel_shader.addShader("Shaders/currVel_vp.glsl",GLSLObject::VERTEX_SHADER);
  currVel_shader.addShader("Shaders/currVel_fp.glsl",GLSLObject::FRAGMENT_SHADER);
  currVel_shader.createProgram();

  uniform_currentPrime = currVel_shader.createUniform("currPrime");
  uniform_windVelocity = currVel_shader.createUniform("windVel");
  uniform_partPos = currVel_shader.createUniform("position");
  uniform_prevPartPos = currVel_shader.createUniform("position_prev");
  
  GLint unir = currVel_shader.createUniform("numInRow");

  GLint unx = currVel_shader.createUniform("nx");
  GLint uny = currVel_shader.createUniform("ny");
  GLint unz = currVel_shader.createUniform("nz");

  currVel_shader.activate();

  glUniform1iARB(unx, nx);
  glUniform1iARB(uny, ny);
  glUniform1iARB(unz, nz);
  glUniform1fARB(unir,numInRow);

  currVel_shader.deactivate();

}
void ParticleControl::updateCurrVel(bool odd, GLuint prime0,GLuint prime1,GLuint windField,
				    GLuint positions0,GLuint positions1){

  glDrawBuffer(currVelBuffer);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glViewport(0, 0, twidth, theight);
    
  glEnable(texType);
  currVel_shader.activate();

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
    glTexCoord2f(0, 0);            glVertex3f(-1, -1, -0.5f);
    glTexCoord2f(twidth, 0);       glVertex3f( 1, -1, -0.5f);
    glTexCoord2f(twidth, theight); glVertex3f( 1,  1, -0.5f);
    glTexCoord2f(0, theight);      glVertex3f(-1,  1, -0.5f);
  }
  glEnd();

  currVel_shader.deactivate();

  glBindTexture(texType,0);
}
void ParticleControl::setupMeanVel_shader(int numInRow){
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

  glUniform1iARB(unx, nx);
  glUniform1iARB(uny, ny);
  glUniform1iARB(unz, nz);
  glUniform1fARB(unir,numInRow);

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
    glTexCoord2f(0, 0);            glVertex3f(-1, -1, -0.5f);
    glTexCoord2f(twidth, 0);       glVertex3f( 1, -1, -0.5f);
    glTexCoord2f(twidth, theight); glVertex3f( 1,  1, -0.5f);
    glTexCoord2f(0, theight);      glVertex3f(-1,  1, -0.5f);
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

void ParticleControl::initWindTex(GLuint windField, int* numInRow,
				  int dataSet){
  // Create wind velocity data texture
  wind_vel = new wind[nx*ny*nz];
  cellQuic = new cellType[nx*ny*nz];

  //Matrix tau11,tau22,tau33,and tau13
  tau = new Matrix[nx*ny*nz];

  //sigU,sigV,and sigW
  sig = new wind[nx*ny*nz];

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

  }

  /////////////////////////////////////////////////////////
  //Calculate width and height for wind texture
  //
  //This tries to minimize the width and height values
  //to try and fit the wind field into a 2D texture without
  //wasting too much space.  
  /////////////////////////////////////////////////////////
  int total = nx*ny*nz;
  width = (int)sqrt((float)total);
  
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
  height = width;

  ////////////////////////////////////////////////////////
  //Convert this to data array for a texture
  //
  //This will directly put the 3D data into an array
  //that is used to make the 2D texture.
  ///////////////////////////////////////////////////////
  (*numInRow) = (width - (width % nx))/nx;
  //std::cout << width << " " << *numInRow << std::endl;

  if(dataSet != 5){

    int qi, qj, qk;
    int p2idx = 0, texidx = 0;
    int row = 0;
  
    GLfloat *data = new GLfloat[ width * height * 4 ];

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
	  
	    data[texidx] = wind_vel[p2idx].u;
	    data[texidx+1] = wind_vel[p2idx].v;
	    data[texidx+2] = wind_vel[p2idx].w;	 
	    //if(wind_vel[p2idx].id == -1.0)
	    data[texidx+3] = (0.5*5.7)*(ustar*ustar*ustar)/(0.4*(qk+1));//This value is the '0.5*CoEps' value	
	      //else
	      //data[texidx+3] = wind_vel[p2idx].id;
	  }

    createTexture(windField, GL_RGBA32F_ARB, width, height, data);

    delete [] data;
  }

}
void ParticleControl::initLambda_and_TauTex(GLuint lambda, GLuint tau_dz, GLuint duvw_dz, int numInRow){
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
  float dz = 1.0;
  //Surface roughness
  float znaut = 0.01;
  
  for (qk=0; qk<nz; qk++) 
    for (qi=0; qi<ny; qi++)
      for (qj=0; qj<nx; qj++)
	{
	  p2idx = qk*ny*nx + qi*nx + qj;
	    
	  row = qk / (numInRow);
	  texidx = row * width * ny * 4 +
	  qi * width * 4 +
	  qk % (numInRow) * nx * 4 +
	  qj * 4;

	  rowBelow = (qk-1) / (numInRow);
	  texidxBelow = rowBelow * width * ny * 4 +
	    qi * width * 4 +
	    (qk-1) % (numInRow) * nx * 4 +
	    qj * 4;
	  
	  //The cells right above and below the current cell
	  idxBelow = (qk-1)*ny*nx + qi*nx + qj;
	  idxAbove = (qk+1)*ny*nx + qi*nx + qj;
	  //The cell two below the current cell
	  idx2Below = (qk-2)*ny*nx + qi*nx + qj;

	  //Calculate du/dz

	  //Cell just above the ground
	  //du/dz = Uat 1st cell/0.5dz(log(dz/znaut))
	  if(qk == 0){
	    du_dz = wind_vel[p2idx].u/(dz*log(dz/znaut));
	  }
	  //Cell at the top of the domain
	  else if(qk == (nz-1)){
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
	    du_dz = (wind_vel[idxAbove].u - wind_vel[idxBelow].u)/(2.0*dz);
	  }
	  	  
	  data3[texidx] = du_dz;    //du_dz
	  data3[texidx+1] = 0.0;    //eventually dv_dz
	  data3[texidx+2] = 0.0;    //eventually dw_dz
	  data3[texidx+3] = 0.0;
	  
	  ustar = 0.4*(qk+1)*du_dz;

	  sigU = 2.0*ustar;
	  sigV = 2.0*ustar;
	  sigW = 1.3*ustar;

	  sig[p2idx].u = sigU;   //sigU
	  sig[p2idx].v = sigV;   //sigV
	  sig[p2idx].w = sigW;   //sigW

	  float tau11=sigU*sigU;
	  float tau22=sigV*sigV;
	  float tau33=sigW*sigW;
	  float tau13=ustar*ustar;
	  float tauDetInv=1.0/((tau11*tau22*tau33)-(tau13*tau13*tau22));
	  
	  tau[p2idx].t11   = tau11;             //Tau11
	  tau[p2idx+1].t22 = tau22;             //Tau22
	  tau[p2idx+2].t33 = tau33;             //Tau33
	  tau[p2idx+3].t13 = tau13;             //Tau13

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
  
  for (qk=0; qk<nz; qk++) 
    for (qi=0; qi<ny; qi++)
      for (qj=0; qj<nx; qj++)
	{
	  p2idx = qk*ny*nx + qi*nx + qj;
	    
	  row = qk / (numInRow);
	  texidx = row * width * ny * 4 +
	  qi * width * 4 +
	  qk % (numInRow) * nx * 4 +
	  qj * 4;

	  rowBelow = (qk-1) / (numInRow);
	  texidxBelow = rowBelow * width * ny * 4 +
	    qi * width * 4 +
	    (qk-1) % (numInRow) * nx * 4 +
	    qj * 4;
	  
	  //The cells right above and below the current cell
	  idxBelow = (qk-1)*ny*nx + qi*nx + qj;
	  idxAbove = (qk+1)*ny*nx + qi*nx + qj;

	  //The cell two below the current cell
	  idx2Below = (qk-2)*ny*nx + qi*nx + qj;

	  //Calculate Tau/dz
	  
	  //Cell just above the ground
	  //dT/dz = -2(T/(0.5dz)log((0.5dz/znaut)))
	  if(qk == 0){
	    tau11_dz = -2.0*(tau[p2idx].t11/(dz*log(dz/znaut)));
	    tau22_dz = -2.0*(tau[p2idx].t22/(dz*log(dz/znaut)));
	    tau33_dz = -2.0*(tau[p2idx].t33/(dz*log(dz/znaut)));
	    tau13_dz = -2.0*(tau[p2idx].t13/(dz*log(dz/znaut)));
	  }
	  //Cell at the top of the domain
	  else if(qk == (nz-1)){
	    tau11_dz = data[texidxBelow];
	    tau22_dz = data[texidxBelow+1];
	    tau33_dz = data[texidxBelow+2];
	    tau13_dz = data[texidxBelow+3];
	  }
	  /*
	  //Cell at the top of the domain
	  //dT/dz = dT/dz at k-1th cell
	  //which is the same as dT/dz = Tk-Tk-2/2dz if domain is higher than 2		     
	  else if((qk == (nz-1)) && qk > 1){
	    tau11_dz = (tau[p2idx].t11 - tau[idx2Below].t11)/(2.0*dz);
	    tau22_dz = (tau[p2idx].t22 - tau[idx2Below].t22)/(2.0*dz);
	    tau33_dz = (tau[p2idx].t33 - tau[idx2Below].t33)/(2.0*dz);
	    tau13_dz = (tau[p2idx].t13 - tau[idx2Below].t13)/(2.0*dz);
	  }
	  //Cell at the top of the domain
	  //which is the same as cell just above ground if domain is only 2 high
	  else if((qk == (nz-1)) && qk == 1){
	    tau11_dz = -2.0*(tau[idxBelow].t11/(0.5*dz*log((0.5*dz)/znaut)));
	    tau22_dz = -2.0*(tau[idxBelow].t22/(0.5*dz*log((0.5*dz)/znaut)));
	    tau33_dz = -2.0*(tau[idxBelow].t33/(0.5*dz*log((0.5*dz)/znaut)));
	    tau13_dz = -2.0*(tau[idxBelow].t13/(0.5*dz*log((0.5*dz)/znaut)));
	    }*/


	  //All other cells
	  //dT/dz = (Tk+1 - Tk-1)/2dz
	  else{
	    tau11_dz = (tau[idxAbove].t11 - tau[idxBelow].t11)/(2.0*dz);
	    tau22_dz = (tau[idxAbove].t22 - tau[idxBelow].t22)/(2.0*dz);
	    tau33_dz = (tau[idxAbove].t33 - tau[idxBelow].t33)/(2.0*dz);
	    tau13_dz = (tau[idxAbove].t13 - tau[idxBelow].t13)/(2.0*dz);
	  }

	  data[texidx] = tau11_dz;
	  data[texidx+1] = tau22_dz;
	  data[texidx+2] = tau33_dz;
	  data[texidx+3] = tau13_dz;

	}

  createTexture(tau_dz, GL_RGBA32F_ARB, width,height, data);
  
  delete [] data;

}
void ParticleControl::initLambda_and_TauTex_fromQUICFILES(GLuint windField,GLuint lambda, GLuint tau_dz, GLuint duvw_dz, int numInRow){
  GLfloat *data = new GLfloat[ width * height * 4 ];
  GLfloat *dataWind = new GLfloat[ width * height * 4 ];
  GLfloat *dataTwo = new GLfloat[ width * height * 4 ];
  GLfloat *data3 = new GLfloat[width*height*4];
  GLfloat *data4 = new GLfloat[width*height*4];
  GLfloat *data5 = new GLfloat[width*height*4];

  int qi, qj, qk;
  int p2idx = 0, texidx = 0, texidxBelow = 0;
  int row = 0, rowBelow = 0;
  
  //cell below, cell above, cell two below
  int idxBelow, idxAbove, idx2Below,idxRight,idxLeft,idxFront,idxBehind;

  float du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,dw_dx,dw_dy,dw_dz;
  double S11,S22,S33,S21,S12,S31,S13,S23,S32;
  //Cell height
  float dz = 1.0,dx=1.0,dy=1.0;
  //Surface roughness
  float znaut = 0.1;

  //initializes the array cellQuic[]
  initCellType();
 
  //Reading turbulence data from the QUIC generated file
  std::ifstream turbulence;
  turbulence.open("Settings/QP_turbfield.dat");
  
  std::string header;
  
  turbulence>>header>>header>>header>>header>>header>>header>>header>>header>>header>>header;
  turbulence>>header>>header>>header>>header>>header>>header>>header>>header>>header>>header;
  turbulence>>header>>header>>header>>header>>header>>header>>header>>header>>header>>header;
  turbulence>>header>>header>>header>>header>>header>>header>>header>>header>>header>>header;
  turbulence>>header>>header>>header;
  //End reading the header, Now the values are read in the loop

  double indexVal,extraVal,elz,eps; // indexVal is for x,y and z locations , eps and elz are required for ustar calculations.

  for (qk=0; qk<nz; qk++) 
    for (qi=0; qi<ny; qi++)
      for (qj=0; qj<nx; qj++)
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

	  p2idx = qk*ny*nx + qi*nx + qj;
		
	  if(cellQuic[p2idx].c!=0){ //Calculate gradients ONLY if its not a building or ground cell    

	    // Calculating distance to the wall in all directions.

	    //to look for distance to a building in negative Z-direction (towards the ground)
	    /*int disNegZSurf=qk,qkk=qk;
	    while(qkk>=0){
	      int p2idx_k = qkk*ny*nx + qi*nx + qj;
	      if(cellQuic[p2idx_k].c==0){
		disNegZSurf=(qk-qkk)+1;
		break;
	      }
	      if(qkk==0){
		disNegZSurf=qk+1;
		break;
	      }
	      --qkk;
	    }
           
	    //For Y-direction.
	    int p2idx_i;
	    //to look for a building in negative Y-direction
	    int disNegYSurf=1000*ny , qii=qi;
	    while(qii>=0){
	      p2idx_i = qk*ny*nx + qii*nx + qj;
	      if(cellQuic[p2idx_i].c==0){
		disNegYSurf=qi-qii;
		break;
	      }
	      --qii;
	    }
	    //to look for a building in postive Y-direction
	    int disPosYSurf=1000*ny;
	    qii=qi;
	    while(qii<=ny){
	      p2idx_i = qk*ny*nx + qii*nx + qj;
	      if(cellQuic[p2idx_i].c==0){
		disPosYSurf=qii-qi;
		break;
	      }
	      ++qii;
	    }

	    //For X-direction.
	    int p2idx_j;
	    //to look for a building in negative X-direction
	    int disNegXSurf=1000*nx , qjj=qj;
	    while(qjj>=0){
	      p2idx_j = qk*ny*nx + qi*nx + qjj;
	      if(cellQuic[p2idx_j].c==0){
		disNegXSurf=qj-qjj;
		break;
	      }
	      --qjj;
	    }
	    //to look for a building in postive X-direction
	    int disPosXSurf=1000*nx;
	    qjj=qj;
	    while(qjj<=nx){
	      p2idx_j = qk*ny*nx + qi*nx + qjj;
	      if(cellQuic[p2idx_j].c==0){
		disPosXSurf=qjj-qj;
		break;
	      }
	      ++qjj;
	    }
	    //Calculating "minimum" distance to the wall.
	    int disArray[]={disNegXSurf,disPosXSurf,disNegYSurf,disPosYSurf,disNegZSurf};

	    int minDistance = disNegXSurf;  //let the lowest value be the first in the array

	    for(int i=1; i<(int)(sizeof(disArray)/sizeof*(disArray));i++){//loop for getting the actual minimum distance
	      if(disArray[i] <= minDistance){
		minDistance = disArray[i];
	      }
	    }*/

	    row = qk / (numInRow);
	    texidx = row * width * ny * 4 +
	      qi * width * 4 +
	      qk % (numInRow) * nx * 4 +
	      qj * 4;

	    rowBelow = (qk-1) / (numInRow);
	    texidxBelow = rowBelow * width * ny * 4 +
	      qi * width * 4 +
	      (qk-1) % (numInRow) * nx * 4 +
	      qj * 4;
	         
	    //The cells right above and below the current cell
	    idxBelow = (qk-1)*ny*nx + qi*nx + qj;
	    idxAbove = (qk+1)*ny*nx + qi*nx + qj;         

	    //The cells at front and behind of the current cell   
	    idxBehind = qk*ny*nx + qi*nx + (qj-1);
	    idxFront  = qk*ny*nx + qi*nx + (qj+1);
	    //The cell two below the current cell
	    idx2Below = (qk-2)*ny*nx + qi*nx + qj;
            
	    //Calculate du/dz

	    //Cell just above the ground
	    //du/dz = Uat 1st cell/0.5dz(log(dz/znaut))
	    /*if(qk == 0 || (cellQuic[idxBelow].c==0 && cellQuic[p2idx].c==1)){
	      du_dz = wind_vel[p2idx].u/(dz*log(dz/znaut));
	      dv_dz = wind_vel[p2idx].v/(dz*log(dz/znaut));
	      dw_dz = wind_vel[p2idx].w/(dz*log(dz/znaut));
	      ustar = 0.4*(minDistance+1)*du_dz;
	    }
	    //Cell at the top of the domain
	    else if(qk == (nz-1)){
	      du_dz = data3[texidxBelow];  // du_dz at k-1th cell
	      dv_dz = data3[texidxBelow+1];
	      dw_dz = data3[texidxBelow+2];
	      ustar = 0.4*(minDistance+1)*du_dz;
	    }

	    //All other cells. i.e not boundary cells
	    //du/dz = (Uk+1 - Uk-1)/2dz
	    else{
	      du_dz = (wind_vel[idxAbove].u - wind_vel[idxBelow].u)/(2.0*dz);
	      dv_dz = (wind_vel[idxAbove].v - wind_vel[idxBelow].v)/(2.0*dz);
	      dw_dz = (wind_vel[idxAbove].w - wind_vel[idxBelow].w)/(2.0*dz);
	      ustar = 0.4*(minDistance+1)*du_dz;
	    }*/
	       				  	  
	    data3[texidx] = (sigU/2.5)/elz;  //du_dz;      //du_dz
	    data3[texidx+1] = 0.0;    //dv_dz
	    data3[texidx+2] = 0.0;    //dw_dz
	    data3[texidx+3] = 0.0;
		  
	    //For gradient in X-direction
		   
	    //Cell at the top of the domain
	    /*if(qj == (nx-1) || qj==0){
	      du_dx = 0;  // du_dz at k-1th cell
	      dv_dx = 0;
	      dw_dx = 0;
	    }
	    else if((cellQuic[idxBehind].c==0 && cellQuic[p2idx].c==1) || (cellQuic[idxFront].c==0 && cellQuic[p2idx].c==1)){
	      du_dx = wind_vel[p2idx].u/(dx*log(dx/znaut));
	      dv_dx = wind_vel[p2idx].v/(dx*log(dx/znaut));
	      dw_dx = wind_vel[p2idx].w/(dx*log(dx/znaut));
	    }
	    else{
	      du_dx = (wind_vel[idxFront].u - wind_vel[idxBehind].u)/(2.0*dx);
	      dv_dx = (wind_vel[idxFront].v - wind_vel[idxBehind].v)/(2.0*dx);
	      dw_dx = (wind_vel[idxFront].w - wind_vel[idxBehind].w)/(2.0*dx);
	    }

	    data4[texidx]   = du_dx;    //du_dx
	    data4[texidx+1] = dv_dx;    //eventually dv_dx
	    data4[texidx+2] = dw_dx;    //eventually dw_dx
	    data4[texidx+3] = 0.0;	      

	    //For gradient in Y-direction
	    //The cells at right and left of the current cell 
	    idxRight = qk*ny*nx + (qi-1)*nx + qj;
	    idxLeft  = qk*ny*nx + (qi+1)*nx + qj;
		   
	    //Cell at the top of the domain
	    if(qi == (ny-1) || qi==0){
	      du_dy = 0;  // du_dz at k-1th cell
	      dv_dy = 0;
	      dw_dy = 0;
	    }
	    else if((cellQuic[idxLeft].c==0 && cellQuic[p2idx].c==1) || (cellQuic[idxRight].c==0 && cellQuic[p2idx].c==1)){
	      du_dy = wind_vel[p2idx].u/(dy*log(dy/znaut));
	      dv_dy = wind_vel[p2idx].v/(dy*log(dy/znaut));
	      dw_dy = wind_vel[p2idx].w/(dy*log(dy/znaut));
	    }
	    else{
	      du_dy = (wind_vel[idxLeft].u - wind_vel[idxRight].u)/(2.0*dy);
	      dv_dy = (wind_vel[idxLeft].v - wind_vel[idxRight].v)/(2.0*dy);
	      dw_dy = (wind_vel[idxLeft].w - wind_vel[idxRight].w)/(2.0*dy);
	    }

	    data5[texidx]   = du_dy;    //du_dy
	    data5[texidx+1] = dv_dy;    //eventually dv_dy
	    data5[texidx+2] = dw_dy;    //eventually dw_dy
	    data5[texidx+3] = 0.0;
		

	    double velGrad[]={du_dx,dv_dx,dw_dx,du_dy,dv_dy,dw_dy,du_dz,dv_dz,dw_dz};
	    double maxVelGradAbs=fabs(du_dx);
	    double maxVelGrad;

	    for(int i=1; i<(int)(sizeof(velGrad)/sizeof*(velGrad));i++){//loop for getting the actual max velocity gradient
	      if(fabs(velGrad[i]) >= maxVelGradAbs){
		maxVelGradAbs = fabs(velGrad[i]);
		maxVelGrad = velGrad[i];
	      }
	    }
	    
	    double VertGradFactor=pow( (1-(minDistance/20)) ,3.0/2.0); 
	    ustar=0.4*minDistance*maxVelGrad*VertGradFactor;
        */
	   
        ustar=sigU/2.5;
		data3[texidx+3] = ustar;//temporarily used for the ustar.


	    /*sigU = 2.0*ustar;
	    sigV = 2.0*ustar;
	    sigW = 1.3*ustar;*/

	    sig[p2idx].u = sigU;   //sigU
	    sig[p2idx].v = sigV;   //sigV
	    sig[p2idx].w = sigW;   //sigW

	    float tau11=sigU*sigU;
	    float tau22=sigV*sigV;
	    float tau33=sigW*sigW;
	    float tau13=ustar*ustar;
	    float tauDetInv=1.0/((tau11*tau22*tau33)-(tau13*tau13*tau22));

	  
	    tau[p2idx].t11   = tau11;             //Tau11
	    tau[p2idx+1].t22 = tau22;             //Tau22
	    tau[p2idx+2].t33 = tau33;             //Tau33
	    tau[p2idx+3].t13 = tau13;             //Tau13


	    dataTwo[texidx]   =  1.0/(tau11-tau13*tau13/tau33);// tauDetInv*(tau22*tau33);            //Lam11
	    dataTwo[texidx+1] =  1.0/tau22;// tauDetInv*(tau11*tau33-tau13*tau13);//Lam22
	    dataTwo[texidx+2] =  1.0/(tau33-tau13*tau13/tau11);//tauDetInv*(tau11*tau22);	          //Lam33
	    dataTwo[texidx+3] =  -tau13/(tau11*tau33-tau13*tau13);//tauDetInv*(-tau13*tau22);           //Lam13
        
	    dataWind[texidx] = wind_vel[p2idx].u;
	    dataWind[texidx+1] = wind_vel[p2idx].v;
	    dataWind[texidx+2] = wind_vel[p2idx].w;	  
	    dataWind[texidx+3] = (0.5*5.7)*eps;//(0.5*5.7)*(ustar*ustar*ustar)/(0.4*(minDistance));//This value is the '0.5*CoEps' value	
	  }

	}
  
  createTexture(lambda, GL_RGBA32F_ARB, width,height, dataTwo);
  createTexture(duvw_dz, GL_RGBA32F_ARB, width,height, data3);
  createTexture(windField, GL_RGBA32F_ARB, width, height, dataWind);
  
  delete [] dataWind;
  delete [] data3;
  delete [] dataTwo;

  //Create the Tau/dz texture
  float tau11_dz, tau22_dz, tau33_dz, tau13_dz;
  
  for (qk=0; qk<nz; qk++) 
    for (qi=0; qi<ny; qi++)
      for (qj=0; qj<nx; qj++)
	{
	  p2idx = qk*ny*nx + qi*nx + qj;
	    
	  row = qk / (numInRow);
	  texidx = row * width * ny * 4 +
	    qi * width * 4 +
	    qk % (numInRow) * nx * 4 +
	    qj * 4;

	  rowBelow = (qk-1) / (numInRow);
	  texidxBelow = rowBelow * width * ny * 4 +
	    qi * width * 4 +
	    (qk-1) % (numInRow) * nx * 4 +
	    qj * 4;
	  
	  //The cells right above and below the current cell
	  idxBelow = (qk-1)*ny*nx + qi*nx + qj;
	  idxAbove = (qk+1)*ny*nx + qi*nx + qj;

	  //The cell two below the current cell
	  idx2Below = (qk-2)*ny*nx + qi*nx + qj;

	  //Calculate Tau/dz
	  
	  //Cell just above the ground
	  //dT/dz = -2(T/(0.5dz)log((0.5dz/znaut)))
	  if(qk == 0 || (cellQuic[idxBelow].c==0 && cellQuic[p2idx].c==1)){
	    tau11_dz = -2.0*(tau[p2idx].t11/(dz*log(dz/znaut)));
	    tau22_dz = -2.0*(tau[p2idx].t22/(dz*log(dz/znaut)));
	    tau33_dz = -2.0*(tau[p2idx].t33/(dz*log(dz/znaut)));
	    tau13_dz = -2.0*(tau[p2idx].t13/(dz*log(dz/znaut)));
	  }
	  //Cell at the top of the domain
	  else if(qk == (nz-1)){
	    tau11_dz = data[texidxBelow];
	    tau22_dz = data[texidxBelow+1];
	    tau33_dz = data[texidxBelow+2];
	    tau13_dz = data[texidxBelow+3];
	  }
	  //All other cells
	  //dT/dz = (Tk+1 - Tk-1)/2dz
	  else{
	    tau11_dz = (tau[idxAbove].t11 - tau[idxBelow].t11)/(2.0*dz);
	    tau22_dz = (tau[idxAbove].t22 - tau[idxBelow].t22)/(2.0*dz);
	    tau33_dz = (tau[idxAbove].t33 - tau[idxBelow].t33)/(2.0*dz);
	    tau13_dz = (tau[idxAbove].t13 - tau[idxBelow].t13)/(2.0*dz);
	  }

	  data[texidx] = tau11_dz;
	  data[texidx+1] = tau22_dz;
	  data[texidx+2] = tau33_dz;
	  data[texidx+3] = tau13_dz;

	}

  createTexture(tau_dz, GL_RGBA32F_ARB, width,height, data);
  
  delete [] data;

}
/*void ParticleControl::initLambda_and_TauTex_fromQUICFILES(GLuint lambda, GLuint tau_dz, GLuint duvw_dz, int numInRow){
  GLfloat *data = new GLfloat[ width * height * 4 ];
  GLfloat *dataTwo = new GLfloat[ width * height * 4 ];
  GLfloat *data3 = new GLfloat[width*height*4];
  GLfloat *data4 = new GLfloat[width*height*4];
  GLfloat *data5 = new GLfloat[width*height*4];

  int qi, qj, qk;
  int p2idx = 0, texidx = 0, texidxBelow = 0;
  int row = 0, rowBelow = 0;
  
  //cell below, cell above, cell two below
  int idxBelow, idxAbove, idx2Below,idxRight,idxLeft,idxFront,idxBehind;

  float du_dx,du_dy,du_dz,dv_dx,dv_dy,dv_dz,dw_dx,dw_dy,dw_dz;
  double S11,S22,S33,S21,S12,S31,S13,S23,S32;
  //Cell height
  float dz = 1.0,dx=1.0,dy=1.0;
  //Surface roughness
  float znaut = 0.01;

  //initializes the array cellQuic[]
  initCellType();
  
  for (qk=0; qk<nz; qk++) 
    for (qi=0; qi<ny; qi++)
      for (qj=0; qj<nx; qj++)
	{	
	  p2idx = qk*ny*nx + qi*nx + qj;
		
	  if(cellQuic[p2idx].c!=0){ //Calculate gradients ONLY if its not a building or ground cell    

	    // Calculating distance to the wall in all directions.

	    //to look for distance to a building in negative Z-direction (towards the ground)
	    int disNegZSurf=qk,qkk=qk;
	    while(qkk>=0){
	      int p2idx_k = qkk*ny*nx + qi*nx + qj;
	      if(cellQuic[p2idx_k].c==0){
		disNegZSurf=qk-qkk;
		break;
	      }
	      --qkk;
	    }
           
	    //For Y-direction.
	    int p2idx_i;
	    //to look for a building in negative Y-direction
	    int disNegYSurf=1000*ny , qii=qi;
	    while(qii>=0){
	      p2idx_i = qk*ny*nx + qii*nx + qj;
	      if(cellQuic[p2idx_i].c==0){
		disNegYSurf=qi-qii;
		break;
	      }
	      --qii;
	    }
	    //to look for a building in postive Y-direction
	    int disPosYSurf=1000*ny;
	    qii=qi;
	    while(qii<=ny){
	      p2idx_i = qk*ny*nx + qii*nx + qj;
	      if(cellQuic[p2idx_i].c==0){
		disPosYSurf=qii-qi;
		break;
	      }
	      ++qii;
	    }

	    //For X-direction.
	    int p2idx_j;
	    //to look for a building in negative X-direction
	    int disNegXSurf=1000*nx , qjj=qj;
	    while(qjj>=0){
	      p2idx_j = qk*ny*nx + qi*nx + qjj;
	      if(cellQuic[p2idx_j].c==0){
		disNegXSurf=qj-qjj;
		break;
	      }
	      --qjj;
	    }
	    //to look for a building in postive X-direction
	    int disPosXSurf=1000*nx;
	    qjj=qj;
	    while(qjj<=nx){
	      p2idx_j = qk*ny*nx + qi*nx + qjj;
	      if(cellQuic[p2idx_j].c==0){
		disPosXSurf=qjj-qj;
		break;
	      }
	      ++qjj;
	    }
	    //Calculating "minimum" distance to the wall.
	    int disArray[]={disNegXSurf,disPosXSurf,disNegYSurf,disPosYSurf,disNegZSurf};

	    int minDistance = disNegXSurf;  //let the lowest value be the first in the array

	    for(int i=1; i<(int)(sizeof(disArray)/sizeof*(disArray));i++){//loop for getting the actual minimum distance
	      if(disArray[i] <= minDistance){
		minDistance = disArray[i];
	      }
	    }

	    row = qk / (numInRow);
	    texidx = row * width * ny * 4 +
	      qi * width * 4 +
	      qk % (numInRow) * nx * 4 +
	      qj * 4;

	    rowBelow = (qk-1) / (numInRow);
	    texidxBelow = rowBelow * width * ny * 4 +
	      qi * width * 4 +
	      (qk-1) % (numInRow) * nx * 4 +
	      qj * 4;
	         
	    //The cells right above and below the current cell
	    idxBelow = (qk-1)*ny*nx + qi*nx + qj;
	    idxAbove = (qk+1)*ny*nx + qi*nx + qj;         

	    //The cells at front and behind of the current cell   
	    idxBehind = qk*ny*nx + qi*nx + (qj-1);
	    idxFront  = qk*ny*nx + qi*nx + (qj+1);
	    //The cell two below the current cell
	    idx2Below = (qk-2)*ny*nx + qi*nx + qj;
            
	    //Calculate du/dz

	    //Cell just above the ground
	    //du/dz = Uat 1st cell/0.5dz(log(dz/znaut))
	    if(qk == 0 || (cellQuic[idxBelow].c==0 && cellQuic[p2idx].c==1)){
	      du_dz = wind_vel[p2idx].u/(dz*log(dz/znaut));
	      dv_dz = wind_vel[p2idx].v/(dz*log(dz/znaut));
	      dw_dz = wind_vel[p2idx].w/(dz*log(dz/znaut));
	      ustar = 0.4*(minDistance+1)*du_dz;
	    }
	    //Cell at the top of the domain
	    else if(qk == (nz-1)){
	      du_dz = data3[texidxBelow];  // du_dz at k-1th cell
	      dv_dz = data3[texidxBelow+1];
	      dw_dz = data3[texidxBelow+2];
	      ustar = 0.4*(minDistance+1)*du_dz;
	    }

	    //All other cells. i.e not boundary cells
	    //du/dz = (Uk+1 - Uk-1)/2dz
	    else{
	      du_dz = (wind_vel[idxAbove].u - wind_vel[idxBelow].u)/(2.0*dz);
	      dv_dz = (wind_vel[idxAbove].v - wind_vel[idxBelow].v)/(2.0*dz);
	      dw_dz = (wind_vel[idxAbove].w - wind_vel[idxBelow].w)/(2.0*dz);
	      ustar = 0.4*(minDistance+1)*du_dz;
	    }
	       				  	  
	    data3[texidx] = du_dz;      //du_dz
	    data3[texidx+1] = dv_dz;    //dv_dz
	    data3[texidx+2] = dw_dz;    //dw_dz
	    data3[texidx+3] = 0.0;
		  
	    //For gradient in X-direction
		   
	    //Cell at the top of the domain
	    if(qj == (nx-1) || qj==0){
	      du_dx = 0;  // du_dz at k-1th cell
	      dv_dx = 0;
	      dw_dx = 0;
	    }
	    else if((cellQuic[idxBehind].c==0 && cellQuic[p2idx].c==1) || (cellQuic[idxFront].c==0 && cellQuic[p2idx].c==1)){
	      du_dx = wind_vel[p2idx].u/(dx*log(dx/znaut));
	      dv_dx = wind_vel[p2idx].v/(dx*log(dx/znaut));
	      dw_dx = wind_vel[p2idx].w/(dx*log(dx/znaut));
	    }
	    else{
	      du_dx = (wind_vel[idxFront].u - wind_vel[idxBehind].u)/(2.0*dx);
	      dv_dx = (wind_vel[idxFront].v - wind_vel[idxBehind].v)/(2.0*dx);
	      dw_dx = (wind_vel[idxFront].w - wind_vel[idxBehind].w)/(2.0*dx);
	    }

	    data4[texidx]   = du_dx;    //du_dx
	    data4[texidx+1] = dv_dx;    //eventually dv_dx
	    data4[texidx+2] = dw_dx;    //eventually dw_dx
	    data4[texidx+3] = 0.0;	      

	    //For gradient in Y-direction
	    //The cells at right and left of the current cell 
	    idxRight = qk*ny*nx + (qi-1)*nx + qj;
	    idxLeft  = qk*ny*nx + (qi+1)*nx + qj;
		   
	    //Cell at the top of the domain
	    if(qi == (ny-1) || qi==0){
	      du_dy = 0;  // du_dz at k-1th cell
	      dv_dy = 0;
	      dw_dy = 0;
	    }
	    else if((cellQuic[idxLeft].c==0 && cellQuic[p2idx].c==1) || (cellQuic[idxRight].c==0 && cellQuic[p2idx].c==1)){
	      du_dy = wind_vel[p2idx].u/(dy*log(dy/znaut));
	      dv_dy = wind_vel[p2idx].v/(dy*log(dy/znaut));
	      dw_dy = wind_vel[p2idx].w/(dy*log(dy/znaut));
	    }
	    else{
	      du_dy = (wind_vel[idxLeft].u - wind_vel[idxRight].u)/(2.0*dy);
	      dv_dy = (wind_vel[idxLeft].v - wind_vel[idxRight].v)/(2.0*dy);
	      dw_dy = (wind_vel[idxLeft].w - wind_vel[idxRight].w)/(2.0*dy);
	    }

	    data5[texidx]   = du_dy;    //du_dy
	    data5[texidx+1] = dv_dy;    //eventually dv_dy
	    data5[texidx+2] = dw_dy;    //eventually dw_dy
	    data5[texidx+3] = 0.0;


	    S11=du_dx;
	    S22=dv_dy;
	    S33=dw_dz;

	    S12=0.5 * ( du_dy + dv_dz );
	    S13=0.5 * ( du_dz + dw_dx );
	    S32=0.5 * ( dw_dy + dv_dz );
		   
	    S21=S12;
	    S31=S13;
	    S23=S32;
              
	    sigU = 2.0*ustar;
	    sigV = 2.0*ustar;
	    sigW = 1.3*ustar;

	    sig[p2idx].u = sigU;   //sigU
	    sig[p2idx].v = sigV;   //sigV
	    sig[p2idx].w = sigW;   //sigW

	    float tau11=pow(0.4*(minDistance+1)*S11,2);//sigU*sigU;
	    float tau22=pow(0.4*(minDistance+1)*S22,2);//sigV*sigV;
	    float tau33=pow(0.4*(minDistance+1)*S33,2);//sigW*sigW;
	    float tau13=pow(0.4*(minDistance+1)*S13,2);//ustar*ustar;
	    float tauDetInv=1.0/((tau11*tau22*tau33)-(tau13*tau13*tau22));

	  
	    tau[p2idx].t11   = tau11;             //Tau11
	    tau[p2idx+1].t22 = tau22;             //Tau22
	    tau[p2idx+2].t33 = tau33;             //Tau33
	    tau[p2idx+3].t13 = tau13;             //Tau13

	    dataTwo[texidx]   = tauDetInv*(tau22*tau33);            //Lam11
	    dataTwo[texidx+1] = tauDetInv*(tau11*tau33-tau13*tau13);//Lam22
	    dataTwo[texidx+2] = tauDetInv*(tau11*tau22);	          //Lam33
	    dataTwo[texidx+3] = tauDetInv*(-tau13*tau22);           //Lam13
	  }
	}
  createTexture(lambda, GL_RGBA32F_ARB, width,height, dataTwo);
  createTexture(duvw_dz, GL_RGBA32F_ARB, width,height, data3);
  
  delete [] data3;
  delete [] dataTwo;

  //Create the Tau/dz texture
  float tau11_dz, tau22_dz, tau33_dz, tau13_dz;
  
  for (qk=0; qk<nz; qk++) 
    for (qi=0; qi<ny; qi++)
      for (qj=0; qj<nx; qj++)
	{
	  p2idx = qk*ny*nx + qi*nx + qj;
	    
	  row = qk / (numInRow);
	  texidx = row * width * ny * 4 +
	    qi * width * 4 +
	    qk % (numInRow) * nx * 4 +
	    qj * 4;

	  rowBelow = (qk-1) / (numInRow);
	  texidxBelow = rowBelow * width * ny * 4 +
	    qi * width * 4 +
	    (qk-1) % (numInRow) * nx * 4 +
	    qj * 4;
	  
	  //The cells right above and below the current cell
	  idxBelow = (qk-1)*ny*nx + qi*nx + qj;
	  idxAbove = (qk+1)*ny*nx + qi*nx + qj;

	  //The cell two below the current cell
	  idx2Below = (qk-2)*ny*nx + qi*nx + qj;

	  //Calculate Tau/dz
	  
	  //Cell just above the ground
	  //dT/dz = -2(T/(0.5dz)log((0.5dz/znaut)))
	  if(qk == 0 || (cellQuic[idxBelow].c==0 && cellQuic[p2idx].c==1)){
	    tau11_dz = -2.0*(tau[p2idx].t11/(dz*log(dz/znaut)));
	    tau22_dz = -2.0*(tau[p2idx].t22/(dz*log(dz/znaut)));
	    tau33_dz = -2.0*(tau[p2idx].t33/(dz*log(dz/znaut)));
	    tau13_dz = -2.0*(tau[p2idx].t13/(dz*log(dz/znaut)));
	  }
	  //Cell at the top of the domain
	  else if(qk == (nz-1)){
	    tau11_dz = data[texidxBelow];
	    tau22_dz = data[texidxBelow+1];
	    tau33_dz = data[texidxBelow+2];
	    tau13_dz = data[texidxBelow+3];
	  }
	  //All other cells
	  //dT/dz = (Tk+1 - Tk-1)/2dz
	  else{
	    tau11_dz = (tau[idxAbove].t11 - tau[idxBelow].t11)/(2.0*dz);
	    tau22_dz = (tau[idxAbove].t22 - tau[idxBelow].t22)/(2.0*dz);
	    tau33_dz = (tau[idxAbove].t33 - tau[idxBelow].t33)/(2.0*dz);
	    tau13_dz = (tau[idxAbove].t13 - tau[idxBelow].t13)/(2.0*dz);
	  }

	  data[texidx] = tau11_dz;
	  data[texidx+1] = tau22_dz;
	  data[texidx+2] = tau33_dz;
	  data[texidx+3] = tau13_dz;

	}

  createTexture(tau_dz, GL_RGBA32F_ARB, width,height, data);
  
  delete [] data;

  }*/
void ParticleControl::initLambdaTex(GLuint lambda, int numInRow){

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

   for (qk=0; qk<nz; qk++) 
    for (qi=0; qi<ny; qi++)
      for (qj=0; qj<nx; qj++)
	{
	  p2idx = qk*ny*nx + qi*nx + qj;
	    
	  row = qk / (numInRow);
	  texidx = row * width * ny * 4 +
	  qi * width * 4 +
	  qk % (numInRow) * nx * 4 +
	  qj * 4;
	  
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
  for(int k = 0; k < nz; k++){   
    for(int i = 0; i < ny; i++){
      for(int j = 0; j < nx; j++){
	int p2idx = k*nx*ny + i*nx + j;
  	if(i%2 == 0){
	  wind_vel[p2idx].u = 0;
	  wind_vel[p2idx].v = 1.0;
	  wind_vel[p2idx].w = 0;
	  wind_vel[p2idx].id = -1.0;
	 }
	else {
	  wind_vel[p2idx].u = 0;
	  wind_vel[p2idx].v = 0;
	  wind_vel[p2idx].w = 0;
	  wind_vel[p2idx].id = -1.0;
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
	wind_vel[p2idx].u = 0.0;//randVal();
	wind_vel[p2idx].v = 1.0;//randVal();
	wind_vel[p2idx].w = 0.0;//randVal();
	wind_vel[p2idx].id = -1.0;
      }
    }
  }
}

void ParticleControl::uniformUWindField(){
  for(int k = 0; k < nz; k++){   
    for(int i = 0; i < ny; i++){
      for(int j = 0; j < nx; j++){
	int p2idx = k*nx*ny + i*nx + j;
	wind_vel[p2idx].u = 1.0;
	wind_vel[p2idx].v = 0.0;
	wind_vel[p2idx].w = 0.0;
	wind_vel[p2idx].id = -1.0;
      }
    }
  }
}
void ParticleControl::variedUWindField(){
   for(int k = 0; k < nz; k++){   
    for(int i = 0; i < ny; i++){
      for(int j = 0; j < nx; j++){
	int p2idx = k*nx*ny + i*nx + j;
	wind_vel[p2idx].u = 7.52*pow(((k+1)/20.0),0.15);
	wind_vel[p2idx].v = 0.0;
	wind_vel[p2idx].w = 0.0;
	//wind_vel[p2idx].id = -1.0;
      }
    }
   }
   /*if(xfo != NULL)
     addBuildingsInWindField();
   else
   std::cout << "NO BUILDINGS ADDED TO WIND FIELD" << std::endl;*/
}
void ParticleControl::addBuildingsInWindField(GLuint cellType,int numInRow){

  GLfloat* data = new GLfloat[width*height*4];
  wind* cell_type = new wind[nx*ny*nz];
  
  for(int k = 0; k < nz; k++){   
    for(int i = 0; i < ny; i++){
      for(int j = 0; j < nx; j++){
	int p2idx = k*nx*ny + i*nx + j;
	cell_type[p2idx].u = 1.0;
	cell_type[p2idx].v = 1.0;
	cell_type[p2idx].w = 1.0;
	cell_type[p2idx].id = 0.0;
      }
    }
   }
  

  for(int n=0; n < numBuild; n++){

     for(int k=(int)zfo[n]; k < (int)(zfo[n]+ht[n]); k++){
       for(int i=(int)(yfo[n]-wti[n]/2.0); i < (int)(yfo[n]+wti[n]/2.0); i++){
	 for(int j=(int)xfo[n]; j < (int)(xfo[n]+lti[n]); j++){
	   int p2idx = k*nx*ny + i*nx + j;
	   
	   //wind_vel[p2idx].u = 0.0;
	   //wind_vel[p2idx].v = 0.0;
	   //wind_vel[p2idx].w = 0.0;
	   //wind_vel[p2idx].id = n;
	   cell_type[p2idx].u = 0.0;
	   cell_type[p2idx].v = 0.0;
	   cell_type[p2idx].w = 0.0;
	   cell_type[p2idx].id = n;

	 }
       }
     }
   }

   int qk,qi,qj;
   int row,texidx,p2idx;

   for (qk=0; qk<nz; qk++) 
      for (qi=0; qi<ny; qi++)
	for (qj=0; qj<nx; qj++)
	  {
	    p2idx = qk*ny*nx + qi*nx + qj;
	    
	    row = qk / (numInRow);
	    texidx = row * width * ny * 4 +
	      qi * width * 4 +
	      qk % (numInRow) * nx * 4 +
	      qj * 4;
	  
	    data[texidx] = cell_type[p2idx].u;
	    data[texidx+1] = cell_type[p2idx].v;
	    data[texidx+2] = cell_type[p2idx].w;	    
	    data[texidx+3] = cell_type[p2idx].id;
	    
	  }

    createTexture(cellType, GL_RGBA32F_ARB, width, height, data);

    delete [] data;
    delete [] cell_type;

}
void ParticleControl::QUICWindField(){
  std::ifstream QUICWindField;//,QUICCellType;
	
  QUICWindField.open("Settings/QU_velocity.dat"); //opening the wind file  to read
  if(!QUICWindField){
    std::cerr<<"Unable to open QUIC Windfield file : QU_velocity.dat ";
    exit(1);
  }

  std::string header;  //I am just using a very crude method to read the header of the wind file
  QUICWindField>>header>>header>>header>>header>>header>>header>>header>>header>>header>>header>>header>>header>>header>>header>>header>>header;
  QUICWindField>>header>>header>>header>>header>>header;
    
  double groundVal; // ignoring the ground values, so that k=0 have the first cell having non-zero velocity 
  //Balli had ++k?

  for(int k=0;k<(6*nx*ny);++k){ // there are 6 columns in the wind file 
    QUICWindField>>groundVal;
  }

  double quicIndex;

  for(int k = 0; k < nz; k++){   
    for(int i = 0; i < ny; i++){
      for(int j = 0; j < nx; j++){
	int p2idx = k*nx*ny + i*nx + j;
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
  
}
void ParticleControl::initCellType(){
  std::ifstream QUICCellType;
  QUICCellType.open("Settings/QU_celltype.dat");//opening the Celltype file  to read
  if(!QUICCellType){
    std::cerr<<"Unable to open QUIC Celltype file : QU_celltype.dat ";
    exit(1);
  }
  double groundVal;
  // ignoring the ground values for the cellype also
  //Balli had ++k ?

  for(int k=0;k<(4*nx*ny);++k){// there are 4 columns in the wind file
    QUICCellType>>groundVal;
  }
  double quicIndex;
  for(int k = 0; k < nz; k++){   
    for(int i = 0; i < ny; i++){
      for(int j = 0; j < nx; j++){
	int p2idx = k*nx*ny + i*nx + j;
	
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
