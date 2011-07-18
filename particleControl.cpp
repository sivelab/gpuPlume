#include <cassert>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <exception>
#include <algorithm>

#include <math.h>
#include <float.h>
#include "particleControl.h"
#include "Random.h"
#include "glErrorUtil.h"
#include "util.h"


ParticleControl::ParticleControl()
{
}
ParticleControl::~ParticleControl()
{
  delete windFieldLookupFBO;
}
ParticleControl::ParticleControl(GLenum type,int width,int height,
				 int x, int y, int z,
				 float c_dx, float c_dy, float c_dz, Util* util_ptr){

  m_util_ptr = util_ptr;

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

#if 0
  // Bad 
  nzdz = (int)(nz*(1.0/c_dz));
  nydy = (int)(ny*(1.0/c_dy));
  nxdx = (int)(nx*(1.0/c_dx));
#endif

  nzdz = (int)(nz);
  nydy = (int)(ny);
  nxdx = (int)(nx);

  cell_dx = c_dx;
  cell_dy = c_dy;
  cell_dz = c_dz;
}
void ParticleControl::setBuildingParameters(int nB,int* ns,float* x,float* y,float* z,
					   float* h,float* w,float* l,float* g){
  numBuild = nB;
  numSides = ns;
  xfo = x;
  yfo = y;
  zfo = z;
  ht = h;
  wti = w;
  lti = l; 
  gamma =g;
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
 uniform_dxyz_wall = multipleBuildings_shader.createUniform("dxy_wall");//Balli: new texture for storing distance to a surface in all directions  
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
					      GLuint duvw_dz, GLuint dxyz_wall,float time_step, GLuint buildings,
					      GLuint cellType, GLuint advect_terms)
{
  CheckErrorsGL("PC: multipleBuildingsAdvect() - entered");

  //Prints out the previous prime values
  if(outputPrime)
    printPrime(odd, true);


  //If coloring particles using advect_terms, advect_terms will not be
  //NULL here.

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

  CheckErrorsGL("PC: multipleBuildingsAdvect() - activated advection shader");

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

 //Bind the dxyz_wall texture to TEXTURE UNIT 9 (was 18 but that might not exist on some hardware!)
  glActiveTexture(GL_TEXTURE9);
  glBindTexture(texType, dxyz_wall);
  glUniform1iARB(uniform_dxyz_wall, 9);
    

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

  CheckErrorsGL("PC: multipleBuildingsAdvect() - shader advection completion");
 
  glBindTexture(texType, 0);
  
  //Prints out the updated prime values
  if(outputPrime)
    printPrime(odd,false);

}

void ParticleControl::setupWindFieldLookupShader() {
  windFieldLookup_shader.addShader("Shaders/windFieldLookup_vp.glsl", GLSLObject::VERTEX_SHADER);
  windFieldLookup_shader.addShader("Shaders/windFieldLookup_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  windFieldLookup_shader.createProgram();
  
  // Get location of sampler uniforms
  uniform_wfl_wind = windFieldLookup_shader.createUniform("wind_texunit");
  uniform_wfl_pos = windFieldLookup_shader.createUniform("pos");
  
  GLint unxdx = windFieldLookup_shader.createUniform("nxdx");
  GLint unydy = windFieldLookup_shader.createUniform("nydy");
  GLint uNumInRow = windFieldLookup_shader.createUniform("numInRow");
  
  windFieldLookup_shader.activate();
  
  glUniform1iARB(unxdx, nxdx);
  glUniform1iARB(unydy, nydy);
  glUniform1iARB(uNumInRow, numInRow);
  
  windFieldLookup_shader.deactivate();
  
  // Setup the result texture
  glGenTextures(1, &windFieldLookupTex);
  glBindTexture(GL_TEXTURE_2D, windFieldLookupTex);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, 1, 1, 0, GL_RGBA, GL_FLOAT, NULL);

  // Bind it to the FBO
  windFieldLookupFBO = new FramebufferObject();
  windFieldLookupFBO->AttachTexture(GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, windFieldLookupTex);
  windFieldLookupFBO->IsValid();
  FramebufferObject::Disable();
}

void ParticleControl::lookupWindField(GLfloat pos[3], GLuint windField, GLfloat &wx, GLfloat &wy, GLfloat &wz) {
  windFieldLookupFBO->Bind();
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(-1, 1, -1, 1);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glViewport(0, 0, 1, 1);
  
  glEnable(texType);
  windFieldLookup_shader.activate();
  
  CheckErrorsGL("PC: lookupWindField() - activated lookup shader");
  
  // pass the position into the shader
  glUniform3fARB(uniform_wfl_pos, pos[0], pos[1], pos[2]);
  
  // wind field can be stored here
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(texType, windField);
  glUniform1iARB(uniform_wfl_wind, 0);

  glBegin(GL_QUADS);
  {
    glVertex3f(-1, -1, -0.5f);
    glVertex3f( 1, -1, -0.5f);
    glVertex3f( 1,  1, -0.5f);
    glVertex3f(-1,  1, -0.5f);
  }
  glEnd();
  

  windFieldLookup_shader.deactivate();

  CheckErrorsGL("PC: lookupWindField() - shader advection completion");
  
  GLfloat tmp[4];
  
  glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
  glReadPixels(0, 0, 1, 1, GL_RGBA, GL_FLOAT, &tmp);
  
  FramebufferObject::Disable();
  
  wx = tmp[0];
  wy = tmp[1];
  wz = tmp[2];
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
     // Balli had set this to true...

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

  // These values can get large, much larger than ints, so this value needs to be a long at the least!
  long arrSize = nxdx*nydy*nzdz;
  try
    {
      wind_vel = new wind[arrSize];
      cellQuic = new cellType[arrSize];
    }
  catch (std::bad_alloc&)
    {
      std::cout << std::endl;
      std::cout << "*****************************************************************************************************" << std::endl;
      std::cout << "Error allocating memory for the wind and cell type textures.  Attempted to allocate " << arrSize << " floats." << std::endl;
      std::cout << "\tDomain Specs: nx=" << nx << ", ny=" << ny << ", nz=" << nz << std::endl;
      std::cout << "\t\tcell spacing = (" << cell_dx << ", " << cell_dy << ", " << cell_dz << ")" << std::endl;
      std::cout << "\t\total cells = " << nxdx * nydy * nzdz << std::endl;
      std::cout << "*****************************************************************************************************" << std::endl;
      std::cout << std::endl;
    }

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
  float znaut = 0.1f;
  
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
    {
      std::cerr << "gpuPlume Error: QP_turbield.dat file must exist in the project!" << std::endl;
      exit(EXIT_FAILURE);
    }

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

    for (qi=0; qi<nydy; qi++)
      for (qj=0; qj<nxdx; qj++)
  for (qk=0; qk<nzdz; qk++) 
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

	    // The following is unused...
	    // float tauDetInv=1.0f/((tau11*tau22*tau33)-(tau13*tau13*tau22));
	  
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

	    dataTwo[texidx]   =  1.0f/(tau11-tau13*tau13/tau33);// tauDetInv*(tau22*tau33);      //Lam11 //changed Balli
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

    //createTexture(lambda, GL_RGBA32F_ARB, width,height, dataTwo);
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

  std::cout << "---------------------------------" << std::endl;
  std::cout << "initLambda_and_Taus_withCalculations!" << std::endl;
  std::cout << "---------------------------------" << std::endl;

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

  float dz = 1.0;
  // The following should no longer be needed... not sure why dz is
  // still needed; need to check.
  // ,dx=1.0,dy=1.0; //Hardwired** should be in input file

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

	  // The following is unused
	  // float tauDetInv=1.0f/((tau11*tau22*tau33)-(tau13*tau13*tau22));
		
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
	    pow((1.0f-0.85f*minDistance/25.0f),(1.5f))/(0.4*(minDistance));//This value is the '0.5*CoEps' value	

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
  
  std::cout << "Calling getMinDistance!!!!!!!!!!!!!!!!!!!!!!!!! should be fixed!" << std::endl;

  //Following is hardwired, but we need to read the building info mentioned in the input file
  int totBuild=6;
  int xfo[]={20,20,35,35,50,50};//{25};
  int yfo[]={18,32,18,32,18,32};//{25};
  // int zfo[]={0,0,0,0,0,0};//{0};
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
    // int z=zfo[build];
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
	wind_vel[p2idx].u = 7.52f*pow(((k+1)/25.0f),0.15f);
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
  in.open(path.c_str(),std::ios::binary|std::ios::in);
    	
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
                  if(numSides[n]==4)
                  {
                     float lk = zfo[n];
                     float uk = zfo[n]+ht[n];
                     float li = yfo[n]-(wti[n]/2.0);
                     float ui = yfo[n]+(wti[n]/2.0);
                     float lj = xfo[n];
                     float uj = xfo[n]+lti[n];

                     if(x >= lj && x <= uj && y >= li && y <= ui && z >= lk && z <= uk) cell_type[p2idx].id = float(n);
                      
                  }
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

 /*for(int n=0; n < numBuild; n++){
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
   }*/
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
	
  // for now, only read the ascii files... binary isn't yet done.
  bool asciiFile = true;

  std::string path;
  if(quicFilesPath.c_str() != NULL){
    if (asciiFile)
      path = quicFilesPath + "QU_velocity.dat";
    else 
      path = quicFilesPath + "QU_velocity.bin";
  }
  else
    path = "Settings/QU_velocity.dat";

  if (asciiFile)
    QUICWindField.open(path.c_str()); //opening the wind file  to read
  else
    QUICWindField.open(path.c_str(), std::ios::binary); //opening the wind file  to read

  if(!QUICWindField){
    std::cerr<<"Unable to open QUIC Windfield file : QU_velocity.dat ";
    exit(1);
  }

  if (asciiFile)
    {
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
	    
	    // std::cout << "v = (" << wind_vel[p2idx].u << ", " << wind_vel[p2idx].v << ", " << wind_vel[p2idx].w << ")" << std::endl;
	    /*QUICCellType>>quicIndex;// ignoring the X,Y and Z values
	      QUICCellType>>quicIndex;
	      QUICCellType>>quicIndex;

	      QUICCellType>>cellQuic[p2idx].c ;//storing the Celltype values in the Cell structure*/
		
	  }
	}
      }
    }
  else 
    {
      // Binary read!  In the binary file, comments are NOT written.
      // read nxXnyXnz floats All u, All v, All w is the format.

      QUICWindField.seekg(0, std::ios::end);
      int numBytes = QUICWindField.tellg();
      QUICWindField.seekg(0, std::ios::beg);

      // allocate memory for whole file and read everything in at once
      char *fileBuffer = new char[numBytes];
      QUICWindField.read(fileBuffer, numBytes);
      QUICWindField.close();

      // skip over 4 or 8 bytes since that's how fortran dumps binary
      // writes (depends on 32-bit versus 64-bit, respectively and I'm not
      // sure how to determine that yet)  ???
      int cd = 0;
      cd += sizeof(int);

      int dataN = nxdx * nydy * nzdz;
      int dataSz = dataN * sizeof(double);

      double *dval;
      int p2idx = 0;

#if 0
      for(int k = 0; k < nzdz; k++)
	for(int i = 0; i < nydy; i++)
	  for(int j = 0; j < nxdx; j++)
#endif
	    for(int p2idx = 0; p2idx<dataN; p2idx++)
	    {
	      // p2idx = k*nxdx*nydy + i*nxdx + j;
	      //for (int i=0; i<(dataN); i++)
	      //	{
	      dval = reinterpret_cast<double*>(&fileBuffer[cd]);
	      // u[i] = *dval;
	      wind_vel[p2idx].u = *dval;
	  
	      dval = reinterpret_cast<double*>(&fileBuffer[cd + dataSz]);
	      // v[i] = *dval;
	      wind_vel[p2idx].v = *dval;
	  
	      dval = reinterpret_cast<double*>(&fileBuffer[cd + dataSz*2]);
	      // w[i] = *dval;
	      wind_vel[p2idx].w = *dval;
	  
	      // std::cout << "v = (" << wind_vel[p2idx].u << ", " << wind_vel[p2idx].v << ", " << wind_vel[p2idx].w << ")" << std::endl;

	      cd += sizeof(double);
	}

      // skip over 4 bytes since that's how fortran dumps binary writes 
      // cd += 8;
      
      delete [] fileBuffer;
    }

  long double tot_vel=0.0;
  long double vel=0.0;
  long double mx=0.0;
  long double my=0.0;
  long double mz=0.0;
  long double avg_vel=0.0;
  long double sd=0.0;
  long double sum_of_squares=0.0;
  long double local_max_vel = -1.0;

  // BIG QUESTION from Pete: why do I need to have nzdz - 1 here to get numbers that are correct for the wind field.   I didn't need this on
  // Linux so there must be some bug in how we're interpreting the nzdz K layer... someone needs to look into this... note that this is only
  // affecting the visualization at the moment, but I'd like to understand what's happening!!!
  // ???????????
  // for(int k = 0; k < nzdz-1; k++){   

  for(int k = 0; k < nzdz; k++){   
    for(int i = 0; i < nydy; i++){
      for(int j = 0; j < nxdx; j++){
	int p2idx = k*nxdx*nydy + i*nxdx + j;

	assert(p2idx < (nxdx * nydy * nzdz));
	assert(p2idx >= 0);

	double usq = (double)wind_vel[p2idx].u * (double)wind_vel[p2idx].u;
	double vsq = (double)wind_vel[p2idx].v * (double)wind_vel[p2idx].v;
	double wsq = (double)wind_vel[p2idx].w * (double)wind_vel[p2idx].w;

	vel = sqrt(usq + vsq + wsq);

	if (vel > local_max_vel) 
	  {
	    // std::cout << "nxdx = " << nxdx << ", nydy = " << nydy << ", msize = " << nxdx*nydy*nzdz << ", p2idx = " << p2idx << ", usq=" << usq << ", vsq=" << vsq << ", wsq=" << wsq << ", vel = " << vel << ", max vel = " << local_max_vel << std::endl;
	    local_max_vel = vel;
	    mx = wind_vel[p2idx].u;
	    my = wind_vel[p2idx].v; 
	    mz = wind_vel[p2idx].w;
	  }
	tot_vel += vel;

      }
    }
  }

  avg_vel = tot_vel/(double)(nxdx*nydy*nzdz);
  // for(int k = 0; k < nzdz-1; k++){   
  for(int k = 0; k < nzdz; k++){   
    for(int i = 0; i < nydy; i++){
      for(int j = 0; j < nxdx; j++){
	int p2idx = k*nxdx*nydy + i*nxdx + j;
	
	vel = sqrt(wind_vel[p2idx].u * wind_vel[p2idx].u + 
		   wind_vel[p2idx].v * wind_vel[p2idx].v + 
		   wind_vel[p2idx].w * wind_vel[p2idx].w);
	sum_of_squares += ((vel-avg_vel) * (vel-avg_vel));
	}
    }
   }

  sd=sqrt(sum_of_squares/(double)(nxdx*nydy*nzdz));

  std::cout << "Calculated maximum wind velocity=" << local_max_vel << "\n\taverage wind velocity=" << avg_vel << "\n\tstandard deviation (SD)=" << sd << std::endl;
  if (avg_vel + 3.0*sd < local_max_vel)
    {
      m_util_ptr->calculatedMaxVel = avg_vel + 3.0*sd;
      std::cout << "\tsetting maximum wind velocity to 3*SD from mean = " << m_util_ptr->calculatedMaxVel << std::endl;
    }
  else 
    {
      m_util_ptr->calculatedMaxVel = local_max_vel;
    }

  QUICWindField.close();
}

void ParticleControl::initCellType(){
  std::ifstream QUICCellType;
  std::string path;

  // NOTE!  celltypes are smallish integers!  Might be able to save
  // memory by using basic byte, not 32-bit floating textures

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

#if 0
  // this appears to do NOTHING! Commenting out for now... -Pete
  std::cout << "Local Max t11 values" << std::endl;
  for(int i=0; i < nzdz; i++){
    //std::cout << tauLocalMax[(i*4)+3] << std::endl;
    
  }
#endif
}

void ParticleControl::nonLocalMixing(GLuint windField,GLuint lambda, GLuint tau_dz, GLuint duvw_dz, GLuint dxyz_wall, 
                                     GLuint tauTex){

  std::cout << "Running with Non-Local Mixing Calculations" << std::endl;

    GLfloat *data     = new GLfloat[ width * height * 4 ];
    GLfloat *dataWind = new GLfloat[ width * height * 4 ];
    GLfloat *dataTwo  = new GLfloat[ width * height * 4 ];
    GLfloat *dataTau  = new GLfloat[width*height*4];
    GLfloat *data3    = new GLfloat[width*height*4];
    GLfloat *data4    = new GLfloat[width*height*4];

    
#if 0
    //
    // Create a VBO to hold the data for drawing arrows
    //
    
    GLfloat *windVBOdata = new GLfloat[ width * height * 3];  // Enough to hold the wind field direction vectors AND the center position of each cell
#endif
    
    
    //Balli's new additions this is essentially a direct copy of the FORTRAN
    initCellType();
    
    //Balli: Substracting 1 from nzdz as it is increased by 1  after reading from QU_simparams.inp in Util.cpp
    nzdz=nzdz-1;
    float nxnynz=nxdx*nydy*nzdz;
    std::string s;
    float dx = m_util_ptr->dx;
    float dy = m_util_ptr->dy;
    float dz = m_util_ptr->dz;
    std::vector<float> dz_array,z,zm;
    dz_array.resize(nzdz,dz);
    z.resize(nzdz,0.0f);
    zm.resize(nzdz,0.0f);
    //Balli: Initialized first element of z and zm array before begining the loop as in GPU Plume
    //we do not store values below the ground, which are zero anyways
    z.at(0)  = dz_array.at(0);
    zm.at(0) = z.at(0)-0.5*dz_array.at(0);
    
    for(int k=1;k<nzdz;k++){
        z.at(k)  = z.at(k-1)+dz_array.at(k);
        zm.at(k) = z.at(k)-0.5*dz_array.at(k);
    }
    //**************************ATTENTION********************************************************
    //Balli: Following should be grabbed from the input files: Hardwired here !!!!     IMPORTANT!!!
    float rcl  = 0.0f;  //Monin-obo length, should be in the input file-Balli-06/10/09
    float z0   = 0.1f;  //Should be in the input file-Balli-06/10/09
    int roofflag = 2;   //should be in the input file-Balli-06/10/09
    float h      = 25.f; //Boundary Layer Height-should be in the input file-Balli-06/10/09
    //**************************ATTENTION********************************************************

    
    //Balli: declaring few constants
    const float kkar = 0.4f;           //von karman constant
    const float pi   = 4.f*atan(1.0f);
    const float knlc   = 0.113f;         
    const float ctau13 = 1.f;           
    const float cusq   = 2.5f*2.5f;     
    const float cvsq   = 2.f*2.f;         
    const float cwsq   = 1.3f*1.3f;    


    //Balli: "theta" is never read into or initilized in the FORTRAN code, but used for calculating ualoft and valoft.
    float theta  = 0.f;          // This varibale is not used anymore in QP but the legacy code still uses it-Balli-06/10/09
    float ualoft = 0.f;          
    float valoft = 0.f;             

    int time_idx=1;
    int check1=0;
    int check2=0;
    int check3=0;
    int check=0;

    //Balli: For writing turbulence data- Can be removed later
    std::ofstream turbfield;
    turbfield.open("GPU_turbfield.dat");
    
    //Balli : Declaring local vectors
    std::vector<float> elz,ustarz,sigwi,sigvi,ustarij,xi,yi,zi,hgt,hgtveg,eleff,xcb,ycb,icb,jcb,phib,weff,leff,lfr,zcorf,lr;
    std::vector<float>uref,urefu,urefv,urefw, utotktp,uktop,vktop,wktop,deluc,ustargz,elzg,ustarg;
    std::vector<float>utotcl1,utotmax,gamma,atten,Sx,Sy;
    std::vector<int>bldtype;
    
    //Balli : Vectors needs to be resized before they are  used otherwise they sometime give runtime errors
    eleff.resize(nxnynz,0.0);// efective length scale-initialized with zero values.
    ustarg.resize(nxnynz,0.0);


    
    //Balli : Reading "QU_buildout.inp"; THis should be handled in Util.cpp
    //This file is required for getiign effective lenths of downstream and upstream(along with other parameters) cavities from QUICURB
    std::ifstream QUbldout;
    char char18[1024];
    QUbldout.open("QP_buildout.inp");
    int inumveg,inumgarage;
    //Balli: IMPORTANT!!!!
    
    //OVERWRITING!!! the building parameters here and reading them from QP_bldout file. The values of the key parameters(xfo,yfo,zfo,ht,lti,wti)
    //will remain same as bldout file has all the information contained in QU_buildings with some additional information needed for
    //non-local mixing
    QUbldout>>numBuild;
    bldtype.resize(numBuild);
    gamma.resize(numBuild);
    atten.resize(numBuild);
    Sx.resize(numBuild);
    Sy.resize(numBuild);
    weff.resize(numBuild);
    leff.resize(numBuild);
    lfr.resize(numBuild);
    lr.resize(numBuild);
    std::getline(QUbldout,s);
    QUbldout>>inumveg;
    std::getline(QUbldout,s);
    
    for(int i=0;i<numBuild-inumveg;i++){
        QUbldout>>char18>>char18>>char18;
        int buildnum;
        QUbldout>>buildnum;
        QUbldout>>char18>>char18;
        QUbldout>>bldtype.at(i);
        QUbldout>>char18>>char18;
        QUbldout>>gamma.at(i);
        QUbldout>>char18>>char18;
        QUbldout>>ht[i];
        QUbldout>>char18>>char18;
        QUbldout>>wti[i];
        QUbldout>>char18>>char18;
        QUbldout>>lti[i];
        QUbldout>>char18>>char18;
        QUbldout>>xfo[i];
        QUbldout>>char18>>char18;
        QUbldout>>yfo[i];
        QUbldout>>char18>>char18;
        QUbldout>>zfo[i];
        QUbldout>>char18>>char18;
        QUbldout>>weff.at(i);
        QUbldout>>char18>>char18;
        QUbldout>>leff.at(i);
        QUbldout>>char18>>char18;
        QUbldout>>lfr.at(i);
        QUbldout>>char18>>char18;
        QUbldout>>lr.at(i);
        QUbldout>>char18>>char18;
        QUbldout>>atten.at(i);
        QUbldout>>char18>>char18;
        QUbldout>>Sx.at(i);
        QUbldout>>char18>>char18;
        QUbldout>>Sy.at(i);
        QUbldout>>char18>>char18>>char18>>char18;
    }
    QUbldout.close();
    
    //Balli: IMPORTANT!!
    //QP differs from GPU in the indices i,j,k of all the arrays (u,v,w etc.)
    //Following calculations are for Boundary layer case only and it provides an insight into the coordinate sytem differences
    //QP velocity vectors, for example u(i,j,k), i goes from 1(0.5) to nx(9.5), j goes from 1(0.5) to ny(9.5)
    //and k goes from 1(-0.5) to nz+1(29.5)
    //[Note: nx,ny,nz above are what QP reads from input file, QP adds 1 to nx and ny, and 2 to nz after reading them from input file]
    
    //QP's k goes from 1 to nz+2, therefore for dz=1, zi goes from -0.5(k=1) to 30.5(k=32)
    //GPU k goes from 0 to nz-1, therefore for dz=1, zi goes from 0.5(k=0) to 29.5 (k=29)

    zi.resize(nzdz);
    for(int k=0;k<nzdz;k++){ 
        zi.at(k)=.5*dz+dz*k; //this expression is different from that used in the QP, but values are same
    }
    //QP's j goes from 1 to ny+1, therefore for dy=1, yi goes from 0.5(j=1) to 10.5(j=11)
    //GPU j goes from 0 to ny-1, therefore for dy=1, yi goes from 0.5(j=0) to 9.5 (j=9)
    
    yi.resize(nydy);
    for(int j=0;j<nydy;j++){
        yi.at(j)=.5*dy+dy*j;
    }
    //QP's i goes from 1 to nx+1, therefore for dx=1, xi goes from 0.5(i=1) to 10.5(i=11)
    //GPU i goes from 0 to nx-1, therefore for dx=1, xi goes from 0.5(i=0) to 9.5 (i=9)
    
    xi.resize(nxdx);
    for(int i=0;i<nxdx;i++){
        xi.at(i)=.5*dx+dx*i;
    }
    
    float ht_avg = 0.0;
    int k=0;
    if(numBuild > 0 && numBuild != inumveg){
        for(int  i_b=0;i_b<numBuild;i_b++){
            if(bldtype.at(i_b)==9){
                continue;
            }
            ht_avg=ht[i_b]+zfo[i_b]+ht_avg;
        }
        ht_avg=ht_avg/float(numBuild-inumveg);
        float temp=ht_avg/dz;
        for(int kk=0;kk<nzdz;kk++){
            k=kk;
            if(ht_avg<z.at(kk))break;
        }
    }
    else{
        //BL Flow case: Control comes here as we have no buildings for this test case
        k=0; //altered to comply with GPU
    }
    dz=dz_array.at(k);
    
    //Obtain avg. velocity from the boundary of the domain at above cal avg ht of the buildings
    int i=0;//altered to comply with GPU
    int j=0;
    float u_left=0;
    for(j=0;j<nydy;j++){
        int p2idx = k*nxdx*nydy + j*nxdx + i;
        u_left=sqrt(wind_vel[p2idx].u*wind_vel[p2idx].u + wind_vel[p2idx].v*wind_vel[p2idx].v + wind_vel[p2idx].w*wind_vel[p2idx].w) +u_left;
    }
    u_left=u_left/(nydy);//altered to comply GPU
    // in QP, total number of cells in y is read from input file as ny and then QP adds 1to ny, therefore, 1 is substracted from ny above in QP.
    
    j=nydy-1;// substracted 1 as edge of the domain in y is nydy-1
    float u_top=0;
    for(i=0;i<nxdx;i++){
        int p2idx = k*nxdx*nydy + j*nxdx + i;
        u_top=sqrt(wind_vel[p2idx].u*wind_vel[p2idx].u + wind_vel[p2idx].v*wind_vel[p2idx].v + wind_vel[p2idx].w*wind_vel[p2idx].w) +u_top;
    }
    u_top=u_top/nxdx;
    i=nxdx-1;//same explanation as in case of j above
    float u_right=0;
    for(j=0;j<nydy;j++){
        int p2idx = k*nxdx*nydy + j*nxdx + i;
        u_right=sqrt(wind_vel[p2idx].u*wind_vel[p2idx].u + wind_vel[p2idx].v*wind_vel[p2idx].v + wind_vel[p2idx].w*wind_vel[p2idx].w) +u_right;
    }
    u_right=u_right/nydy;
    j=0;//alterted for GPU
    
    float u_bottom=0;
    for(i=0;i<nxdx;i++){
        int p2idx = k*nxdx*nydy + j*nxdx + i;
        u_bottom=sqrt(wind_vel[p2idx].u*wind_vel[p2idx].u + wind_vel[p2idx].v*wind_vel[p2idx].v + wind_vel[p2idx].w*wind_vel[p2idx].w) +u_bottom;
    }
    u_bottom=u_bottom/(nxdx);
    
    float u_b=0.25*(u_left+ u_top+ u_right+ u_bottom);//average velocity
    float nu_b=1.5e-5; //nu for air
    float del_b=(0.328* pow(nu_b/u_b,.2f) ) * pow(ht_avg,.8f);//expression for BL layer thickness (growth) (i think)
    // above expression is used for obtaining turbulence close to the walls of the buildings
    
    hgt.resize(nxnynz,0.);
    hgtveg.resize(nxnynz,0.);
    if(time_idx == 1){
        for(int k=0; k<nzdz;k++){
            for(int j=0; j<nydy;j++){
                for(int i=0; i<nxdx;i++){
                    int p2idx = k*nxdx*nydy + j*nxdx + i;
                    int ij = j*nxdx + i;
                    if(cellQuic[p2idx].c == 0)hgt.at(ij)=std::max(hgt.at(ij),z.at(k));
                    if(cellQuic[p2idx].c == 8)hgtveg.at(ij)=std::max(hgtveg.at(ij),z.at(k));
                }
            }
        }
    }    
    
    elz.resize(nxnynz);
    ustarz.resize(nxnynz);
    sigwi.resize(nxnynz);
    sigvi.resize(nxnynz);
    ustarij.resize(nxnynz);
    ustarz.resize(nxnynz);

    
    //Balli:Allocating global arrays, declared in header file
    dutotdxi.resize(nxnynz);
    dutotdyi.resize(nxnynz);
    dutotdzi.resize(nxnynz,0.0);
    dutotdni.resize(nxnynz);
    dutotdsi.resize(nxnynz);
    alph1ij.resize(nxnynz);
    alph2ij.resize(nxnynz);
    alph3ij.resize(nxnynz);
    bet1ij.resize(nxnynz);
    bet2ij.resize(nxnynz);
    bet3ij.resize(nxnynz);
    gam1ij.resize(nxnynz);
    gam2ij.resize(nxnynz);
    gam3ij.resize(nxnynz);
    alphn1ij.resize(nxnynz);
    alphn2ij.resize(nxnynz);
    alphn3ij.resize(nxnynz);
    betn1ij.resize(nxnynz);
    betn2ij.resize(nxnynz);
    betn3ij.resize(nxnynz);
    gamn1ij.resize(nxnynz);
    gamn2ij.resize(nxnynz);
    gamn3ij.resize(nxnynz);
    ani.resize(nxnynz);
    bni.resize(nxnynz);
    cni.resize(nxnynz);
    ufsqgi.resize(nxnynz);
    vfsqgi.resize(nxnynz);
    wfsqgi.resize(nxnynz);
    ufvfgi.resize(nxnynz);
    ufwfgi.resize(nxnynz);
    vfwfgi.resize(nxnynz);

    
    std::vector<float> dzm,dzp,dym,dyp,dxm,dxp,ufwfi,ufvfi,vfwfi,sigui,upwpi,epsi;
    
    dzm.resize(nzdz*nydy*nxdx);
    dzp.resize(nzdz*nydy*nxdx);
    dxm.resize(nzdz*nydy*nxdx);
    dxp.resize(nzdz*nydy*nxdx);
    dym.resize(nzdz*nydy*nxdx);
    dyp.resize(nzdz*nydy*nxdx);
    ufwfi.resize(nzdz*nydy*nxdx);
    ufvfi.resize(nzdz*nydy*nxdx);
    vfwfi.resize(nzdz*nydy*nxdx);
    sigui.resize(nzdz*nydy*nxdx);
    upwpi.resize(nzdz*nydy*nxdx);
    epsi.resize(nzdz*nydy*nxdx);
    
    int kcantop=0;
    float ucantop=0.;
    float elcanopy=0.;
    float utotl=0.;
    float utotu=0.;
    float phim=0.;
    float psim=0.;
    float dutotl=0.;
    float dutotu=0.f;
    float dutot=0.f;
    float utot=0.f;
    float dutotdzc=0.;
    float dutotdzp=0.;
    float dutotdzm=0.;
    float dutotdza=0.;
   
    for(int j=0;j<nydy;j++){
        for(int i=0;i<nxdx;i++){
            int ij = j*nxdx + i;
            if(hgtveg.at(ij) > 0.){
                for(int kk=0;kk<nzdz;k++){
                    kcantop=kk;
                    if(hgtveg.at(ij) <= z.at(kk))break;
                }
                int idcan  = kcantop*nxdx*nydy + j*nxdx + i;
                int id1can = (kcantop+1)*nxdx*nydy + j*nxdx + i;
                
                ucantop=.5*sqrt(wind_vel[idcan].u*wind_vel[idcan].u+wind_vel[idcan].v*wind_vel[idcan].v+wind_vel[idcan].w*wind_vel[idcan].w)
                    + .5*sqrt(wind_vel[id1can].u*wind_vel[id1can].u+wind_vel[id1can].v*wind_vel[id1can].v+wind_vel[id1can].w*wind_vel[id1can].w); 
            }
            for(int k=0;k<nzdz;k++){

                dz=dz_array.at(k);
                int km1   = (k-1)*nxdx*nydy + j*nxdx + i;
                int kp1   = (k+1)*nxdx*nydy + j*nxdx + i;
                int knz1 = (nzdz-1)*nxdx*nydy + j*nxdx + i;
                int p2idx = k*nxdx*nydy + j*nxdx + i;
                int ij = j*nxdx + i;
                int idklow=0;
                //new changes from QUIC
                if(cellQuic[p2idx].c != 0){
                    dzm.at(p2idx)=zm.at(k)-hgt.at(ij);
                    eleff.at(p2idx)=dzm.at(p2idx);
                }
                else{
                    dzm.at(p2idx)=0.f;
                    eleff.at(p2idx)=0.f;
                }
                elcanopy=0.f;
                if(((cellQuic[km1].c == 0) || (cellQuic[km1].c==8)) && 
                   (cellQuic[p2idx].c != 0 && cellQuic[p2idx].c != 8) || k == 0){//altered k
                    utotl=0.f;
                    utotu=sqrt(wind_vel[p2idx].u*wind_vel[p2idx].u+wind_vel[p2idx].v*wind_vel[p2idx].v+wind_vel[p2idx].w*wind_vel[p2idx].w);
                    //MDW 7-01-2005 changed the way vertical gradients are calculated to avoid inaccuracies
                    // in the representation of the gradients of a log-law term
                    
                    if(rcl>0){
                        phim=1.f+4.7f*rcl*.5f*dz;
                        psim=-4.7f*rcl*.5f*dz;
                    }
                    else{
                        phim=pow( (1.f-15.f*rcl*.5f*dz),(-.25f));
                        psim=2.f*log((1.f+1.f/phim)/2.f)+log((1.f+1.f/pow(phim,2.f))/2.f)-2.f*atan(1.f/phim)+pi/2.f;
                    }
                    if(hgtveg.at(ij)>zi.at(k)){                           
                        phim=1.f;
                        psim=0.f;
                        elz.at(p2idx)=elcanopy*std::min(1.f,(zi.at(k)-z0)/(.3f*hgtveg.at(ij)));
                        ustar=elz.at(p2idx)*utotu/(.5f*dz);
                        dutotdzi.at(p2idx)=utotu/(.5f*dz);
                        ustarz.at(p2idx)=ustar;
                    }
                    else{
                        if(cellQuic[km1].c!=8){
                            ustar=kkar*utotu/(log(.5f*dz/z0)-psim);
                            elz.at(p2idx)=kkar*.5f*dz;
                            ustarz.at(p2idx)=ustar;
                            dutotdzi.at(p2idx)=ustar*phim/elz.at(p2idx);
                        }
                        else{
                            utotl=sqrt(wind_vel[km1].u*wind_vel[km1].u+wind_vel[km1].v*wind_vel[km1].v+wind_vel[km1].w*wind_vel[km1].w);
                            dutotdzi.at(p2idx)=2.f*(utotu-utotl)/(dz_array.at(k-1)+dz_array.at(k));
                            elz.at(p2idx)=kkar*.5f*dz;
                            ustar=elz.at(p2idx)*dutotdzi.at(p2idx);
                            ustarz.at(p2idx)=ustar;
                        }
                    }
                    if(cellQuic[km1].c!=8 && k!=0){
                        sigwi.at(km1)=0.f;
                        sigvi.at(km1)=0.f;
                        ustarij.at(km1)=0.f;
                        ustarz.at(km1)=0.f;
                    }
                }
                else{
                    if(k==nzdz-1){ // find gradient using a non-CDD approach
                        utotl=sqrt(wind_vel[km1].u*wind_vel[km1].u+wind_vel[km1].v*wind_vel[km1].v+wind_vel[km1].w*wind_vel[km1].w);
                        utotu=sqrt(wind_vel[p2idx].u*wind_vel[p2idx].u+wind_vel[p2idx].v*wind_vel[p2idx].v+wind_vel[p2idx].w*wind_vel[p2idx].w);
                        dutotdzi.at(knz1)=dutotdzi.at(km1)*zm.at(k-1)/zm.at(k);
                        elz.at(p2idx)=kkar*(eleff.at(p2idx)-hgtveg.at(ij));
                    }
                    else{ // find gradient using a CDD approach
                        utotl=sqrt(wind_vel[km1].u*wind_vel[km1].u+wind_vel[km1].v*wind_vel[km1].v+wind_vel[km1].w*wind_vel[km1].w);
                        utotu=sqrt(wind_vel[kp1].u*wind_vel[kp1].u+wind_vel[kp1].v*wind_vel[kp1].v+wind_vel[kp1].w*wind_vel[kp1].w);
                        // mdw 7-08-2005 changed the way vertical gradients are calculated to better represent
                        // log-law behavior
                        if(cellQuic[p2idx].c==8){
                            dutotdzi.at(p2idx)=(utotu-utotl)/(dz_array.at(k)+.5*dz_array.at(k-1)+.5*dz_array.at(k+1));
                            ustarz.at(p2idx)=elcanopy*dutotdzi.at(p2idx);
                            elz.at(p2idx)=elcanopy*std::min(1.f,(zi.at(k)-z0)/(0.3f*hgtveg.at(ij)));
                        }
                        else{
                            int klow=0;
                            for (int kk=0;kk<nzdz;kk++){
                                klow=kk;
                                if(std::max(hgt.at(ij),hgtveg.at(ij))<z.at(kk))break;
                            }

                            idklow = klow*nxdx*nydy + j*nxdx + i;
                            if(rcl>0){
                                phim=1.f+4.7f*rcl*(eleff.at(km1)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f));
                                psim=-4.7f*rcl*(eleff.at(km1)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f));
                            }
                            else{
                                phim=pow( (1.f-15.f*rcl*(eleff.at(km1)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f))),(-.25f));
                                psim=2.f*log((1.f+1.f/phim)/2.f)+log((1.f+1.f/pow(phim,2.f))/2.f)-2.f*atan(1.f/phim)+pi/2.f;
                            }
                            dutotl=utotl-ustarz.at(idklow)*(log((zi.at(k-1)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f))/z0)-psim)/kkar;
                            if(rcl>0){
                                phim=1.f+4.7f*rcl*(eleff.at(kp1)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f));
                                psim=-4.7f*rcl*(eleff.at(kp1)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f));
                            }
                            else{
                                phim=pow( (1.f-15.f*rcl*(eleff.at(kp1)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f))),(-.25f));
                                psim=2.f*log((1.f+1.f/phim)/2.f)+log((1.f+1.f/pow(phim,2.f))/2.f)-2.f*atan(1.f/phim)+pi/2.f;
                            }
                            dutotu=utotu-ustarz.at(idklow)*(log((zi.at(k+1)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f))/z0)-psim)/kkar;
                            dutotdzi.at(p2idx)=(dutotu-dutotl)/(dz_array.at(k)+.5*dz_array.at(k-1)+.5*dz_array.at(k+1))
                                +ustarz.at(idklow)*psim/(kkar*zi.at(k));
                            elz.at(p2idx)=kkar*(eleff.at(p2idx)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f));
                              if(cellQuic[kp1].c != 0 && cellQuic[p2idx].c != 0  && cellQuic[km1].c != 0){
                                  // mdw 7-01-2005 centered around k instead of k-1 and ajusted for log-law behavior
                                  utot=sqrt(wind_vel[p2idx].u*wind_vel[p2idx].u+wind_vel[p2idx].v*wind_vel[p2idx].v+wind_vel[p2idx].w*wind_vel[p2idx].w);
                                  utotl=sqrt(wind_vel[km1].u*wind_vel[km1].u+wind_vel[km1].v*wind_vel[km1].v+wind_vel[km1].w*wind_vel[km1].w);
                                  utotu=sqrt(wind_vel[kp1].u*wind_vel[kp1].u+wind_vel[kp1].v*wind_vel[kp1].v+wind_vel[kp1].w*wind_vel[kp1].w);
                                  if(rcl>0){
                                      phim=1.f+4.7f*rcl*(eleff.at(km1)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f));
                                      psim=-4.7f*rcl*(eleff.at(km1)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f));
                                  }
                                  else{
                                      phim=pow( (1.-15.*rcl*(eleff.at(km1)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f))),(-.25));
                                      psim=2.*log((1.+1./phim)/2.)+log((1.+1./pow(phim,2.f))/2.)-2.*atan(1./phim)+pi/2.;
                                  }
                                  dutotl=utotl-ustarz.at(idklow)*(log((zi.at(k-1)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f))/z0)-psim)/kkar;
                                  if(rcl>0){
                                      phim=1.+4.7*rcl*(eleff.at(kp1)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f));
                                      psim=-4.7*rcl*(eleff.at(kp1)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f));
                                  }
                                  else{
                                      phim=pow( (1.-15.*rcl*(eleff.at(kp1)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f))),(-.25));
                                      psim=2.*log((1.+1./phim)/2.)+log((1.+1./pow(phim,2.f))/2.)-2.*atan(1./phim)+pi/2.;
                                  }
                                  dutotu=utotu-ustarz.at(idklow)*(log((zi.at(k+1)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f))/z0)-psim)/kkar;
                                  // mdw 3-08-2004 begin changes for highest gradient rather than centered diff gradient
                                  if(rcl>0){
                                      phim=1.+4.7*rcl*(eleff.at(p2idx)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f));
                                      psim=-4.7*rcl*(eleff.at(p2idx)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f));
                                  }
                                  else{
                                      phim=pow( (1.-15.*rcl*(eleff.at(p2idx)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f))),(-.25) );
                                      psim=2.*log((1.+1./phim)/2.)+log((1.+1./pow(phim,2.f))/2.)-2.*atan(1./phim)+pi/2.;
                                  }
                                  dutot=utot-ustarz.at(idklow)*(log((zi.at(k)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f))/z0)-psim)/kkar;
                                  dutotdzc=(dutotu-dutotl)/(dz_array.at(k)+.5*dz_array.at(k+1)+.5*dz_array.at(k-1));
                                  dutotdzp=(dutotu-dutot)/(.5*dz_array.at(k+1)+.5*dz_array.at(k));
                                  dutotdzm=(dutot-dutotl)/(.5*dz_array.at(k)+.5*dz_array.at(k-1));
                                  dutotdza=0.5*(fabs(dutotdzp+ustarz.at(idklow)*phim/(kkar*(zi.at(k) 
                                                                                          -std::max((hgtveg.at(ij)-elcanopy/kkar),0.f)))) 
                                                +fabs(dutotdzm+ustarz.at(idklow)*phim/(kkar*(zi.at(k) 
                                                                                           -std::max((hgtveg.at(ij)-elcanopy/kkar),0.f)))));
                                  if(abs(dutotdzp+ustarz.at(idklow)*phim/(kkar*(zi.at(k)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f))))> 
                                     fabs(dutotdzm+ustarz.at(idklow)*phim/(kkar*(zi.at(k)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f))))){
                                      dutotdzi.at(p2idx)=dutotdzp+ustarz.at(idklow)*phim/(kkar*(zi.at(k) 
                                                                                               -std::max((hgtveg.at(ij)-elcanopy/kkar),0.f)));
                                  }
                                  else{
                                      dutotdzi.at(p2idx)=dutotdzm+ustarz.at(idklow)*phim/(kkar*(zi.at(k) 
                                                                                                -std::max((hgtveg.at(ij)-elcanopy/kkar),0.f)));
                                      
                                      
                                  }
                              }
                              // use centered differences away from the boundaries
                        }
                    }
                }
            }
        }
    }//end for loops

    //Balli: Above loop is working as expected. I have matched every variable value with QP (Balli-06/14/09)
    // IMPORTANT!!! Theta is never read from input file. so its values will be zero always.
    //Therefore following few lines do not effect the final solution at all
    // phi will be calculated again by taking into account the actual wind angle at each building.

    float phi=270.-theta;
    phi=phi*pi/180.;
    float cosphi=cos(phi);
    int iupc=0;
    if(cosphi>=0)
        iupc=0;//altered for GPU
    else
        iupc=nxdx;
    
    float sinphi=sin(phi);
    int jupc=0;
    if(sinphi>=0.f)
        jupc=0;//altered for GPU
    else
        jupc=nydy;
    
    float phit=phi+0.5*pi;
    float cosphit=cos(phit);
    float sinphit=sin(phit);
    
    //Following variables are required for non-local mixing.
    float xcelt=0.f;
    float ycelt=0.f;
    int icelt=0;
    int jcelt=0;
    
    float xceln=0.f;
    float yceln=0.f;
    int iceln=0;
    int jceln=0;
    float utott=0.f;
    float delut=0.f;
    float delutz=0.f;
    xcb.resize(numBuild);
    ycb.resize(numBuild);
    icb.resize(numBuild);
    jcb.resize(numBuild);
    phib.resize(numBuild);
    zcorf.resize(nxnynz);
    uref.resize(nxnynz);
    urefu.resize(nxnynz);
    urefv.resize(nxnynz);
    urefw.resize(nxnynz);
    utotktp.resize(nxnynz);
    uktop.resize(nxnynz);
    vktop.resize(nxnynz);
    wktop.resize(nxnynz);
    deluc.resize(nxnynz);
    ustargz.resize(nxnynz);
    elzg.resize(nxnynz,(nxdx+1.)*dx);// initialized it with a value similar to QP
    utotcl1.resize(nxnynz);
    utotmax.resize(nxnynz);

    for(int i=0;i<numBuild;i++){
        if(bldtype.at(i)==9)continue;
        //! mdw 4-16-2004 added proper treatment of zfo
        float temp=0.f;
        int ktop=0;
        for(int k=0;k<nzdz;k++){
            ktop=k;
            if(ht[i]+zfo[i]<z.at(k))break;
        }
        int kmid=0;
        for(int k=0;k<nzdz;k++){        
            kmid=k;
            if(0.5*ht[i]+zfo[i]<z.at(k))break;
        }
        if(bldtype.at(i)==3){
            xcb.at(i)=xfo[i];
        }
        else{
            xcb.at(i)=xfo[i]+.5*lti[i];
        }
        ycb.at(i)=yfo[i];
        temp=(xcb.at(i)-dx)/dx;//substracted dx to comply with GPU
        icb.at(i)=nint(temp);
        temp=(ycb.at(i)-dy)/dy;//substracted dy to comply with GPU
        jcb.at(i)=nint(temp);
        //!mdw 6-05-2005 put in procedure to calculate phi & phit
        int kendv=0;
        if(roofflag==2){
            float Bs=ht[i];
            float BL=wti[i];
            
            if(wti[i]<ht[i]){
                Bs=wti[i];
                BL=ht[i];
            }
            float Rscale = ((pow(Bs,(2.f/3.f)))*(pow(BL,(1.f/3.f))));
            float temp=std::max(.22f*Rscale,.11f*wti[i]);
            float zclim  =std::max(temp,.11f*lti[i]);
            for(int k=0;k<nzdz;k++){        
                kendv=k;
                if(zclim<z.at(k))break;
            }
        }
        else{
            for(int k=0;k<nzdz;k++){        
                kendv=k;
                if(ht[i]+zfo[i]<z.at(k))break;
            }
        }
        kendv=std::min(kendv,nzdz);
        
        int idvel=kendv*nxdx*nydy + jcb.at(i)*nxdx +icb.at(i);
        double tempv=wind_vel[idvel].v;
        double tempu=wind_vel[idvel].u;
        phib.at(i)=atan2(tempv,tempu);
        phi=phib.at(i);
        cosphi=cos(phi);
        int iupc=0;
        if(cosphi>=0)
            iupc=0;//altered for GPU
        else
            iupc=nxdx;
        
        sinphi=sin(phi);
        int jupc=0;
        if(sinphi>=0)
            jupc=0;//altered for GPU
        else
            jupc=nydy;
        
        float phit=phi+0.5*pi;
        cosphit=cos(phit);
        sinphit=sin(phit);
        
        //! ycbp3, and xcbp3 give points 1.5 units outside
        //! of the bldg boundaries to compute reference utot
        float ycbp3=0.f;
        float xcbp3=0.f;
        float ycbm3=0.f;
        float xcbm3=0.f;
        int icbp3=0;
        int icbm3=0;
        int jcbp3=0;
        int jcbm3=0;
        float dycbp3=0.f;
        float dycbm3=0.f;
        float dxcbp3=0.f;
        float dxcbm3=0.f;
        float ycbp=0.f;
        float xcbp=0.f;
        float ycbm=0.f;
        float xcbm=0.f;
        float xcd,ycd,xcu,ycu,xcul,ycul,cosfac;


        if(fabs(sinphit)>=fabs(cosphit)){
            ycbp3=ycb.at(i)+(.5*weff.at(i)+.33*weff.at(i))*sinphit;// ! Get reference values for x,y for non-local mixing
            xcbp3=xcb.at(i)+(.5*weff.at(i)+.33*weff.at(i))*cosphit;// ! 1/3 bldg width outside of building is the boundary for the non-local mixing
            ycbm3=ycb.at(i)-(.5*weff.at(i)+.33*weff.at(i))*sinphit;
            xcbm3=xcb.at(i)-(.5*weff.at(i)+.33*weff.at(i))*cosphit;
            temp=(xcbp3-dx)/dx;
            icbp3=nint(temp);//substracted dx to comply gpu
            temp=(xcbm3-dx)/dx;
            icbm3=nint(temp);//substracted dx to comply gpu
            temp=(ycbp3-dy)/dy;
            jcbp3=nint(temp);//substracted dx to comply gpu
            temp=(ycbm3-dy)/dy;
            jcbm3=nint(temp);//substracted dx to comply gpu
            jcbp3=std::min(jcbp3,nydy-1);
            jcbm3=std::min(jcbm3,nydy-1);
            icbp3=std::min(icbp3,nxdx-1);
            icbm3=std::min(icbm3,nxdx-1);
            jcbp3=std::max(0,jcbp3);//changed from 1 to zeros to comply with gpu
            jcbm3=std::max(0,jcbm3);
            icbp3=std::max(0,icbp3);
            icbm3=std::max(0,icbm3);
            //! searching in the plus y direction for building free flow
            int id=kmid*nxdx*nydy + jcbp3*nxdx +icbp3;
            int jp1=0;
            int jp2=0;
            int isign=0;
            if(cellQuic[id].c == 0){
                if(sinphit>0.f){
                    jp1=jcbp3;
                    jp2=nydy-1;
                    isign=1;
                }
                else{
                    jp1=jcbp3;
                    jp2=0;//altered for GPU
                    isign=-1;
                }
            
                for(int ji=jp1;ji<=jp2;ji=ji+isign){
                    jcbp3=jcbp3+isign;
                    jcbp3=std::min(nydy-1,jcbp3);
                    dycbp3=dy*(jcbp3-1)-ycbp3;
                    ycbp3=dy*(jcbp3-1);
                    xcbp3=xcbp3+cosphit*dycbp3/sinphit;
                    icbp3=int(xcbp3/dx)+1-dx;
                    icbp3=std::min(nx-1,icbp3);
                    //!mdw 34/01/2004 forced indices to be within domain
                    int idMid=kmid*nxdx*nydy + jcbp3*nxdx +icbp3;
                    if(cellQuic[idMid].c!= 0) break;
                }
            }

            //! searching in the minus y direction for building free flow
            int id2=kmid*nxdx*nydy + jcbm3*nxdx +icbm3;
            int jm2=0;
            int jm1=0;
            isign=0;
            if(cellQuic[id2].c == 0){
                if(sinphit>0.f){
                    jm2=0;//altered for GPU;
                    jm1=jcbm3;
                    isign=1;
                }
                else{
                    jm2=nydy-1;
                    jm1=jcbm3;
                    isign=-1;
                }
                for(int ji=jm1;ji>=jm2;ji=ji-isign){// do ji=jm1,jm2,-isign 
                    jcbm3=jcbm3-isign;
                    dycbm3=dy*(jcbm3-1)-ycbm3;
                    ycbm3=dy*(jcbm3-1);
                    xcbm3=xcbm3+cosphit*dycbm3/sinphit;
                    temp=(xcbm3-dx)/dx;
                    icbm3=nint(temp);
                                        
                    jcbp3=std::min(jcbp3,ny-1);
                    jcbm3=std::min(jcbm3,ny-1);
                    icbp3=std::min(icbp3,nx-1);
                    icbm3=std::min(icbm3,nx-1);
                    jcbp3=std::max(0,jcbp3);
                    jcbm3=std::max(0,jcbm3);
                    icbp3=std::max(0,icbp3);
                    icbm3=std::max(0,icbm3);
                    int idMid2=kmid*nxdx*nydy + jcbm3*nxdx +icbm3;
                    if(cellQuic[idMid2].c != 0) break;
                }
            }
            ycbp=ycb.at(i)+(.5*leff.at(i))*sinphi;
            xcbp=xcb.at(i)+(.5*leff.at(i))*cosphi;
            ycbm=ycb.at(i)-(.5*leff.at(i))*sinphi;
            xcbm=xcb.at(i)-(.5*leff.at(i))*cosphi;

            
            if(cosphi>=0.f){
                //! Note the current upstream and downstream limits for the wake non-local mixing
                //! are 3*lr in the downstream direction and lfx upstream in the x direction
                //! and lfy upstream in the y direction
                xcd=xcb.at(i)+(.5*leff.at(i)+.1*dx)*cosphi; // ! get the first point on the center line outside of the building (downstream)
                ycd=ycb.at(i)+(.5*leff.at(i)+.1*dx)*sinphi;// !
                
                //!mdw 7-10-2006 made changes to xcd, ycd,xcu, & ycu - formerly used .5 dx
                if(bldtype.at(i)==3){
                    xcu=xcb.at(i)-(.4*leff.at(i)+dx)*cosphi;// ! (upstream)
                    ycu=ycb.at(i)-(.4*leff.at(i)+dx)*sinphi; //!
                }
                else{
                    xcu=xcb.at(i)-(.5*leff.at(i)+0.1*dx)*cosphi;// ! (upstream)
                    ycu=ycb.at(i)-(.5*leff.at(i)+0.1*dx)*sinphi;// !
                }
                 
                //!mdw 7-05-2006 made changes to xcul & ycul - formerly used .5 dx
                xcul=xcu-(lfr.at(i)+dx)*cosphi;// ! get upper limit of the eddie
                ycul=ycu-(lfr.at(i)+dy)*sinphi;
                xcul=std::max(xcul,0.f);
                xcul=std::min(xcul,dx*(nxdx-1));
                ycul=std::max(ycul,0.f);
                ycul=std::min(ycul,dy*(nydy-1));
                cosfac=1.;
            }
            else{
                //!mdw 7-10-2006 made changes to xcd, ycd,xcu, & ycu - formerly used .5 dx
                xcd=xcb.at(i)+(.5*leff.at(i)+.1*dx)*cosphi;
                ycd=ycb.at(i)+(.5*leff.at(i)+.1*dx)*sinphi;
                if(bldtype.at(i)==3){
                    xcu=xcb.at(i)-(.4*leff.at(i)+dx)*cosphi;// ! (upstream)
                    ycu=ycb.at(i)-(.4*leff.at(i)+dx)*sinphi;// !
                }
                else{  
                    xcu=xcb.at(i)-(.5*leff.at(i)+0.1*dx)*cosphi;// ! (upstream)
                    ycu=ycb.at(i)-(.5*leff.at(i)+0.1*dx)*sinphi;// !
                }
                //!mdw 7-05-2006 made changes to xcul & ycul - formerly used .5 dx
                xcul=xcu-(lfr.at(i)+dx)*cosphi;// ! get upstream limit on the front cavity
                ycul=ycu-(lfr.at(i)+dy)*sinphi;// !
                xcul=std::max(xcul,0.f);
                xcul=std::min(xcul,dx*(nxdx-1));
                ycul=std::max(ycul,0.f);
                ycul=std::min(ycul,dy*(nydy-1));
                cosfac=-1.;
            }
        }
        else{// ! if you are more aligned with y than x
            //! MAN 9/15/2005 use weff and leff appropriately
            ycbp3=ycb.at(i)+(.5*weff.at(i)+.33*weff.at(i))*sinphit;// ! get the effective length of the building
            xcbp3=xcb.at(i)+(.5*weff.at(i)+.33*weff.at(i))*cosphit;
            ycbm3=ycb.at(i)-(.5*weff.at(i)+.33*weff.at(i))*sinphit;
            xcbm3=xcb.at(i)-(.5*weff.at(i)+.33*weff.at(i))*cosphit;
            //! end MAN 9/15/2005
            temp=(xcbp3-dx)/dx;
            icbp3=nint(temp);
            temp=(xcbm3-dx)/dx;
            icbm3=nint(temp);
            temp=(ycbp3-dy)/dy;
            jcbp3=nint(temp);
            temp=(ycbm3-dy)/dy;
            jcbm3=nint(temp);
            jcbp3=std::min(jcbp3,nydy-1);
            jcbm3=std::min(jcbm3,nydy-1);
            icbp3=std::min(icbp3,nxdx-1);
            icbm3=std::min(icbm3,nxdx-1);
            jcbp3=std::max(0,jcbp3);//altered from 1 to zero to comply GPU
            jcbm3=std::max(0,jcbm3);
            icbp3=std::max(0,icbp3);
            icbm3=std::max(0,icbm3);
            //! make sure you are outside of the building !
            int id=kmid*nxdx*nydy + jcbp3*nxdx + icbp3;
            int ip1=0;
            int ip2=0;
            int isign=0;
            
            if(cellQuic[id].c== 0){
                if(cosphit>0){
                    ip1=icbp3;
                    ip2=ny-1;
                    isign=1;
                }
                else{
                    ip1=icbp3;
                    ip2=0;//altered for GPU
                    isign=-1;
                }
                // ! decide which is closest building/floor
				
                for(int ip=ip1;ip<=ip2;ip=ip+isign){//do ip=ip1,ip2,isign 
                    icbp3=icbp3+isign;
                    dxcbp3=dx*(icbp3-1)-xcbp3;
                    xcbp3=dx*((icbp3-1));
                    ycbp3=ycbp3+dxcbp3*sinphit/cosphit;
                    temp=(ycbp3-dy)/dy;
                    jcbp3=nint(temp);
                    jcbp3=std::min(jcbp3,nydy-1);
                    jcbm3=std::min(jcbm3,nydy-1);
                    icbp3=std::min(icbp3,nxdx-1);
                    icbm3=std::min(icbm3,nxdx-1);
                    jcbp3=std::max(0,jcbp3);//altered for GPU
                    jcbm3=std::max(0,jcbm3);//altered for GPU
                    icbp3=std::max(0,icbp3);//altered for GPU
                    icbm3=std::max(0,icbm3);//altered for GPU
                    int idMid=kmid*nxdx*nydy + jcbp3*nxdx + icbp3;
                    if(cellQuic[idMid].c!= 0) break;
                }
            }
            int id2=kmid*nxdx*nydy +jcbm3*nxdx + icbm3;
            
            if(cellQuic[id2].c == 0){
                int im1=0;
                int im2=0;
                isign=0;
                if(cosphit>0.f){
                    im1=icbm3;
                    im2=0;//altered for GPU
                    isign=1;
                }
                else{
                    im1=icbm3;
                    im2=nx-icbm3+1;
                    isign=-1;
                }
                for(int im=im1;im<=im2;im=im+isign){//do im=im1,im2,-isign 
                    icbm3=icbm3-isign;
                    dxcbm3=dx*((icbm3-1))-xcbm3;
                    xcbm3=dx*((icbm3-1));
                    jcbm3=jcbm3+dxcbm3*sinphit/cosphit;
                    jcbp3=std::min(jcbp3,ny-1);
                    jcbm3=std::min(jcbm3,ny-1);
                    icbp3=std::min(icbp3,nx-1);
                    icbm3=std::min(icbm3,nx-1);
                    jcbp3=std::max(0,jcbp3);
                    jcbm3=std::max(0,jcbm3);
                    icbp3=std::max(0,icbp3);
                    icbm3=std::max(0,icbm3);
                    int idMid2=kmid*nxdx*nydy + jcbm3*nxdx +icbm3;
                    if(cellQuic[idMid2].c != 0) break;
                }
            }
            ycbp=ycb.at(i)+(.5*leff.at(i))*sinphit;// !  get back of the building
            xcbp=xcb.at(i)+(.5*leff.at(i))*cosphit;// !
            ycbm=ycb.at(i)-(.5*leff.at(i))*sinphit;// !  get front of the building
            xcbm=xcb.at(i)-(.5*leff.at(i))*cosphit;// !
            if(sinphi>=0.f){
                //! Note the current upstream and downstream limits for the wake non-local mixing
                //    ! are 3*lr in the downstream direction and lfx upstream in the x direction
                //  ! and lfy upstream in the y direction
                //! MAN 9/15/2005 use weff and leff appropriately
                //!mdw 7-05-2006 made changes to xcu,ycu, xcd & ycd - formerly used .5 dy or .5 dx
                xcd=xcb.at(i)+(.5*leff.at(i)+dy)*cosphi;// ! get the first point on the center line outside of the building (downstream)
                ycd=ycb.at(i)+(.5*leff.at(i)+dy)*sinphi;// !
                if(bldtype.at(i)==3){
                    xcu=xcb.at(i)-(.4*leff.at(i)+dx)*cosphi;// ! (upstream)
                    ycu=ycb.at(i)-(.4*leff.at(i)+dx)*sinphi; //!
                }
                else{
                    xcu=xcb.at(i)-(.5*leff.at(i)+0.1*dx)*cosphi;// ! (upstream) 
                    ycu=ycb.at(i)-(.5*leff.at(i)+0.1*dx)*sinphi;// !
                }
                //! end MAN 9/15/2005
                //! mdw 7-05-2006 eliminated .5 dx  or .5 dy in favor of dx & dy
                xcul=xcu-(lfr.at(i)+dx)*cosphi;// ! get upper limit of the eddie
                ycul=ycu-(lfr.at(i)+dy)*sinphi;
                xcul=std::max(xcul,0.f);
                xcul=std::min(xcul,dx*(nxdx-1));
                ycul=std::max(ycul,0.f);
                ycul=std::min(ycul,dy*(nydy-1));
                cosfac=1.f;
            }
            else{
                   //! MAN 9/15/2005 use weff and leff appropriately
                xcd=xcb.at(i)+(.5*leff.at(i)+dy)*cosphi;
                ycd=ycb.at(i)+(.5*leff.at(i)+dy)*sinphi;
                if(bldtype.at(i)==3){
                    xcu=xcb.at(i)-(.4*leff.at(i)+dx)*cosphi;// ! (upstream)
                    ycu=ycb.at(i)-(.4*leff.at(i)+dx)*sinphi;// !
                }
                else{
                    xcu=xcb.at(i)-(.5*leff.at(i)+dx)*cosphi;// ! (upstream) 
                    ycu=ycb.at(i)-(.5*leff.at(i)+dx)*sinphi;// !
                }
                //! end MAN 9/15/2005
                
                xcul=xcu+(lfr.at(i)+dx)*cosphi;// ! get upstream limit on the front cavity
                ycul=ycu+(lfr.at(i)+dy)*sinphi;// !
                xcul=std::max(xcul,0.f);
                xcul=std::min(xcul,dx*(nxdx-1));
                ycul=std::max(ycul,0.f);
                ycul=std::min(ycul,dy*(nydy-1));
                cosfac=-1.f;
            }
        }
        //!mdw 7-05-2006 change form to ixxx or jxxx =nint()+1
        temp=(xcd-dx)/dx;//altered to comply with GPU
        int icd=nint(temp)+1;// ! get indicies for the downstream center line to back of the building
        temp=(ycd-dy)/dy;//altered to comply with GPU
        int jcd=nint(temp)+1;
        //!mdw 4-16-2004 added correction for ktop+3 > nz-1
        int ktp=std::min(ktop,nzdz-1);//didn't alter here as ktop is already aligned with GPU coordinates
        float zk=0.f;
        float zbrac=0.f;
        float zkfac=0.f;
        float xcdl=0.f;
        float ycdl=0.f;
        int icdl=0;
        int jcdl=0;
        int icu=0;
        int jcu=0;
        int icul=0;
        int jcul=0;
        float urefz=0.f;
        float ds=0.f;
        float sdown=0.f;
        float sup=0.f;
        float stin=0.f;
        float istinf=0.f;
        float st=0.f;
        int istf=0;
        int isf=0;
        int isfu=0;
        float utotp=0.f;
        float utotm=0.f;
        float cosu=0.f;
        float sinv=0.f;
        int isini=0;
        float cosl=0.f;
        float sinl=0.f;
        float delutz=0.f;
        float upvpg=0.f;
        float upwpg=0.f;
        float upsqg=0.f;
        float vpsqg=0.f;
        float vpwpg=0.f;
        float wpsqg=0.f;
        float duy=0.f;
        
        for(int k=ktp;k>=0;k--){//do k=ktp,2,-1  ! Account for wake difference in the cavity
            
            zk=zm.at(k);
            if(zi.at(k)<.99*h){
                zbrac=pow( (1.f-zi.at(k)/h) , 1.5f);
            }
            else{
                zbrac=pow( (1.f-.99f),1.5f);
            }
            //zbrac=pow( (1.f-zi.at(k)/h) , 1.5f);
            //!mdw 4-16-2004 added correction for ktop+3 > nz-1
            int idupc=k*nxdx*nydy + jupc*nxdx +iupc;
            int idupcktop=(ktop+3)*nxdx*nydy + jupc*nxdx +iupc;
            int idupcnzm1=(nzdz-1)*nxdx*nydy + jupc*nxdx +iupc;
            if(ktop+3<=nzdz-1){
                zcorf.at(k)=sqrt(wind_vel[idupc].u*wind_vel[idupc].u + wind_vel[idupc].v*wind_vel[idupc].v + wind_vel[idupc].w*wind_vel[idupc].w)/
                    sqrt(wind_vel[idupcktop].u*wind_vel[idupcktop].u + wind_vel[idupcktop].v*wind_vel[idupcktop].v
                         + wind_vel[idupcktop].w*wind_vel[idupcktop].w);
            }
            else{
                zcorf.at(k)=sqrt(wind_vel[idupc].u*wind_vel[idupc].u + wind_vel[idupc].v*wind_vel[idupc].v + wind_vel[idupc].w*wind_vel[idupc].w)/
                    sqrt(wind_vel[idupcnzm1].u*wind_vel[idupcnzm1].u + wind_vel[idupcnzm1].v*wind_vel[idupcnzm1].v
                         + wind_vel[idupcnzm1].w*wind_vel[idupcnzm1].w);
            }
                        
            
            //! mdw 4-16-2004 added proper treatment of zfo
            if(zk<ht[i]+zfo[i]){
                zkfac=sqrt(1.-pow((zk/(ht[i]+zfo[i])),2));
            }
            else{
                if(k==ktp){
                    zkfac=1.;
                }
                else{
                    zkfac=0.f;
                }
                
            }
            //! mdw 7-05-2006 changed from .5 dx or .5 dy to dx & dy to be consistent with nint
            xcdl=xcd+(3.*lr.at(i)+dx)*zkfac*cosphi;// ! calculate the x,y limit of the wake as a function of height
            ycdl=ycd+(3.*lr.at(i)+dy)*zkfac*sinphi;// !
            xcdl=std::min(xcdl,dx*(nxdx));
            ycdl=std::min(ycdl,dy*(nydy));
            xcdl=std::max(xcdl,0.f);
            ycdl=std::max(ycdl,0.f);
            
            temp=(xcdl-dx)/dx;//altered for GPU, substracted dx, same below
            icdl=nint(temp)+1;// ! Calculate the indicies for i,j according to xcdl,ycdl
            temp=(ycdl-dy)/dy;
            jcdl=nint(temp)+1;
            temp=(xcu-dx)/dx;
            icu=nint(temp)+1;//   ! indicies for the upstream cavity (building)
            temp=(ycu-dy)/dy;
            jcu=nint(temp)+1;//   !
            temp=(xcul-dx)/dx;
            icul=nint(temp)+1;// ! (furthest upstream)
            temp=(ycul-dy)/dy;
            jcul=nint(temp)+1;// !!!
            //!mdw 4-16-2004 added correction for ktop+3 > nz-1
            int idktop=(ktop+3)*nxdx*nydy + jcb.at(i)*nxdx +icb.at(i);
            if(ktop+3<=nzdz-1){
                //! calculating the reference wind un-disturbed by the building
                urefz=sqrt(wind_vel[idktop].u*wind_vel[idktop].u + wind_vel[idktop].v*wind_vel[idktop].v + wind_vel[idktop].w*wind_vel[idktop].w);
            }
            else{
                urefz=sqrt(pow(ualoft,2.f)+pow(valoft,2.f));
            }
            ds=0.7*std::min(dx,dy);// ! pick a step that is small enough to not skip grid cells
            sdown=sqrt((xcdl-xcd)*(xcdl-xcd)+(ycdl-ycd)*(ycdl-ycd))+2.*ds;// ! calculate the limits for the distance measured along the centerline (rear)
            sup=sqrt((xcul-xcu)*(xcul-xcu)+(ycul-ycu)*(ycul-ycu))+2.*ds;//   ! same for the front eddy
            stin=.5*leff.at(i);//
            temp=stin/ds;
            istinf=nint(temp)+1.f;
            //!mdw 7-11-2006 changed istinf to allow replacement to center of bldg
            //!mdw 5-14-2004 corrected expression for st; older versions gave errors for wide blds
            st=sqrt((xcbp3-xcb.at(i))*(xcbp3-xcb.at(i))+(ycbp3-ycb.at(i))*(ycbp3-ycb.at(i)))+1.*ds;// ! total distance to point
			temp=(st+.333*leff.at(i))/ds;
            istf=nint(temp)+1.f;//   ! (transverse direction) 
            //!mdw 6-9-2004 extended the transverse integration to st+.333*leff
            temp=sdown/ds;
            isf=nint(temp)+1;// ! setup limits of calculations (for do loops) (along cneterline down)  
            temp=sup/ds;
            isfu=nint(temp)+1;//  ! (along centerline up) 
            if(lfr.at(i) < 0.f)isfu=0;
            
            //!mdw 4-16-2004 added correction for ktop+3 > nz-1
            
            //! Select the largest reference wind of the plus or minus side of the building
            int id1=k*nxdx*nydy + jcbp3*nxdx +icbp3;
            int id2=k*nxdx*nydy + jcbm3*nxdx +icbm3;

            utotp=sqrt(wind_vel[id1].u*wind_vel[id1].u + wind_vel[id1].v*wind_vel[id1].v + wind_vel[id1].w*wind_vel[id1].w);
            utotm=sqrt(wind_vel[id2].u*wind_vel[id2].u + wind_vel[id2].v*wind_vel[id2].v + wind_vel[id2].w*wind_vel[id2].w);
            int ik=k*nxdx*nydy + i;
            int idp=k*nxdx*nydy + jcbp3*nxdx +icbp3;
            int idm=k*nxdx*nydy + jcbm3*nxdx +icbm3;
            if(utotp>=utotm){
                uref.at(ik)=utotp+.000001;
                urefu.at(ik)=uref.at(ik)*cos(phib.at(i));
                urefv.at(ik)=uref.at(ik)*sin(phib.at(i));
                urefw.at(ik)=wind_vel[idp].w;
            }
            else{
                uref.at(ik)=utotm+.000001;
                urefu.at(ik)=uref.at(ik)*cos(phib.at(i));
                urefv.at(ik)=uref.at(ik)*sin(phib.at(i));
                urefw.at(ik)=wind_vel[idm].w;
            }
            //!!!!!!!
            cosu=(urefu.at(ik)+.000001)/uref.at(ik);
            sinv=urefv.at(ik)/uref.at(ik);
            //! downstream wake  along axis do loop for delta u
            isini=1;
            float xcell=0.f;
            float ycell=0.f;
            int icel=0;
            int jcel=0;
            float utot=0.f;
            for(int is=1;is<=isf;is++){//   do is=1,isf 
                xcell=xcd+ds*(is-1)*cosphi;
                ycell=ycd+ds*(is-1)*sinphi;
                temp=(xcell-dx)/dx;//substracted dx for GPU
                icel=nint(temp)+1;
                temp=(ycell-dy)/dy;//substracted dy for GPU
                jcel=nint(temp)+1;
                icel=std::min(nxdx-1,icel);
                icel=std::max(1,icel);//altered for GPU (2 to 1)
                jcel=std::min(nydy-1,jcel);
                jcel=std::max(1,jcel);//altered for GPU (2 to 1)
                int id=k*nxdx*nydy + jcel*nxdx +icel;
                if(cellQuic[id].c == 0 && is==1){
                    isini=2;
                }
                utot=sqrt(wind_vel[id].u*wind_vel[id].u + wind_vel[id].v*wind_vel[id].v + wind_vel[id].w*wind_vel[id].w);
                //!mdw 4-16-2004 added correction for ktop+3 > nz-1
                int iceljcel=jcel*nxdx +icel;
                int idcel=ktop*nxdx*nydy + jcel*nxdx +icel;
                
                if(k==ktp){
                    if(ktop<=nzdz-1){
                        utotktp.at(iceljcel)=utot;
                        uktop.at(iceljcel)=wind_vel[idcel].u;
                        vktop.at(iceljcel)=wind_vel[idcel].v;
                        wktop.at(iceljcel)=wind_vel[idcel].w;
                    }
                    else{
                        utotktp.at(iceljcel)=sqrt(ualoft*ualoft+valoft*valoft); //check:: compare with QP, may be a bug in QP
                        uktop.at(iceljcel)=ualoft;
                        vktop.at(iceljcel)=valoft;
                        wktop.at(iceljcel)=0.f;
                    }
                }
                //! this sets reference for vertical transfer
                utot=utot+.000001;
                int idcelk=k*nxdx*nydy +jcel*nxdx +icel;
                int ik=k*nxdx*nydy+i;
                cosl=wind_vel[idcelk].u/utot;
                sinl=wind_vel[idcelk].v/utot;
                if(cellQuic[idcelk].c > 0){
                    delutz=sqrt( pow( (wind_vel[idcelk].u-zcorf.at(k)*uktop.at(iceljcel)),2.f)
                                 +pow( (wind_vel[idcelk].v -zcorf.at(k)*vktop.at(iceljcel)),2.f)
                                 +pow( (wind_vel[idcelk].w -zcorf.at(k)*wktop.at(iceljcel)),2.f) );
                    deluc.at(ik)=sqrt( pow( (urefu.at(ik)-wind_vel[idcelk].u),2.f)
                                       +pow( (urefv.at(ik)-wind_vel[idcelk].v),2.f)
                                       +pow( (urefw.at(ik)-wind_vel[idcelk].w),2.f));
                    //!mdw 4-16-2004 added correction for ktop+3 > nz-1
                    if(k!=ktp){
                        //! Selects the largest gradient (vert or horiz transfer)
                        //! mdw 4-16-2004 added proper treatment of zfo
                        
                        if((2.*deluc.at(ik)/weff.at(i))<(utotktp.at(iceljcel)/(ht[i]+zfo[i])) &&
                           delutz>.2*zcorf.at(k)*utotktp.at(iceljcel)){// ! vertical dominates
                            ustargz.at(idcelk)=std::max(knlc*utotktp.at(iceljcel),ustargz.at(idcelk)); 
                            if(fabs(ustargz.at(idcelk)-knlc*utotktp.at(iceljcel))<1.e-05*ustargz.at(idcelk)){//!This value dominates over prev. buildings.
                                elzg.at(idcelk)=ht[i]+zfo[i];
                                upvpg=0.f;
                                upwpg=-ctau13*zbrac*ustargz.at(idcelk)*ustargz.at(idcelk);
                                upsqg=cusq*zbrac*ustargz.at(idcelk)*ustargz.at(idcelk);
                                vpsqg=cvsq*zbrac*ustargz.at(idcelk)*ustargz.at(idcelk);
                                vpwpg=0.f;
                                wpsqg=cwsq*zbrac*ustargz.at(idcelk)*ustargz.at(idcelk);
                                ustarg.at(idcelk)=ustargz.at(idcelk);
                                rotate2d(idcelk,cosphi,sinphi,upsqg,upvpg,vpsqg,wpsqg,upwpg,vpwpg);
                            }
                        }
                        else{
                            //! We use the vertical gradient as dominant if it is sharper than the horizontal
                            duy=-deluc.at(ik)*sinl*cosu+deluc.at(ik)*sinv*cosl;
                            //! we now have the delta u between the outside of the bldg and the center of the wake
                            //! mdw 6-10-2004 removed
                            if(deluc.at(ik)>.2*uref.at(ik)){
                                ustarg.at(idcelk)=std::max(ustarg.at(idcelk),knlc*deluc.at(ik));
                                if(fabs(ustarg.at(idcelk)-knlc*deluc.at(ik))<1.e-05*ustarg.at(idcelk)){// ! if the horiz is dominant calculate sigmas
                                    upvpg=0.f;
                                    //! on axis u prime v prime is zero
                                    upwpg=0.f;
                                    //! for eddy transport in uv we dont consider uw
                                    upsqg=cusq*zbrac*ustarg.at(idcelk)*ustarg.at(idcelk);
                                    wpsqg=cvsq*zbrac*ustarg.at(idcelk)*ustarg.at(idcelk);
                                    vpwpg=0.f;
                                    elzg.at(idcelk)=0.5*weff.at(i);
                                    vpsqg=cwsq*zbrac*ustarg.at(idcelk)*ustarg.at(idcelk);
                                    rotate2d(idcelk,cosphi,sinphi,upsqg,upvpg,vpsqg,wpsqg,upwpg,vpwpg);
                                }
                            }
                        }
                    }
                }
                else{
                    deluc.at(ik)=0.;
                    delutz=0.;
                }
                //! transverse do loop in downstream wake
                
                for(int ist=2;ist<=istf;ist++){//do ist=2,istf 
                    //! first direction in the transverse of the wake
                    xcelt=xcell+ds*(ist-1.f)*cosphit;
                    ycelt=ycell+ds*(ist-1.f)*sinphit;
                    temp=(xcelt-dx)/dx;
                    icelt=nint(temp)+1;
                    temp=(ycelt-dy)/dy;
                    jcelt=nint(temp)+1;
                    if(fabs(xcelt-xcell)<.5f*ds)icelt=icel;
                    if(fabs(ycelt-ycell)<.5f*ds)jcelt=jcel;
                    icelt=std::min(nxdx-1,icelt);
                    icelt=std::max(1,icelt);
                    jcelt=std::min(nydy-1,jcelt);
                    jcelt=std::max(1,jcelt);
                    int iceltjcelt=jcelt*nxdx + icelt;
                    int idceltk= k*nxdx*nydy + jcelt*nxdx +icelt;
                    if(cellQuic[idceltk].c > 0){
                        utott=sqrt(wind_vel[idceltk].u*wind_vel[idceltk].u + wind_vel[idceltk].v*wind_vel[idceltk].v
                                   + wind_vel[idceltk].w*wind_vel[idceltk].w);
                        utott=utott+.000001f;
                        //!mdw 4-16-2004 added correction for ktop+3 > nz-1

                        int idceltktop=ktop*nxdx*nydy + jcelt*nxdx +icelt;
                        if(k==ktp){
                            if(ktop<nzdz-1){
                                utotktp.at(iceltjcelt)=utott;
                                uktop.at(iceltjcelt)=wind_vel[idceltktop].u;
                                vktop.at(iceltjcelt)=wind_vel[idceltktop].v;
                                wktop.at(iceltjcelt)=wind_vel[idceltktop].w;
                            }
                            else{
                                utotktp.at(iceltjcelt)=sqrt(ualoft*ualoft+valoft*valoft);
                                uktop.at(iceltjcelt)=ualoft;
                                vktop.at(iceltjcelt)=valoft;  
                                wktop.at(iceltjcelt)=0.;
                            }
                        }
                        int ik=k*nxdx*nydy +i;
                        delut=sqrt(pow( (urefu.at(ik)-wind_vel[idceltk].u),2.f)+
                                   pow( (urefv.at(ik)-wind_vel[idceltk].v),2.f)+
                                   pow( (urefw.at(ik)-wind_vel[idceltk].w),2.f));
                        delutz=sqrt(pow( (wind_vel[idceltk].u-zcorf.at(k)*uktop.at(iceltjcelt)),2.f)
                                    +pow( (wind_vel[idceltk].v-zcorf.at(k)*vktop.at(iceltjcelt)),2.f)
                                    +pow( (wind_vel[idceltk].w-zcorf.at(k)*wktop.at(iceltjcelt)),2.f));
                        //!mdw 4-16-2004 added correction for ktop+3 > nz-1
                        if(k!=ktp){
                            //! mdw 4-16-2004 added proper treatment of zfo
                            //! mdw 6-10-2004 changed to make check on centerline rather than local value
                            int ik=k*nxdx*nydy +i;
                            if((2.*deluc.at(ik)/weff.at(i))<(utotktp.at(iceltjcelt)/(ht[i]+zfo[i])) 
                               && delutz>.2*zcorf.at(k)*utotktp.at(iceltjcelt)){
                                if(ustargz.at(idceltk)<knlc*utotktp.at(iceltjcelt)){
                                    ustargz.at(idceltk)=knlc*utotktp.at(iceltjcelt);
                                    elzg.at(idceltk)=ht[i]+zfo[i];
                                    upvpg=0.;
                                    upwpg=-ctau13*zbrac*ustargz.at(idceltk)*ustargz.at(idceltk);
                                    upsqg=cusq*zbrac*ustargz.at(idceltk)*ustargz.at(idceltk);
                                    vpsqg=cvsq*zbrac*ustargz.at(idceltk)*ustargz.at(idceltk);
                                    vpwpg=0.;
                                    wpsqg=cwsq*zbrac*ustargz.at(idceltk)*ustargz.at(idceltk);
                                    ustarg.at(idceltk)=ustargz.at(idceltk);
                                    rotate2d(idceltk,cosphi,sinphi,upsqg,upvpg,vpsqg,wpsqg,upwpg,vpwpg);
                                }
                            }
                            else{
                                // We use the vertical gradient as dominant if it is sharper than the horizontal
                                cosl=wind_vel[idceltk].u/utott;
                                sinl=wind_vel[idceltk].v/utott;
                                duy=-deluc.at(ik)*sinl*cosu+deluc.at(ik)*sinv*cosl;
                                // mdw 6-10-2004 changed check from delut (local value) to deluc.at(ik); centerline
                                if(delut>.2*uref.at(ik)){
                                    if(ustarg.at(idceltk)<knlc*deluc.at(ik)){
                                        ustarg.at(idceltk)=knlc*deluc.at(ik);
                                        upvpg=-((ist-1.f)/(istf-1.f))*ustarg.at(idceltk)*ustarg.at(idceltk);
                                        upwpg=0.;
                                        // for eddy transport in uv we dont consider uw
                                        upsqg=cusq*zbrac*ustarg.at(idceltk)*ustarg.at(idceltk);
                                        wpsqg=cvsq*zbrac*ustarg.at(idceltk)*ustarg.at(idceltk);
                                        vpwpg=0.;
                                        vpsqg=cwsq*zbrac*ustarg.at(idceltk)*ustarg.at(idceltk);
                                        elzg.at(idceltk)=.5*weff.at(i);
                                        rotate2d(idceltk,cosphi,sinphi,upsqg,upvpg,vpsqg,wpsqg,upwpg,vpwpg);
                                    }
                                }
                            }
                            if(is==isini){
                                for(int isin=isini+1;isin<=istinf;isin++){//do isin=isini+1,istinf
                                    xceln=xcelt-ds*(isin-1)*cosphi;
                                    yceln=ycelt-ds*(isin-1)*sinphi;
                                    temp=(xceln-dx)/dx;
                                    iceln=nint(temp)+1;
                                    temp=(yceln-dy)/dy;
                                    jceln=nint(temp)+1;
                                    iceln=std::min(nxdx-1,iceln);
                                    iceln=std::max(1,iceln);
                                    jceln=std::min(nydy-1,jceln);
                                    jceln=std::max(1,jceln);
                                    // mdw 3/22/2004PM added if statement to avoid replacing non-zero ustarg stuff
                                    // with zero values
                                    int idcelnk= k*nxdx*nydy + jceln*nxdx +iceln;
                                    if(ustarg.at(idceltk)>ustarg.at(idcelnk)){
                                        ustarg.at(idcelnk)=ustarg.at(idceltk);
                                        elzg.at(idcelnk)=elzg.at(idceltk);
                                        ufsqgi.at(idcelnk)=ufsqgi.at(idceltk);
                                        vfsqgi.at(idcelnk)=vfsqgi.at(idceltk);
                                        wfsqgi.at(idcelnk)=wfsqgi.at(idceltk);
                                        ufvfgi.at(idcelnk)=ufvfgi.at(idceltk);
                                        ufwfgi.at(idcelnk)=ufwfgi.at(idceltk);
                                        vfwfgi.at(idcelnk)=vfwfgi.at(idceltk);
                                    }
                                    // mdw 3/22/2004PM new endif for new if then
                                }//enddo
                            }
                        }
                    }
                    // opposite direction in the transverse of the wake
                    xcelt=xcell-ds*(ist-1.f)*cosphit;
                    ycelt=ycell-ds*(ist-1.f)*sinphit;
                    temp=(xcelt-dx)/dx;
                    icelt=nint(temp)+1; 
                    temp=(ycelt-dy)/dy;
                    jcelt=nint(temp)+1; 
                    if(fabs(xcelt-xcell)<.5*ds)icelt=icel;
                    if(fabs(ycelt-ycell)<.5*ds)jcelt=jcel;
                    icelt=std::min(nxdx-1,icelt); 
                    icelt=std::max(1,icelt);
                    jcelt=std::min(nydy-1,jcelt);
                    jcelt=std::max(1,jcelt);
                    
                    iceltjcelt=jcelt*nxdx + icelt;
                    int idceltktop=ktop*nxdx*nydy + jcelt*nxdx +icelt;
                    idceltk=k*nxdx*nydy + jcelt*nxdx +icelt;
                    if(cellQuic[idceltk].c > 0){
                        utott=sqrt(wind_vel[idceltk].u*wind_vel[idceltk].u+wind_vel[idceltk].v*wind_vel[idceltk].v
                                   +wind_vel[idceltk].w*wind_vel[idceltk].w);
                        utott=utott+.000001;
                        delut=sqrt(pow( (urefu.at(ik)-wind_vel[idceltk].u),2)
                                   + pow( (urefv.at(ik)-wind_vel[idceltk].v),2)
                                   + pow( (urefw.at(ik)-wind_vel[idceltk].w),2));
                        // mdw 4-16-2004 added correction for ktop+3 > nz-1
                        
                        if(k==ktp){
                            if(ktop<=nzdz-1){
                                utotktp.at(iceltjcelt)=utott;
                                uktop.at(iceltjcelt)=wind_vel[idceltktop].u;
                                vktop.at(iceltjcelt)=wind_vel[idceltktop].v;
                                wktop.at(iceltjcelt)=wind_vel[idceltktop].w;
                            }
                            else{
                                utotktp.at(iceltjcelt)=sqrt(ualoft*ualoft+valoft*valoft);
                                uktop.at(iceltjcelt)=ualoft;
                                vktop.at(iceltjcelt)=valoft;
                                wktop.at(iceltjcelt)=0.;
                            }
                        }
                        // mdw 4-16-2004 added correction for ktop+3 > nz-1
                        
                        if(k!=ktp){
                            delutz=sqrt(pow( (wind_vel[idceltk].u-zcorf.at(k)*uktop.at(iceltjcelt)),2)
                                        +pow( (wind_vel[idceltk].v-zcorf.at(k)*vktop.at(iceltjcelt)),2)
                                        +pow( (wind_vel[idceltk].w-zcorf.at(k)*wktop.at(iceltjcelt)),2));
                            // mdw 4-16-2004 added proper treatment of zfo
                            // mdw 6-10-2004 made check on centerline rather than local value
                            if((2.*deluc.at(ik)/weff.at(i))<(utotktp.at(iceltjcelt)/(ht[i]+zfo[i])) 
                               && delutz>.2*zcorf.at(k)*utotktp.at(iceltjcelt)){
                                if(ustargz.at(idceltk)<knlc*utotktp.at(iceltjcelt)){
                                    ustargz.at(idceltk)=knlc*utotktp.at(iceltjcelt);
                                    elzg.at(idceltk)=ht[i]+zfo[i];
                                    upvpg=0.;
                                    upwpg=-ctau13*zbrac*ustargz.at(idceltk)*ustargz.at(idceltk);
                                    upsqg=cusq*zbrac*ustargz.at(idceltk)*ustargz.at(idceltk);
                                    vpsqg=cvsq*zbrac*ustargz.at(idceltk)*ustargz.at(idceltk);
                                    vpwpg=0.;
                                    wpsqg=cwsq*zbrac*ustargz.at(idceltk)*ustargz.at(idceltk);
                                    ustarg.at(idceltk)=ustargz.at(idceltk);
                                    rotate2d(idceltk,cosphi,sinphi,upsqg,upvpg,vpsqg,wpsqg,upwpg,vpwpg);
                                }
                            }
                            else{
                                // We use the vertical gradient as dominant if it is sharper than the horizontal
                                cosl=wind_vel[idceltk].u/utott;
                                sinl=wind_vel[idceltk].v/utott;
                                duy=-deluc.at(ik)*sinl*cosu+deluc.at(ik)*sinv*cosl; 
                                // mdw 6-10-2004 made check on centerline value rather than local value
                                
                                if(delut>.2f*uref.at(ik)){
                                    if(ustarg.at(idceltk)<knlc*deluc.at(ik)){
                                        ustarg.at(idceltk)=knlc*deluc.at(ik);
                                        upvpg=((ist-1.f)/(istf-1.f))*ustarg.at(idceltk)*ustarg.at(idceltk);
                                        upwpg=0.f;
                                        // for eddy transport in uv we dont consider uw
                                        upvpg=ctau13*zbrac*((ist-1.f)/(istf-1.f))*ustarg.at(idceltk)*ustarg.at(idceltk);
                                        upwpg=0.f;
                                        // for eddy transport in uv we dont consider uw
                                        upsqg=cusq*zbrac*ustarg.at(idceltk)*ustarg.at(idceltk);
                                        wpsqg=cvsq*zbrac*ustarg.at(idceltk)*ustarg.at(idceltk);
                                        elzg.at(idceltk)=0.5f*weff.at(i);
                                        vpwpg=0.f;
                                        vpsqg=cwsq*zbrac*ustarg.at(idceltk)*ustarg.at(idceltk);
                                        rotate2d(idceltk,cosphi,sinphi,upsqg,upvpg,vpsqg,wpsqg,upwpg,vpwpg);
                                    }
                                }
                            }
                            if(is==isini){
                                for(int isin=isini+1;isin<=istinf;isin++){
                                    xceln=xcelt-ds*(isin-1)*cosphi;
                                    yceln=ycelt-ds*(isin-1)*sinphi;
                                    temp=(xceln-dx)/dx;
                                    iceln=nint(temp)+1; 
                                    temp=(yceln-dy)/dy;
                                    jceln=nint(temp)+1; 
                                    iceln=std::min(nxdx-1,iceln);
                                    iceln=std::max(1,iceln);
                                    jceln=std::min(nydy-1,jceln);
                                    jceln=std::max(1,jceln);
                                    int idcelnk=k*nxdx*nydy +jceln*nxdx +iceln;
                                    // mdw 3/22/2004pm adding new if then structure to avoid replacing non-zero
                                    // ustarg with zero ones
                                    if(ustarg.at(idceltk)>ustarg.at(idcelnk)){
                                        ustarg.at(idcelnk)=ustarg.at(idceltk);
                                        elzg.at(idcelnk)=elzg.at(idceltk); 
                                        ufsqgi.at(idcelnk)=ufsqgi.at(idceltk);
                                        vfsqgi.at(idcelnk)=vfsqgi.at(idceltk);
                                        wfsqgi.at(idcelnk)=wfsqgi.at(idceltk);
                                        ufvfgi.at(idcelnk)=ufvfgi.at(idceltk);
                                        ufwfgi.at(idcelnk)=ufwfgi.at(idceltk);
                                        vfwfgi.at(idcelnk)=vfwfgi.at(idceltk);
                                    }
                                }//enddo
                                // mdw 3/22/2004 end of new if then structure
                            }//endif
                        }//endif
                    }
                }   //lp021
            }  //lp022

            isini=1;
            for(int is=1; is<=isfu;is++){//do is=1,isfu
                // upstream front eddy along the centerline
                xcell=xcu-ds*(is-1.f)*cosphi;
                ycell=ycu-ds*(is-1.f)*sinphi;
                //mdw 7-05-2006 changed form form =nint( / ) to nint( / )+1
                temp=(xcell-dx)/dx;
                icel=nint(temp)+1; 
                temp=(ycell-dy)/dy;
                jcel=nint(temp)+1; 
                icel=std::min(nxdx-1,icel);
                icel=std::max(1,icel);
                jcel=std::min(nydy-1,jcel);
                jcel=std::max(1,jcel);
                int idcelk=k*nxdx*nydy +jcel*nxdx +icel;
                int iceljcel=jcel*nxdx +icel;
                if(cellQuic[idcelk].c == 0 && is == 1){
                    isini=2;
                }
                int idcelktop=ktop*nxdx*nydy + jcel*nxdx +icel;
                idcelk=k*nxdx*nydy + jcel*nxdx +icel;
                utot=sqrt(wind_vel[idcelk].u*wind_vel[idcelk].u+wind_vel[idcelk].v*wind_vel[idcelk].v+wind_vel[idcelk].w*wind_vel[idcelk].w);
                // mdw 1-22-2004 new lines in support of bldg infiltration
                if((k==kmid)&&(is==1))utotcl1.at(i)=utot;
                if((k==kmid)&&(utot>utotmax.at(i)))utotmax.at(i)=utot; 
                utot=utot+.000001;
                //mdw 4-16-2004 added correction for ktop+3 > nz-1
                if(k==ktp){
                    if(ktop<=nzdz-1){
                        utotktp.at(iceljcel)=utot;
                        uktop.at(iceljcel)=wind_vel[idcelktop].u;
                        vktop.at(iceljcel)=wind_vel[idcelktop].v;
                        wktop.at(iceljcel)=wind_vel[idcelktop].w;
                    }
                    else{
                        utotktp.at(iceljcel)=sqrt(ualoft*ualoft+valoft*valoft);
                        uktop.at(iceljcel)=ualoft;
                        vktop.at(iceljcel)=valoft;
                        wktop.at(iceljcel)=0.;
                    }
                }
                deluc.at(ik)=sqrt(pow( (urefu.at(ik)-wind_vel[idcelk].u),2)+
                                  pow( (urefv.at(ik)-wind_vel[idcelk].v),2)+
                                  pow( (urefw.at(ik)-wind_vel[idcelk].w),2));
                //mdw 4-16-2004 added correction for ktop+3 > nz-1
                if(k!=ktp){
                    delutz=sqrt(pow( (wind_vel[idcelk].u-zcorf.at(k)*uktop.at(iceljcel)),2)+
                                pow( (wind_vel[idcelk].v-zcorf.at(k)*vktop.at(iceljcel)),2)+
                                pow( (wind_vel[idcelk].w-zcorf.at(k)*wktop.at(iceljcel)),2));
                    deluc.at(ik)=sqrt(pow( (urefu.at(ik)-wind_vel[idcelk].u),2)+
                                      pow( (urefv.at(ik)-wind_vel[idcelk].v),2)+
                                      pow( (urefw.at(ik)-wind_vel[idcelk].w),2));
                    // Selects the largest gradient (vert or horiz transfer)
                    // mdw 4-16-2004 added proper treatment of zfo
                    if((2.*deluc.at(ik)/weff.at(i))<(utotktp.at(iceljcel)/(ht[i]+zfo[i])) 
                       && delutz>.2*zcorf.at(k)*utotktp.at(iceljcel)){ // vertical dominates
                        if(ustargz.at(idcelk)<knlc*utotktp.at(iceljcel)){ // This value dominates over prev. buildings.
                            ustargz.at(idcelk)=knlc*utotktp.at(iceljcel);
                            elzg.at(idcelk)=ht[i]+zfo[i];
                            upvpg=0.;
                            upwpg=-ctau13*zbrac*ustargz.at(idcelk)*ustargz.at(idcelk);
                            upsqg=cusq*zbrac*ustargz.at(idcelk)*ustargz.at(idcelk);
                            vpsqg=cvsq*zbrac*ustargz.at(idcelk)*ustargz.at(idcelk);
                            vpwpg=0.;
                            wpsqg=cwsq*zbrac*ustargz.at(idcelk)*ustargz.at(idcelk); // align sigmas with the overall mean wind
                            ustarg.at(idcelk)=ustargz.at(idcelk);
                            rotate2d(idcelk,cosphi,sinphi,upsqg,upvpg,vpsqg,wpsqg,upwpg,vpwpg);
                        }
                    }
                    else{
                        // We use the vertical gradient as dominant if it is sharper than the horizontal
                        cosl=wind_vel[idcelk].u/utot;
                        sinl=wind_vel[idcelk].v/utot;
                        duy=-deluc.at(ik)*sinl*cosu+deluc.at(ik)*sinv*cosl;
                        // we now have the delta u between the outside of the bldg and the center of the wake
                        if(deluc.at(ik)>.2*uref.at(ik)){
                            if(ustarg.at(idcelk)<knlc*deluc.at(ik)){
                                ustarg.at(idcelk)=knlc*deluc.at(ik);
                                upvpg=0.;
                                // on axis u prime v prime is zero
                                upwpg=0.;
                                // for eddy transport in uv we dont consider uw
                                upsqg=cusq*zbrac*ustarg.at(idcelk)*ustarg.at(idcelk);
                                wpsqg=cvsq*zbrac*ustarg.at(idcelk)*ustarg.at(idcelk);
                                vpwpg=0.;
                                elzg.at(idcelk)=0.5*weff.at(i);
                                vpsqg=cwsq*zbrac*ustarg.at(idcelk)*ustarg.at(idcelk);
                                rotate2d(idcelk,cosphi,sinphi,upsqg,upvpg,vpsqg,wpsqg,upwpg,vpwpg);
                            }
                        }
                    }
                }
                for(int ist=2;ist<=istf;ist++){//do ist=2,istf
                    // first direction in the transverse of the front eddy
                    xcelt=xcell+ds*(ist-1.f)*cosphit;
                    ycelt=ycell+ds*(ist-1.f)*sinphit;
                    //mdw 7-05-2006 changed form from nint( / ) to nint( / )+1
                    temp=(xcelt-dx)/dx;
                    icelt=nint(temp)+1;
                    temp=(ycelt-dy)/dy;
                    jcelt=nint(temp)+1;
                    if(fabs(xcelt-xcell)<.5*ds)icelt=icel;
                    if(fabs(ycelt-ycell)<.5*ds)jcelt=jcel;
                    //mdw 7-11-2006 check added to use closest axis cell
                    icelt=std::min(nxdx-1,icelt);
                    icelt=std::max(1,icelt);
                    jcelt=std::min(nydy-1,jcelt);
                    jcelt=std::max(1,jcelt);
                    int idceltk=k*nxdx*nydy + jcelt*nxdx +icelt;
                    int idceltktop=ktop*nxdx*nydy + jcelt*nxdx +icelt;
                    int iceltjcelt=jcelt*nxdx + icelt;
                    utott=sqrt(wind_vel[idceltk].u*wind_vel[idceltk].u+wind_vel[idceltk].v*wind_vel[idceltk].v
                               +wind_vel[idceltk].w*wind_vel[idceltk].w);
                    utott=utott+.000001;
                    delut=sqrt(pow( (urefu.at(ik)-wind_vel[idceltk].u),2)+
                               pow( (urefv.at(ik)-wind_vel[idceltk].v),2)+
                               pow( (urefw.at(ik)-wind_vel[idceltk].w),2));
                    //mdw 4-16-2004 added correction for ktop+3 > nz-1
                    if(k==ktp){
                        if(ktop<=nzdz-1){
                            utotktp.at(iceltjcelt)=utott;
                            uktop.at(iceltjcelt)=wind_vel[idceltktop].u;
                            vktop.at(iceltjcelt)=wind_vel[idceltktop].v;
                            wktop.at(iceltjcelt)=wind_vel[idceltktop].w;
                        }
                        else{
                            utotktp.at(iceltjcelt)=sqrt(ualoft*ualoft+valoft*valoft);
                            uktop.at(iceltjcelt)=ualoft;
                            vktop.at(iceltjcelt)=valoft;
                            wktop.at(iceltjcelt)=0.;
                        }
                    }
                    //mdw 4-16-2004 added correction for ktop+3 > nz-1
                    if(k!=ktp){
                        delutz=sqrt(pow( (wind_vel[idceltk].u-zcorf.at(k)*uktop.at(iceltjcelt)),2)+
                                    pow( (wind_vel[idceltk].v-zcorf.at(k)*vktop.at(iceltjcelt)),2)+
                                    pow( (wind_vel[idceltk].w-zcorf.at(k)*wktop.at(iceltjcelt)),2));
                        // mdw 4-16-2004 added proper treatment of zfo
                        // mdw 6-10-2004 made check on centerline deluc rather than local delut
                        if((2.*deluc.at(ik)/weff.at(i))<(utotktp.at(iceltjcelt)/(ht[i]+zfo[i])) 
                           && delutz>.2*zcorf.at(k)*utotktp.at(iceltjcelt)){
                            if(ustargz.at(idceltk)<knlc*utotktp.at(iceltjcelt)){
                                ustargz.at(idceltk)=knlc*utotktp.at(iceltjcelt);
                                elzg.at(idceltk)=ht[i]+zfo[i];
                                upvpg=0.;
                                upwpg=-ctau13*zbrac*ustargz.at(idceltk)*ustargz.at(idceltk);
                                upsqg=cusq*zbrac*ustargz.at(idceltk)*ustargz.at(idceltk);
                                vpsqg=cvsq*zbrac*ustargz.at(idceltk)*ustargz.at(idceltk);
                                vpwpg=0.;
                                wpsqg=cwsq*zbrac*ustargz.at(idceltk)*ustargz.at(idceltk);
                                ustarg.at(idceltk)=ustargz.at(idceltk);
                                rotate2d(idceltk,cosphi,sinphi,upsqg,upvpg,vpsqg,wpsqg,upwpg,vpwpg);
                            }
                        }
                        else{
                            // We use the vertical gradient as dominant if it is sharper than the horizontal
                            cosl=wind_vel[idceltk].u/utott;
                            sinl=wind_vel[idceltk].v/utott;
                            duy=-deluc.at(ik)*sinl*cosu+deluc.at(ik)*sinv*cosl;
                            // mdw 6-10-2004 made check on centerline rather than local value
                            if(delut>.2*uref.at(ik)){
                                if(ustarg.at(idceltk)<knlc*deluc.at(ik)){
                                    ustarg.at(idceltk)=knlc*deluc.at(ik);
                                    // for eddy transport in uv we dont consider uw
                                    upvpg=-ctau13*zbrac*((ist-1.f)/(istf-1.f))*ustarg.at(idceltk)*ustarg.at(idceltk);
                                    upwpg=0.;
                                    // for eddy transport in uv we dont consider uw
                                    upsqg=cusq*zbrac*ustarg.at(idceltk)*ustarg.at(idceltk);
                                    wpsqg=cvsq*zbrac*ustarg.at(idceltk)*ustarg.at(idceltk);
                                    vpwpg=0.;
                                    vpsqg=cwsq*zbrac*ustarg.at(idceltk)*ustarg.at(idceltk);
                                    elzg.at(idceltk)=0.5*weff.at(i);
                                    rotate2d(idceltk,cosphi,sinphi,upsqg,upvpg,vpsqg,wpsqg,upwpg,vpwpg);
                                }
                            }
                        }
                        if(is==isini){
                            for(int isin=isini+1;isin<=istinf;isin++){//do isin=isini+1,istinf
                                xceln=xcelt+ds*(isin-1)*cosphi;
                                yceln=ycelt+ds*(isin-1)*sinphi;
                                temp=(xceln-dx)/dx;
                                iceln=nint(temp)+1; 
                                temp=(yceln-dy)/dy;
                                jceln=nint(temp)+1; 
                                iceln=std::min(nxdx-1,iceln);
                                iceln=std::max(1,iceln);
                                jceln=std::min(nydy-1,jceln);
                                jceln=std::max(1,jceln);
                                
                                int idcelnk=k*nxdx*nydy + jceln*nxdx +iceln;
                                int idcelnktop=ktop*nxdx*nydy + jceln*nxdx +iceln;
                                int icelnjceln=jceln*nxdx + iceln;
                                // mdw 3/22/2004pm added new if then structure to prevent replacing non-zero
                                // ustarg s with zero ones
                                if(ustarg.at(idceltk)>ustarg.at(idcelnk)){
                                    ustarg.at(idcelnk)=ustarg.at(idceltk);
                                    elzg.at(idcelnk)=elzg.at(idceltk);
                                    ufsqgi.at(idcelnk)=ufsqgi.at(idceltk);
                                    vfsqgi.at(idcelnk)=vfsqgi.at(idceltk);
                                    wfsqgi.at(idcelnk)=wfsqgi.at(idceltk);
                                    ufvfgi.at(idcelnk)=ufvfgi.at(idceltk);
                                    ufwfgi.at(idcelnk)=ufwfgi.at(idceltk);
                                    vfwfgi.at(idcelnk)=vfwfgi.at(idceltk);
                                }
                                // mdw 3/22/2004pm end  of new if then structure
                            }//enddo
                        }
                    }
                    // opposite direction in the transverse of the front eddy
                    xcelt=xcell-ds*(ist-1.f)*cosphit;
                    ycelt=ycell-ds*(ist-1.f)*sinphit;
                    //mdw 7-05-2006 changed form from nint( / ) to nint( / )+1
                    temp=(xcelt-dx)/dx;
                    icelt=nint(temp)+1; 
                    temp=(ycelt-dy)/dy;
                    jcelt=nint(temp)+1; 
                    
                    if(fabs(xcelt-xcell)<.5*ds)icelt=icel;
                    if(fabs(ycelt-ycell)<.5*ds)jcelt=jcel;
                    //mdw 7-11-2006 check added to use closest axis cell
                    icelt=std::min(nxdx-1,icelt);
                    icelt=std::max(1,icelt);
                    jcelt=std::min(nydy-1,jcelt);
                    jcelt=std::max(1,jcelt);
                    iceltjcelt=jcelt*nxdx + icelt;
                    idceltktop=ktop*nxdx*nydy+jcelt*nxdx +icelt;
                    idceltk=k*nxdx*nydy +jcelt*nxdx +icelt;
                    utott=sqrt(wind_vel[idceltk].u*wind_vel[idceltk].u+wind_vel[idceltk].v*wind_vel[idceltk].v
                               +wind_vel[idceltk].w*wind_vel[idceltk].w);
                    
                    utott=utott+.000001;
                    //mdw 4-16-2004 added correction for ktop+3 > nz-1
                    if(k==ktp){
                        if(ktop<=nzdz-1){
                            utotktp.at(iceltjcelt)=utott;
                            uktop.at(iceltjcelt)=wind_vel[idceltktop].u;
                            vktop.at(iceltjcelt)=wind_vel[idceltktop].v;
                            wktop.at(iceltjcelt)=wind_vel[idceltktop].w;
                        }
                        else{
                            utotktp.at(iceltjcelt)=sqrt(ualoft*ualoft+valoft*valoft);
                            uktop.at(iceltjcelt)=ualoft;
                            vktop.at(iceltjcelt)=valoft;
                            wktop.at(iceltjcelt)=0.;
                        }
                    }
                    delut=sqrt(pow( (urefu.at(ik)-wind_vel[idceltk].u),2)
                               +pow( (urefv.at(ik)-wind_vel[idceltk].v),2)
                               +pow( (urefw.at(ik)-wind_vel[idceltk].w),2));
                    //mdw 4-16-2004 added correction for ktop+3 > nz-1
                    if(k!=ktp){
                        delutz=sqrt(pow( (wind_vel[idceltk].u-zcorf.at(k)*uktop.at(iceltjcelt)),2)+
                                    pow( (wind_vel[idceltk].v-zcorf.at(k)*vktop.at(iceltjcelt)),2)+
                                    pow( (wind_vel[idceltk].w-zcorf.at(k)*wktop.at(iceltjcelt)),2));
                        // mdw 4-16-2004 added proper treatment of zfo
                        // mdw 6-10-2004 made check on centerline rather than local value
                        if((2.*deluc.at(ik)/weff.at(i))<(utotktp.at(iceltjcelt)/(ht[i]+zfo[i])) 
                           &&delutz>.2*zcorf.at(k)*utotktp.at(iceltjcelt)){
                            
                            if(ustargz.at(idceltk)<knlc*utotktp.at(iceltjcelt)){
                                ustargz.at(idceltk)=knlc*utotktp.at(iceltjcelt);
                                elzg.at(idceltk)=ht[i]+zfo[i];
                                upvpg=0.;
                                upwpg=-ctau13*zbrac*ustargz.at(idceltk)*ustargz.at(idceltk);
                                upsqg=cusq*zbrac*ustargz.at(idceltk)*ustargz.at(idceltk);
                                vpsqg=cvsq*zbrac*ustargz.at(idceltk)*ustargz.at(idceltk);
                                vpwpg=0.;
                                wpsqg=cwsq*zbrac*ustargz.at(idceltk)*ustargz.at(idceltk);
                                ustarg.at(idceltk)=ustargz.at(idceltk);
                                rotate2d(idceltk,cosphi,sinphi,upsqg,upvpg,vpsqg,wpsqg,upwpg,vpwpg);
                            }
                        }
                        else{
                            // We use the vertical gradient as dominant if it is sharper than the horizontal
                            cosl=wind_vel[idceltk].u/utott;
                            sinl=wind_vel[idceltk].v/utott;
                            duy=-deluc.at(ik)*sinl*cosu+deluc.at(ik)*sinv*cosl;
                            // mdw 6-10-2004 made check on centerline rather than local (delut) value
                            if(delut>.2*uref.at(k)){
                                if(ustarg.at(idceltk)<knlc*deluc.at(ik)&&ustargz.at(idceltk)<knlc*deluc.at(ik)){
                                    ustarg.at(idceltk)=knlc*deluc.at(ik);
                                    
                                    // for eddy transport in uv we dont consider uw
                                    float tau13=0.;
                                    upvpg=-tau13*zbrac*((ist-1.f)/(istf-1.f))*ustarg.at(idceltk)*ustarg.at(idceltk);//check:: might be bug in QP
                                    upwpg=0.;
                                    // for eddy transport in uv we dont consider uw
                                    upsqg=cusq*zbrac*ustarg.at(idceltk)*ustarg.at(idceltk);
                                    wpsqg=cvsq*zbrac*ustarg.at(idceltk)*ustarg.at(idceltk);
                                    vpwpg=0.;
                                    elzg.at(idceltk)=0.5*weff.at(i);
                                    vpsqg=cwsq*zbrac*ustarg.at(idceltk)*ustarg.at(idceltk);
                                    rotate2d(idceltk,cosphi,sinphi,upsqg,upvpg,vpsqg,wpsqg,upwpg,vpwpg);
                                }
                            }
                        }
                        if(is==isini){
                            for(int isin=isini+1;isin<=istinf;isin++){//do isin=isini+1,istinf{
                                xceln=xcelt+ds*(isin-1)*cosphi;
                                yceln=ycelt+ds*(isin-1)*sinphi;
                                temp=(xceln-dx)/dx;
                                iceln=nint(temp)+1; 
                                temp=(yceln-dy)/dy;
                                jceln=nint(temp)+1; 
                                iceln=std::min(nxdx-1,iceln);
                                iceln=std::max(1,iceln);
                                jceln=std::min(nydy-1,jceln);
                                jceln=std::max(1,jceln);
                                int idcelnk=k*nxdx*nydy + jceln*nxdx +iceln;
                                int idcelnktop=ktop*nxdx*nydy + jceln*nxdx +iceln;
                                int icelnjceln=jceln*nxdx + iceln;
                                // mdw 3/22/2004pm added new if then structure to prevent replacing non-zero
                                // ustargs with zero ones
                                if(ustarg.at(idceltk)>ustarg.at(idcelnk)){
                                    ustarg.at(idcelnk)=ustarg.at(idceltk);
                                    
                                    elzg.at(idcelnk)=elzg.at(idceltk);
                                    ufsqgi.at(idcelnk)=ufsqgi.at(idceltk);
                                    vfsqgi.at(idcelnk)=vfsqgi.at(idceltk);
                                    wfsqgi.at(idcelnk)=wfsqgi.at(idceltk);
                                    ufvfgi.at(idcelnk)=ufvfgi.at(idceltk);
                                    ufwfgi.at(idcelnk)=ufwfgi.at(idceltk);
                                    vfwfgi.at(idcelnk)=vfwfgi.at(idceltk);
                                }
                                // mdw 3/22/2004pm end of new if then structure
                            }
                        }
                    }

                }//   lp023
            }//   lp024
        }//   lp025
        float xpent1, ypent1;
        int npentx,npenty,ipent1,ipent2,jpent1,jpent2,ibuild;
        
        ibuild=1;
        switch(bldtype.at(i)){
        case(3):
            xpent1=xfo[i]-wti[i]*.2f-dx;
            temp = ((.4f*wti[i])/dx); 
            npentx=nint(temp)+1; 
            ypent1=yfo[i]-wti[i]*.2-dy;
            temp=((.4f*wti[i])/dy);
            npenty=nint(temp)+1; 
            temp = (xpent1-dx)/dx;
            ipent1=nint(temp)+1; 
            ipent2=ipent1+npentx;
            temp= ((ypent1-dy)/dy);
            jpent1=nint(temp)+1; 
            jpent2=jpent1+npenty;
            
            for(int icel=ipent1;icel<=ipent2;icel++){
                for(int jcel=jpent1;jcel<=jpent2;jcel++){
                    for(int k=ktp;k>=0;k--){
                        int idcelk=k*nxdx*nydy + jcel*nxdx +icel;    
                        int iceljcel=jcel*nxdx +icel;
                        
                        utot=sqrt(wind_vel[idcelk].u*wind_vel[idcelk].u+wind_vel[idcelk].v*wind_vel[idcelk].v+wind_vel[idcelk].w*wind_vel[idcelk].w);
                        utot=utot+.000001;
                        
                        //mdw 4-16-2004 added correction for ktop+3 > nz-1
                        if(k==ktp){
                            if(ktop<=nzdz){
                                int idcelktop=ktop*nxdx*nydy + jcel*nxdx +icel;
                                int iceljcel=jcel*nxdx +icel; 
                                utotktp.at(iceljcel)=utot;
                                uktop.at(iceljcel)=wind_vel[idcelktop].u;
                                vktop.at(iceljcel)=wind_vel[idcelktop].v;
                                wktop.at(iceljcel)=wind_vel[idcelktop].w;
                            }
                            else{
                                utotktp.at(iceljcel)=sqrt(ualoft*ualoft+valoft*valoft);
                                uktop.at(iceljcel)=ualoft;
                                vktop.at(iceljcel)=valoft;
                                wktop.at(iceljcel)=0.;
                            }
                        }
                        if(k!=ktp && cellQuic[idcelk].c != 0){
                            // MAN 9/14/2005 pentagon courtyard nonlocal mixing fix
                            delutz=sqrt((wind_vel[idcelk].u-uktop.at(iceljcel))*(wind_vel[idcelk].u-uktop.at(iceljcel))
                                        +(wind_vel[idcelk].v-vktop.at(iceljcel))*(wind_vel[idcelk].v-vktop.at(iceljcel))+ 
                                        (wind_vel[idcelk].w-wktop.at(iceljcel))*(wind_vel[idcelk].w-wktop.at(iceljcel)));
                            if(delutz>.2*utotktp.at(iceljcel)){ // vertical dominates
                                // end MAN 9/14/2005           
                                if(ustargz.at(idcelk)<knlc*utotktp.at(iceljcel)){ // This value dominates over prev. buildings.
                                    ustargz.at(idcelk)=knlc*utotktp.at(iceljcel);
                                    elzg.at(idcelk)=ht[i]+zfo[i];
                                    upvpg=0.;
                                    upwpg=-ustargz.at(idcelk)*ustargz.at(idcelk);
                                    upsqg=6.25*ustargz.at(idcelk)*ustargz.at(idcelk);
                                    vpsqg=(4./6.25)*upsqg;
                                    vpwpg=0.;
                                    wpsqg=1.69*ustargz.at(idcelk)*ustargz.at(idcelk); // align sigmas with the overall mean wind
                                    upvpg=0.;
                                    upwpg=-ctau13*zbrac*ustargz.at(idcelk)*ustargz.at(idcelk);
                                    upsqg=cusq*zbrac*ustargz.at(idcelk)*ustargz.at(idcelk);
                                    vpsqg=cvsq*zbrac*ustargz.at(idcelk)*ustargz.at(idcelk);
                                    vpwpg=0.;
                                    wpsqg=cwsq*zbrac*ustargz.at(idcelk)*ustargz.at(idcelk); // align sigmas with the overall mean wind
                                    ustarg.at(idcelk)=ustargz.at(idcelk);
                                    rotate2d(idcelk,cosphi,sinphi,upsqg,upvpg,vpsqg,wpsqg,upwpg,vpwpg);
                                }
                            }
                        }
                    }
                    check1=check1+1;
                }
                check2=check2+1;
            }
            check3=check3+1;
            break;
        case(4):
        case(5):
            
            float x0=xfo[i]+0.5*lti[i]*cos(gamma.at(i)*pi/180.);
            float y0=yfo[i]+0.5*lti[i]*sin(gamma.at(i)*pi/180.);
            
            float x1=xfo[i]+0.5*wti[ibuild]*sin(gamma.at(i)*pi/180.);
            float y1=yfo[i]-0.5*wti[ibuild]*cos(gamma.at(i)*pi/180.);
            float x2=x1+lti[i]*cos(gamma.at(i)*pi/180.);
            float y2=y1+lti[i]*sin(gamma.at(i)*pi/180.);
            float x4=xfo[i]-0.5*wti[i]*sin(gamma.at(i)*pi/180.);
            float y4=yfo[i]+0.5*wti[i]*cos(gamma.at(i)*pi/180.);
            float x3=x4+lti[i]*cos(gamma.at(i)*pi/180.);
            float y3=y4+lti[i]*sin(gamma.at(i)*pi/180.);
            float temp1 = std::min(x1,x2);
            float temp2 = std::min(temp1,x3);
            float temp3 = std::min(temp2,x4);
            int icelmin = int(temp3/dx);
            temp1 = std::max(x1,x2);
            temp2 = std::max(temp1,x3);
            temp3 = std::max(temp2,x4);
            int icelmax = int(temp3/dx);

            temp1 = std::min(y1,y2);
            temp2 = std::min(temp1,y3);
            temp3 = std::min(temp2,y4);
            int jcelmin = int(temp3/dy);
            temp1 = std::max(y1,y2);
            temp2 = std::max(temp1,y3);
            temp3 = std::max(temp2,y4);
            int jcelmax = int(temp3/dy);
            
            for(int icel=icelmin;icel<=icelmax+1;icel++){
                for(int jcel=jcelmin; jcel<=jcelmax+1;jcel++){
                    float xc=(((icel)-0.5)*dx-x0)*cos(gamma.at(i)*pi/180.)+
                        (((jcel)-0.5)*dy-y0)*sin(gamma.at(i)*pi/180.);
                    float yc=-(((icel)-0.5)*dx-x0)*sin(gamma.at(i)*pi/180.)+
                        (((jcel)-0.5)*dy-y0)*cos(gamma.at(i)*pi/180.);
					int kk=0;
                    for(int k=1;k<=ktp;k++){
                        kk=k;
                        if(zfo[i]<z.at(k))break;
                    }
                    int kzfo=kk;
                    
                    for(k=ktp;k>=kzfo;k--){ 
                        dz=dz_array.at(k);
                        int incourt=0;
                        int idcelk=k*nxdx*nydy + jcel*nxdx +icel;
                        int idcelktop=ktop*nxdx*nydy + jcel*nxdx +icel;
                        int iceljcel=jcel*nxdx +icel;
                        if(cellQuic[idcelk].c != 0){
                            utot=sqrt(wind_vel[idcelk].u*wind_vel[idcelk].u+wind_vel[idcelk].v*wind_vel[idcelk].v+wind_vel[idcelk].w*wind_vel[idcelk].w)+.000001f;
                            if(bldtype.at(i) == 4){
                                if(xc > -0.5*lti[i] && xc < 0.5*lti[i] && 
                                   yc > -0.5*wti[i] && yc < 0.5*wti[i]){
                                    incourt=1;
                                }
                            }
                            else{
                                float rc=sqrt((xc*xc)+(yc*yc));
                                float tc=atan2(yc,xc);
                                if(rc < 0.25*lti[i]*wti[i]/
                                   sqrt((pow( (0.5f*lti[i]*sin(tc)),2.f))+(pow( (0.5f*wti[i]*cos(tc)),2.f)))){
                                    incourt=1;
                                }
                            }
                        }
                        else{
                            continue;
                        }
                        //mdw 4-16-2004 added correction for ktop+3 > nz-1
                        if(incourt == 1){
                            if(k==ktp){
                                if(ktop<=nz-1){
                                    utotktp.at(iceljcel)=utot;
                                    uktop.at(iceljcel)=wind_vel[idcelktop].u;
                                    vktop.at(iceljcel)=wind_vel[idcelktop].v;
                                    wktop.at(iceljcel)=wind_vel[idcelktop].w;
                                }
                                else{
                                    utotktp.at(iceljcel)=sqrt(ualoft*ualoft+valoft*valoft);
                                    uktop.at(iceljcel)=ualoft;
                                    vktop.at(iceljcel)=valoft;
                                    wktop.at(iceljcel)=0.;
                                }
                            }
                            if(k!=ktp && cellQuic[idcelk].c != 0){
                                // MAN 9/14/2005 pentagon courtyard nonlocal mixing fix
                                delutz=sqrt((wind_vel[idcelk].u-uktop.at(iceljcel))*(wind_vel[idcelk].u-uktop.at(iceljcel))
                                            +(wind_vel[idcelk].v-vktop.at(iceljcel))*(wind_vel[idcelk].v-vktop.at(iceljcel))+ 
                                            (wind_vel[idcelk].w-wktop.at(iceljcel))*(wind_vel[idcelk].w-wktop.at(iceljcel)));
                                if(delutz>.2*utotktp.at(iceljcel)){ // vertical dominates
                                    // end MAN 9/14/2005              
                                    if(ustargz.at(idcelk)<knlc*utotktp.at(iceljcel)){ // This value dominates over prev. buildings.
                                        ustargz.at(idcelk)=knlc*utotktp.at(iceljcel);
                                        elzg.at(idcelk)=ht[i]+zfo[i];
                                        upvpg=0.;
                                        upwpg=-ustargz.at(idcelk)*ustargz.at(idcelk);
                                        upsqg=6.25*ustargz.at(idcelk)*ustargz.at(idcelk);
                                        vpsqg=(4./6.25)*upsqg;
                                        vpwpg=0.;
                                        wpsqg=1.69*ustargz.at(idcelk)*ustargz.at(idcelk); // align sigmas with the overall mean wind
                                        upvpg=0.;
                                        upwpg=-ctau13*zbrac*ustargz.at(idcelk)*ustargz.at(idcelk);
                                        upsqg=cusq*zbrac*ustargz.at(idcelk)*ustargz.at(idcelk);
                                        vpsqg=cvsq*zbrac*ustargz.at(idcelk)*ustargz.at(idcelk);
                                        vpwpg=0.;
                                        wpsqg=cwsq*zbrac*ustargz.at(idcelk)*ustargz.at(idcelk); // align sigmas with the overall mean wind
                                        ustarg.at(idcelk)=ustargz.at(idcelk);
                                        rotate2d(idcelk,cosphi,sinphi,upsqg,upvpg,vpsqg,wpsqg,upwpg,vpwpg);
                                    }
                                }
                            }
                        }
                    }
                    check1=check1+1;
                }
                check2=check2+1;
            }
            check3=check3+1;
        }
        check=check+1;
        
         
    }//for loop for buildings for non-local mixing ends
    
    // Following code is for local mixing-Balli -06/14/09
    // calculate distance to ground and walls if within 2 cells
    float zbrac=0.;
    float m_roof=6.;
    float eps=0.;
    float sigu=0.;
    float sigv=0.;
    float sigw=0.;
    float upwp=0.;
    float delym=0.;
    float delxm=0.;
    float u3psq=0.;
    float v3psq=0.;
    float w3psq=0.;
    float upvp=0.;
    float vpwp=0.;
    float ufwf=0.;
    float ufvf=0.;
    float vfwf=0.;
    float utotm=0.;
    float utotp=0.;
    float dutotdxp=0.;
    float dutotdxm=0.;
    float dutotdxa=0.;
    float dutotdyp=0.;
    float dutotdym=0.;
    float dutotdyc=0.;
    float dutotdya=0.;
    float x_b=0.;
    float y_b=0.;
    float dwallg=0.;
    float elzv=0.;
    float xloc=0.;
    float ufsq=0.;
    float vfsq=0.;
    float wfsq=0.;
    float dwall=0.;
    float ufsqb=0.;
    float wfsqb=0.;
    float vfsqb=0.;
    
    
    
    for(int j=0;j<nydy;j++){//altered for GPU, in QP it was ->do j=1,ny-1 -Balli(06/14/09)
        for(int i=0;i<nxdx;i++){//altered for GPU, in QP it was ->do i=1,nx-1 -Balli(06/14/09)
            int ij=j*nxdx+i;
            if(hgtveg.at(ij) > 0.){
                for(int  kk=1;kk<=nzdz;kk++){
                    kcantop=kk;
                    if(hgtveg.at(ij) <= z.at(kk))break;
                }
                int idcan=kcantop*nxdx*nydy+j*nxdx+i;
                int idcan1=(kcantop+1)*nxdx*nydy+j*nxdx+i;
                ucantop=.5*sqrt(wind_vel[idcan].u*wind_vel[idcan].u+wind_vel[idcan].v*wind_vel[idcan].v+wind_vel[idcan].w*wind_vel[idcan].w)+
                    .5*sqrt(wind_vel[idcan1].u*wind_vel[idcan1].u+wind_vel[idcan1].v*wind_vel[idcan1].v+wind_vel[idcan1].w*wind_vel[idcan1].w);
            }
            for(int k=0;k<nzdz;k++){ //altered for GPU, in QP it was ->do k=2,nz-1 -Balli(06/14/09)
                sigu=0.;
                sigv=0.;
                sigw=0.;
                
                int id=k*nxdx*nydy +j*nxdx +i;
                int idkm1=(k-1)*nxdx*nydy +j*nxdx +i;
                int kp1=(k+1)*nxdx*nydy +j*nxdx +i;
                
                int row = k / (numInRow);
                int texidx = row * width * nydy * 4 +
                    j * width * 4 +
                    k % (numInRow) * nxdx * 4 +
                    i * 4;
                
                sigui.at(id)=0.;
                sigvi.at(id)=0.;
                sigwi.at(id)=0.;


                //Balli:New stuff- vertical varying grid
                dz=dz_array.at(k);
                elcanopy=0.;
                int klim=1;
                
                if(cellQuic[id].c==8){
                    zbrac=1.;
                }
                else{
                    if(zi.at(k)<.99*h){
                        zbrac=pow( (1.f-zi.at(k)/h),1.5f);
                    }
                    else{
                        zbrac=pow( (1.f-.99f),1.5f);
                    }
                }
                dzm.at(id)=std::max(zm.at(k)-hgt.at(ij),0.f);
                if(dzm.at(id)>zi.at(k)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f)&&(zi.at(k)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f)>0.f))
                    dzm.at(id)=zi.at(k)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f);
                // MAN 9/21/2005 roof top mixing length fix
                eleff.at(id)=dzm.at(id);
                if(cellQuic[id].c==8)eleff.at(id)=elcanopy*std::min(1.,(dzm.at(id)-z0)/(.3*hgtveg.at(ij)))/kkar;
                int kdif=k-1;
                eleff.at(id)=std::max(zm.at(k)-hgt.at(ij)*pow( (hgt.at(ij)/zm.at(k)),m_roof),0.0);
                if(zi.at(k)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f)>0.f)
                    eleff.at(id)=std::min(eleff.at(id),zi.at(k)-std::max((hgtveg.at(ij)-elcanopy/kkar),0.f));
                if(cellQuic[id].c==8)eleff.at(id)=elcanopy*std::min(1.,(zi.at(k)-z0)/(.3*hgtveg.at(ij)))/kkar;
                klim=nzdz;
                // calculation of ustar in the vertical
                if(cellQuic[idkm1].c == 0 && cellQuic[id].c != 0){
                    utot=sqrt(wind_vel[id].u*wind_vel[id].u+wind_vel[id].v*wind_vel[id].v*+wind_vel[id].w*wind_vel[id].w);
                    if(cellQuic[id].c!=8){
                        if(rcl>0){
                            phim=1.+4.7*rcl*0.5*dz;
                            psim=-4.7*rcl*0.5*dz;
                        }
                        else{
                            phim=pow( (1.f-15.f*rcl*0.5f*dz),(-.25f));
                            psim=2.*log((1.+1./phim)/2.)+log((1.+1./pow(phim,2.f))/2.)-2.*atan(1./phim)+pi/2.;
                        }
                    }
                    else{
                        phim=1.;
                        psim=0.;
                        elz.at(id)=elcanopy*std::min(1.,(.5*dz-z0)*kkar/(.3*hgtveg.at(ij)));
                        ustar=elz.at(id)*utot/(.5*dz-z0);
                        dutotdzi.at(id)=ustar/elz.at(id);
                        ustarz.at(id)=ustar;
                    }
                }
                else{
                    if(cellQuic[id].c != 0){
                        utotu=sqrt(wind_vel[kp1].u*wind_vel[kp1].u+wind_vel[kp1].v*wind_vel[kp1].v+wind_vel[kp1].w*wind_vel[kp1].w);
                        utot=sqrt(wind_vel[id].u*wind_vel[id].u+wind_vel[id].v*wind_vel[id].v*+wind_vel[id].w*wind_vel[id].w);
                        if(fabs(dutotdzi.at(id))>1.e-06 && cellQuic[id].c!=8 && dzm.at(id)>2.*dz){
                            elz.at(id)=kkar*utot/fabs(dutotdzi.at(id));
                            // MAN 9/21/2005 roof top mixing length fix
                            if((kkar*eleff.at(id))<elz.at(id))
                                elz.at(id)=kkar*eleff.at(id);
                            else
                                elz.at(id)=kkar*eleff.at(id);
                            if(cellQuic[id].c==8){
                                elz.at(id)=elcanopy*std::min(1.,(zi.at(k)-z0)/(.3*hgtveg.at(ij)));
                                eleff.at(id)=elz.at(id)/kkar;
                            }
                        }
                        if(k < nz-1){
                            if((cellQuic[kp1].c!=8)&&(cellQuic[id].c==8)){
                                if(fabs(dutotdzi.at(id))>1.e-06){
                                    elz.at(id)=kkar*utot/fabs(dutotdzi.at(id));
                                    if((kkar*eleff.at(id))<elz.at(id)) elz.at(id)=kkar*eleff.at(id);
                                }
                            }
                        }
                        // We have just put in the vortex mixing length for the last cell in the canopy
                        if(cellQuic[idkm1].c!=8){
                            if(rcl>0){
                                phim=1.+4.7*rcl*eleff.at(id);
                                psim=-4.7*rcl*eleff.at(id);
                            }
                            else{
                                phim=pow( (1.f-15.f*rcl*eleff.at(id)),(-.25f));
                                psim=2.*log((1.+1./phim)/2.)+log((1.+1./pow(phim,2.f))/2.)-2.*atan(1./phim)+pi/2.;
                            }
                            ustar=kkar*eleff.at(id)*fabs(dutotdzi.at(id))/phim;
                            ustarz.at(id)=ustar;
                        }
                        else{
                            utotl=sqrt(wind_vel[idkm1].u*wind_vel[idkm1].u+wind_vel[idkm1].v*wind_vel[idkm1].v+wind_vel[idkm1].w*wind_vel[idkm1].w);
                            dutotdzi.at(id)=(utotu-utotl)/(2.*dz);
                            //corrected wrong gradient at the vegetative canopy top 12/22/2008
                            if(cellQuic[id].c!=8){
                                elz.at(id)=kkar*.5*dz;
                            }
                            else{
                                elz.at(id)=elcanopy*std::min(1.,(zi.at(k)-z0)/(.3*hgtveg.at(ij)));
                            }
                            ustar=elz.at(id)*fabs(dutotdzi.at(id));
                            ustarz.at(id)=ustar;
                        }
                        if(cellQuic[idkm1].c!=8){
                            if(rcl>0){
                                phim=1.+4.7*rcl*eleff.at(id);
                                psim=-4.7*rcl*eleff.at(id);
                            }
                            else{
                                phim=pow( (1.f-15.f*rcl*eleff.at(id)),(-.25f));
                                psim=2.*log((1.+1./phim)/2.)+log((1.+1./pow(phim,2.f))/2.)-2.*atan(1./phim)+pi/2.;
                            }
                            ustar=kkar*eleff.at(id)*fabs(dutotdzi.at(id))/phim;
                            ustarz.at(id)=ustar;
                        }
                    }
                }
                // for neutral conditions sigw is only dependent on ustar
                
                // for vertical downward distance (dzm)
                klim=nzdz-1;
                
                float dutotdn=0.f;
                float dutotds=0.f;
                
                sig[id].u = 0.0f;
                sig[id].v = 0.0f;   
                sig[id].w = 0.0f;   
                
                tau[id].t11   = 0.0f;        
                tau[id+1].t22 = 0.0f;        
                tau[id+2].t33 = 0.0f;
                tau[id+3].t13 = 0.0f;
                //Make tau's a texture so that they can be visualized as horizontal layers in the domain
                dataTau[texidx] = 0.0f;   
                dataTau[texidx+1] =0.0f;  
                dataTau[texidx+2] =0.0f;  
                dataTau[texidx+3] = 0.0f; 
                
                dataTwo[texidx]   = 0.0f;
                dataTwo[texidx+1] =  0.0f;
                dataTwo[texidx+2] =  0.0f;
                dataTwo[texidx+3] =  0.0f;
                
                dataWind[texidx]   = wind_vel[id].u;
                dataWind[texidx+1] = wind_vel[id].v;
                dataWind[texidx+2] = wind_vel[id].w;
                dataWind[texidx+3] = 0.;
                if(k>0){
                    dzp.at(id)=10.f*dz+z.at(klim)-z.at(k-1)+.5f*dz;
                }
                else{
                    dzp.at(id)=10.f*dz+z.at(klim)+.5f*dz;
                }
                for(int kk=k;kk<=klim;kk++){//do kk=k,klim
                    int idkk=kk*nxdx*nydy +j*nxdx +i;
                    int celltypeidkk=0;
                    if(idkk>=0){
                        celltypeidkk=cellQuic[idkk].c;
                    }
                    if(celltypeidkk == 0){
                        dzp.at(id)=.5*dz+(kk-k-1)*dz;
                        break;
                    }
                }
                // for distance to the left (dxm)
                int ilim=0;
                dxm.at(id)=.5*dx+i*dx+(nxdx+1)*dx;//altered for GPU [added 1 to nxdx as in QP nx is nx+1 actually] -Balli(06/14/09)
                for(int ii=i;ii>ilim;ii--){//do ii=i,ilim,-1
                    // calculation of the distance to the wall in the negative x direction
                    int idii=k*nxdx*nydy +j*nxdx +ii;
                    int celltypeidii=1;//as in x or y direction we assume fluid in all directions out of domain
                    if(idii>=0){
                        celltypeidii=cellQuic[idii].c;
                    }                        
                    if(celltypeidii == 0){
                        dxm.at(id)=.5*dx+(i-ii-1)*dx;
                        break;
                    }
                }
                // for distance to the right (dxp)
                ilim=nxdx-1;
                dxp.at(id)=.5*dx+(ilim-i)*dx+(nxdx+1)*dx;
                for(int ii=i;ii<=ilim;ii++){// ii=i,ilim
                    // calculation of the distance to the wall in the positive x direction
                    int idii=k*nxdx*nydy +j*nxdx +ii;
                    int celltypeidii=1;//as in x or y direction we assume fluid in all directions out of domain
                    if(idii>=0){
                        celltypeidii=cellQuic[idii].c;
                    }
                    
                    if(celltypeidii == 0){
                        dxp.at(id)=.5*dx+(ii-i-1)*dx;
                        break;
                    }
                }
                // for distance  from the back (dym)
                int jlim=0;
                dym.at(id)=.5*dy+j*dy+(nydy+1)*dy; //added 1 to nydy ,same reason as in x-directiom, see above-Balli(06/14/09)
                for(int jj=j;jj>jlim;jj--){//do jj=j,jlim,-1
                    // calculation of the distance to the wall in the negative y direction
                    int idjj=k*nxdx*nydy +jj*nxdx +i;
                    int celltypeidjj=1;//as in x or y direction we assume fluid in all directions out of domain
                    if(idjj>=0){
                        celltypeidjj=cellQuic[idjj].c;
                    }
                    if(celltypeidjj == 0){
                        dym.at(id)=.5*dy+(j-jj-1)*dy;
                        break;
                    }
                }
                // for distance to the front  (dyp)
                jlim=nydy-1;
                dyp.at(id)=.5*dy+(jlim-j)*dy+(nydy+1)*dy; //added 1 to nydy ,same reason as in x-directiom, see above-Balli(06/14/09)
                for(int jj=j;jj<=jlim;jj++){//do jj=j,jlim
                    // calculation of the distance to the wall in the positive x direction
                    int idjj=k*nxdx*nydy +jj*nxdx +i;
                    int celltypeidjj=1;//as in x or y direction we assume fluid in all directions out of domain
                    if(idjj>=0){
                        celltypeidjj=cellQuic[idjj].c;
                    }
                    if(celltypeidjj == 0){
                        dyp.at(id)=.5*dy+(jj-j-1)*dy;
                        break;
                    }
                }
                // we need to calculate the largest change in utot
                if(cellQuic[id].c == 0){
                    eps=0.;
                    sigu=0.;
                    sigv=0.;
                    sigw=0.;
                    upwp=0.;
                    elz.at(id)=0.;
                    eleff.at(id)=0.;
                    sig[id].u = 0.0f;
                    sig[id].v = 0.0f;   
                    sig[id].w = 0.0f;
                }
                if(cellQuic[id].c != 0){//for all fuid cells
                    // first we set up parameters for cells near boundary
                    if(j<1||j>=ny-1||i<1||i>=nx-1){//boundary cells
                        // calculation of near-boundary values of u*y, ly, dely, and the
                        // gradients of speed in the x and y directions
                        delym=(nydy+1)*dy;
                        utot=sqrt(wind_vel[id].u*wind_vel[id].u+wind_vel[id].v*wind_vel[id].v+wind_vel[id].w*wind_vel[id].w);
                        sigvi.at(id)=0.;
                        delxm=dx*(nxdx-1);
                        dutotdxi.at(id)=0.;
                        dutotdyi.at(id)=0.;
                        elz.at(id)=kkar*eleff.at(id);
                        dutotdni.at(id)=dutotdzi.at(id);

                        detang(0,id,dutotds,dutotdn,i,j,k);
                        if(rcl>0){
                            phim=1.+4.7*rcl*eleff.at(id);
                            psim=-4.7*rcl*eleff.at(id);
                        }
                        else{
                            phim=pow( (1.-15.*rcl*eleff.at(id)),(-.25) );
                            psim=2.*log((1.+1./phim)/2.)+log((1.+1./pow(phim,2.f))/2.)-2.*atan(1./phim)+pi/2.;
                        }
                        ustarz.at(id)=elz.at(id)*dutotdni.at(id)/phim; // calculate local ustar
                        ustarz.at(id)=std::max(ustarz.at(id),3.e-02f);
                        u3psq=cusq*zbrac*ustarz.at(id)*ustarz.at(id);   //// (u''')^2
                        v3psq=cvsq*zbrac*ustarz.at(id)*ustarz.at(id);   //...
                        w3psq=cwsq*zbrac*ustarz.at(id)*ustarz.at(id); //...
                        upwp=-ctau13*zbrac*ustarz.at(id)*ustarz.at(id); // -tau13
                        upvp=0.;
                        vpwp=0.;
                        if(rcl<0 && zi.at(k)<.99*h){
                            u3psq=u3psq+.6*(ustarz.at(id)*ustarz.at(id))*pow( (-h*rcl),(2.f/3.f));
                            v3psq=v3psq+.6*(ustarz.at(id)*ustarz.at(id))*pow( (-h*rcl),(2.f/3.f));
                            w3psq=w3psq+3.3*(ustarz.at(id)*ustarz.at(id))*pow( (-zi.at(k)*rcl),(2.3f))*pow( (1.f-.8f*zi.at(k)/h),2.f);
                            upwp=upwp*pow( (1.f-zi.at(k)/h),(.5f*rcl*h/(1.0f-rcl*h)) );
                        }
                        else{
                            if(rcl<0){
                                u3psq=u3psq+.6*(ustarz.at(id)*ustarz.at(id))*pow( (-h*rcl),(2.f/3.f));
                                v3psq=v3psq+.6*(ustarz.at(id)*ustarz.at(id))*pow( (-h*rcl),(2.f/3.f));
                                w3psq=w3psq+3.3*(ustarz.at(id)*ustarz.at(id))*pow( (-.99f*h*rcl),(2.3f))*pow( (1.f-.8f*.99f),2.f);
                                upwp=upwp*pow( (1.f-.99f),(.5f*rcl*h/(1.0f-rcl*h)) );  
                            }
                        }
                        rotu3psq(id,u3psq,utot,upvp,upwp,vpwp,v3psq,w3psq,ufsqb,wfsqb,vfsqb,ufvf,ufwf,vfwf); // rotate values back into the orig. grid
                        ufwfi.at(id)=ufwf;
                        ufvfi.at(id)=ufvf;
                        vfwfi.at(id)=vfwf;
                        ustarij.at(id)=ustarz.at(id);
                        sigui.at(id)=sqrt(u3psq);
                        sigvi.at(id)=sqrt(v3psq);
                        sigwi.at(id)=sqrt(w3psq);
                        upwpi.at(id)=upwp;
                        // along the boundaries we make y effects negligible
                    }
                    else{
                        utot=sqrt(wind_vel[id].u*wind_vel[id].u+wind_vel[id].v*wind_vel[id].v+wind_vel[id].w*wind_vel[id].w);
                        // away from boundaries u*y, ly, dely, and gradients
                        int idim1=k*nxdx*nydy +j*nxdx +(i-1);
                        int idip1=k*nxdx*nydy +j*nxdx +(i+1);
                        if(cellQuic[idim1].c != 0 && cellQuic[id].c != 0  && cellQuic[idip1].c != 0){
                            //mdw 3-08-2004 start changes for highest gradient
                            utot=sqrt(wind_vel[id].u*wind_vel[id].u+wind_vel[id].v*wind_vel[id].v+wind_vel[id].w*wind_vel[id].w);
                            utotm=sqrt(wind_vel[idim1].u*wind_vel[idim1].u+wind_vel[idim1].v*wind_vel[idim1].v+wind_vel[idim1].w*wind_vel[idim1].w);
                            utotp=sqrt(wind_vel[idip1].u*wind_vel[idip1].u+wind_vel[idip1].v*wind_vel[idip1].v+wind_vel[idip1].w*wind_vel[idip1].w);

                            dutotdxp=(utotp-utot)/dx;
                            dutotdxm=(utot-utotm)/dx;
                            dutotdxa=std::max(fabs(dutotdxp),fabs(dutotdxm));
                            if(dutotdxa==fabs(dutotdxm)){
                                dutotdxi.at(id)=dutotdxm;
                            }
                            else{
                                dutotdxi.at(id)=dutotdxp;
                            }
                            // mdw 3-08-2004end changes
                        }
                        else{
                            if(cellQuic[id].c == 0){ ////BALLI
                                dutotdxi.at(id)=0.;
                            }
                            else{
                                if(cellQuic[idim1].c == 0){ ////BALLI
                                    dutotdxi.at(id)=2.*sqrt(wind_vel[id].u*wind_vel[id].u+wind_vel[id].v*wind_vel[id].v+wind_vel[id].w*wind_vel[id].w)/dx;
                                    dutotdxi.at(id)=sqrt(wind_vel[id].u*wind_vel[id].u+wind_vel[id].v*wind_vel[id].v+wind_vel[id].w*wind_vel[id].w)
                                        /(log((.5*dx)/z0)*(.5*dx));
                                }
                                else{
                                    dutotdxi.at(id)=-sqrt(wind_vel[id].u*wind_vel[id].u+wind_vel[id].v*wind_vel[id].v+wind_vel[id].w*wind_vel[id].w)
                                        /(log((.5*dx)/z0)*(.5*dx));
                                }
                            }
                        }
                        
                        int idjm1=k*nxdx*nydy +(j-1)*nxdx +i;
                        int idjp1=k*nxdx*nydy +(j+1)*nxdx +i;
                        if(cellQuic[id].c != 0 && cellQuic[idjm1].c != 0 && cellQuic[idjp1].c != 0){
                            //mdw 3-08-2008 start gradient changes
                            utot=sqrt(wind_vel[id].u*wind_vel[id].u+wind_vel[id].v*wind_vel[id].v+wind_vel[id].w*wind_vel[id].w);
                            utotm=sqrt(wind_vel[idjm1].u*wind_vel[idjm1].u+wind_vel[idjm1].v*wind_vel[idjm1].v+wind_vel[idjm1].w*wind_vel[idjm1].w);
                            utotp=sqrt(wind_vel[idjp1].u*wind_vel[idjp1].u+wind_vel[idjp1].v*wind_vel[idjp1].v+wind_vel[idjp1].w*wind_vel[idjp1].w);
                            dutotdyc=0.5*(utotp-utotm)/dy;
                            dutotdyp=(utotp-utot)/dy;
                            dutotdym=(utot-utotm)/dy;
                            dutotdya=std::max(fabs(dutotdyp),fabs(dutotdym));
                            if(dutotdya==fabs(dutotdym)){
                                dutotdyi.at(id)=dutotdym;
                            }
                            else{
                                dutotdyi.at(id)=dutotdyp;
                            }
                            // mdw 3-08-2004end changes
                        }
                        else{
                            if(cellQuic[id].c == 0){
                                dutotdyi.at(id)=0.;
                            }
                            else{
                                if(cellQuic[idjm1].c == 0){
                                    dutotdyi.at(id)=sqrt(wind_vel[id].u*wind_vel[id].u+wind_vel[id].v*wind_vel[id].v+wind_vel[id].w*wind_vel[id].w)
                                        /(log((.5*dy)/z0)*(.5*dy));
                                }
                                else{
                                    dutotdyi.at(id)=-sqrt(wind_vel[id].u*wind_vel[id].u+wind_vel[id].v*wind_vel[id].v+wind_vel[id].w*wind_vel[id].w)
                                        /(log((.5*dy)/z0)*(.5*dy));
                                }
                            }
                        }
                    }
                    detang(0,id,dutotds,dutotdn,i,j,k); // Calculates the parameters fot the triple rotatation of coord sys.
                    dwall=std::min(std::min(eleff.at(id),std::min(dxm.at(id),dxp.at(id))),std::min(dym.at(id),dyp.at(id)));
                    dwall=std::min(dwall,dzm.at(id));
                    elz.at(id)=kkar*dwall; // length scale based on distance to wall

                    if(cellQuic[id].c !=8)elz.at(id)=kkar*dwall; // length scale based on distance to wall
                    if(fabs(dutotdni.at(id))>1.e-6){
                        x_b=std::min(dxm.at(id),dxp.at(id));
                        if(x_b>std::max(del_b,dx)) x_b=0;
                        y_b=std::min(dym.at(id),dyp.at(id));
                        if(y_b>std::max(del_b,dy)) y_b=0;
                        dwallg=fabs(dutotdyi.at(id))*y_b+fabs(dutotdxi.at(id))*x_b;
                        dwallg=dwallg+fabs(dutotdzi.at(id))*eleff.at(id);
                        dwallg=dwallg/dutotdni.at(id);
                        elzv=kkar*utot/dutotdni.at(id); // length scale based on distance to null wind
                        if(dwallg*kkar<elzv && (x_b+y_b)>0.) {
                            // mdw 6-29-2006 changed test so that must be near vertical wall
                            elz.at(id)=kkar*dwallg; // pick the smallest length scale
                        }
                        else{
                            // mdw 6-30-2006 changed test so that vortex test does not override normal stuff
                            if(elzv<=elz.at(id)){
                                elz.at(id)=elzv;
                            }
                        }
                    }
                    if(cellQuic[id].c!=8){
                        if(rcl>0){
                            phim=1.+4.7*rcl*eleff.at(id);
                            psim=-4.7*rcl*eleff.at(id);
                        }
                        else{
                            phim=pow( (1.f-15.f*rcl*eleff.at(id)),(-.25f) );
                            psim=2.*log((1.+1./phim)/2.)+log((1.+1./pow(phim,2.f))/2.f)-2.*atan(1./phim)+pi/2.;
                        }
                    }
                    else{
                        phim=1.;
                        psim=0.;
                        elz.at(id)=elcanopy*std::min(1.f,(dzm.at(id)-z0)*kkar/(.3f*hgtveg.at(ij)));
                    }
                    ustarz.at(id)=elz.at(id)*dutotdni.at(id)/phim; // calculate local ustar
                    ustarz.at(id)=std::max(ustarz.at(id),3.e-02f);
                    //mdw 6-23-2004 adjust for vertical structure
                    u3psq=cusq*zbrac*ustarz.at(id)*ustarz.at(id);   // (u''')^2
                    v3psq=cvsq*zbrac*ustarz.at(id)*ustarz.at(id);   //...
                    w3psq=cwsq*zbrac*ustarz.at(id)*ustarz.at(id); //...
                    upwp=-ctau13*zbrac*ustarz.at(id)*ustarz.at(id); // -tau13
                    upvp=0.;
                    vpwp=0.;
                    xloc=1.;
                    if(ustarz.at(id)>xloc*ustarg.at(id)){
                        if(rcl<0. && zi.at(k)<.99*h){
                            u3psq=u3psq+.6*(ustarz.at(id)*ustarz.at(id))*pow( (-h*rcl),(2.f/3.f));
                            v3psq=v3psq+.6*(ustarz.at(id)*ustarz.at(id))*pow( (-h*rcl),(2.f/3.f));
                            w3psq=w3psq+3.3*(ustarz.at(id)*ustarz.at(id))*pow( (-zi.at(k)*rcl),(2.3f))*pow( (1.f-.8f*zi.at(k)/h),2.f);
                            upwp=upwp*pow( (1.f-zi.at(k)/h),(.5f*rcl*h/(1.f-rcl*h)) );
                        }
                        else{
                            if(rcl<0.){
                                u3psq=u3psq+.6*(ustarz.at(id)*ustarz.at(id))*pow( (-h*rcl),(2.f/3.f));
                                v3psq=v3psq+.6*(ustarz.at(id)*ustarz.at(id))*pow( (-h*rcl),(2.f/3.f));
                                w3psq=w3psq+3.3*(ustarz.at(id)*ustarz.at(id))*pow( (-.99f*h*rcl),(2.3f))*pow( (1.f-.8f*.99f),2.f);
                                upwp=upwp*pow( (1.f-.99f),(.5f*rcl*h/(1.f-rcl*h)) );
                            }
                        }
                        rotu3psq(id,u3psq,utot,upvp,upwp,vpwp,v3psq,w3psq,ufsqb,wfsqb,vfsqb,ufvf,ufwf,vfwf); // rotate values back into the orig. grid
                        ufwfi.at(id)=ufwf;
                        ufvfi.at(id)=ufvf;
                        vfwfi.at(id)=vfwf;
                        ustarij.at(id)=ustarz.at(id);
                        sigui.at(id)=sqrt(u3psq);
                        sigvi.at(id)=sqrt(v3psq);
                        sigwi.at(id)=sqrt(w3psq);
                        upwpi.at(id)=upwp;
                    }
                    else{ // non-local dominates (possible place for into of TPT)
                        ufsq=ufsqgi.at(id);
                        vfsq=vfsqgi.at(id);
                        wfsq=wfsqgi.at(id);
                        ufvf=ufvfgi.at(id);
                        ufwf=ufwfgi.at(id);
                        vfwf=vfwfgi.at(id);
                        sigui.at(id)=sqrt(ufsq);
                        sigvi.at(id)=sqrt(vfsq);
                        sigwi.at(id)=sqrt(wfsq);
                        ufwfi.at(id)=ufwf;
                        ufvf=0.;
                        ufvfi.at(id)=ufvf;
                        vfwfi.at(id)=vfwf;
                        //mdw 7-25-2005 corrections for axis rotation with non-local mixing
                        ustarij.at(id)=ustarg.at(id);
                        rotufsq(id,u3psq,upwp,v3psq,w3psq,ufsq,ufvf,ufwf,vfsq,vfwf,wfsq);
                        sigui.at(id)=sqrt(u3psq);
                        sigvi.at(id)=sqrt(v3psq);
                        sigwi.at(id)=sqrt(w3psq);
                        upwpi.at(id)=upwp;
                    }
                    sigu=sigui.at(id);
                    sigv=sigvi.at(id);
                    sigw=sigwi.at(id);
                    if(zi.at(k)<=.99*h){
                        eps=pow(ustarij.at(id),3.f)*(1.f-.75f*zi.at(k)*rcl)*pow((1.f-.85f*zi.at(k)/h),(1.5f))/eleff.at(id); // calculate epsilon for grid cell centers
                    }
                    else{
                        eps=pow(ustarij.at(id),3.f)*(1.f-.75f*.99f*h*rcl)*pow((1.f-.85f*.99f),(1.5f))/eleff.at(id); 
                    }
                    
                }
                epsi.at(id)=eps;
                
                // erp June 18, 2009  shader needs sigu sigv
                // convert them to taus and load them onto the shaders
                //tau13= ufwf; tau12 = ufvf; tau23 = vfwf need to be written to shaders
                turbfield<<std::setw(7) << std::setprecision(3)<< std::scientific<<xi.at(i)<<"   "
                         <<yi.at(j)<<"   "<<zi.at(k)<<"   "<< sigu<<"  "<<sigv<<"  "<<sigw<<"  "<<elz.at(id)<<"  "<<
                    eleff.at(id)<<"  "<<eps<<std::endl;
                
                sig[id].u = sigu;   //sigU
                sig[id].v = sigv;   //sigV
                sig[id].w = sigw;   //sigW
                
                float tau11=sigu*sigu;
                float tau22=sigv*sigv;
                float tau33=sigw*sigw;
                float tau13=ustarij.at(id)*ustarij.at(id);
                
                tau[id].t11   = tau11;             //Tau11
                tau[id+1].t22 = tau22;             //Tau22
                tau[id+2].t33 = tau33;             //Tau33
                tau[id+3].t13 = tau13;             //Tau13
                //Make tau's a texture so that they can be visualized as horizontal layers in the domain

                
                dataTwo[texidx]   =  sigu;
                dataTwo[texidx+1] =  sigv;
                dataTwo[texidx+2] =  sigw;
                dataTwo[texidx+3] =  ustarij.at(id);
                
                data3[texidx]   = dutotdxi.at(id);
                data3[texidx+1] = dutotdyi.at(id);
                data3[texidx+2] = dutotdzi.at(id);
                data3[texidx+3] = 0.;

                data4[texidx]   = dzm.at(id);
                data4[texidx+1] = dzp.at(id);
                data4[texidx+2] = dym.at(id);
                data4[texidx+3] = dyp.at(id);
                
                data[texidx]   = dxm.at(id);
                data[texidx+1] = dxp.at(id);
                data[texidx+2] = .00987;
                data[texidx+3] = .00986;
                
                ustar=std::max(ustarij.at(id),0.003f);
                float ustar3=ustar*ustar*ustar;
                dataWind[texidx+3] = 0.;
            }//   lp027
        }//      lp028
    }//         lp029
    
    createTexture(lambda, GL_RGBA32F_ARB, width,height, dataTwo);
    createTexture(windField, GL_RGBA32F_ARB, width, height, dataWind);
    createTexture(tauTex, GL_RGBA32F_ARB, width, height, dataTau);
    createTexture(tau_dz, GL_RGBA32F_ARB, width,height, data);
    
    delete[] dataTwo;
    delete[] dataWind;
    delete[] dataTau;
    delete[] data;
    std::vector<float>dsigwdni,dsigvdni,dsigudni,dupwpdni;
    dsigwdni.resize(nzdz*nydy*nxdx);
    dsigvdni.resize(nzdz*nydy*nxdx);
    dsigudni.resize(nzdz*nydy*nxdx);
    dupwpdni.resize(nzdz*nydy*nxdx);
    ani.resize(nzdz*nydy*nxdx);
    bni.resize(nzdz*nydy*nxdx);
    cni.resize(nzdz*nydy*nxdx);

    std::ofstream dsigwdnOut;
    dsigwdnOut.open("dsigwdn.dat");
    float dsigwdx=0.;
    float dsigwdy=0.;
    float dsigudx=0.;
    float dsigudy=0.;
    float dsigvdx=0.;
    float dsigvdy=0.;
    float dupwpdx=0.;
    float dupwpdy=0.;
    float dsigwdz=0.;
    float dsigvdz=0.;
    float dsigudz=0.;
    float dupwpdz=0.;
    float dsigwdn=0.;
    float dsigvdn=0.;
    float dsigudn=0.;
    float dupwpdn=0.;

    for(int k=0; k<nzdz;k++){//do k=2,nz-1
        for(int j=0;j<nydy;j++){//do j=1,ny-1
            for(int i=0;i<nxdx;i++){//do i=1,nx-1
                int id=k*nxdx*nydy + j*nxdx + i;
                
                int row = k / (numInRow);
                int texidx = row * width * nydy * 4 +
                    j * width * 4 +
                    k % (numInRow) * nxdx * 4 +
                    i * 4;
                
                
                if(cellQuic[id].c != 0){
                    int idim1=k*nxdx*nydy +j*nxdx +(i-1);
                    int idip1=k*nxdx*nydy +j*nxdx +(i+1);
                    
                    int idjm1=k*nxdx*nydy +(j-1)*nxdx +i;
                    int idjp1=k*nxdx*nydy +(j+1)*nxdx +i;
                    
                    int idkm1=(k-1)*nxdx*nydy +j*nxdx +i;
                    int idkp1=(k+1)*nxdx*nydy +j*nxdx +i;
                    
                    if(j<1||j>=nydy-1||i<1||i>=nxdx-1){
                        dsigwdx=0.;
                        dsigwdy=0.;
                        dsigudx=0.;
                        dsigudy=0.;
                        dsigvdx=0.;
                        dsigvdy=0.;
                        dupwpdx=0.;
                        dupwpdy=0.;
                    }
                    //
                    // calculate the gradient of sigma w normal to the flow using a CDD
                    //
                    utot=sqrt(wind_vel[id].u*wind_vel[id].u+wind_vel[id].v*wind_vel[id].v+wind_vel[id].w*wind_vel[id].w);
                    if(dxm.at(id)>=dx && dxp.at(id)>=dx && i!=0 && i!=nxdx-1){
                        if(idim1<0){
                            dsigwdx=.5*(sigwi.at(idip1)-sigwi.at(id))/dx;
                            dsigvdx=.5*(sigvi.at(idip1)-sigvi.at(id))/dx;
                            dsigudx=.5*(sigui.at(idip1)-sigui.at(id))/dx;
                            dupwpdx=.5*(upwpi.at(idip1)-upwpi.at(id))/dx;
                        }
                        else{
                            dsigwdx=.5*(sigwi.at(idip1)-sigwi.at(idim1))/dx;
                            dsigvdx=.5*(sigvi.at(idip1)-sigvi.at(idim1))/dx;
                            dsigudx=.5*(sigui.at(idip1)-sigui.at(idim1))/dx;
                            dupwpdx=.5*(upwpi.at(idip1)-upwpi.at(idim1))/dx;							  
                        }
                    }
                    else{
                        if(i==1||i==nxdx-1){
                            dsigwdx=0.;
                            dsigvdx=0.;
                            dsigudx=0.;
                            dupwpdx=0.;
                        }
                        else{
                            if(dxm.at(id)<dx && dxp.at(id)>dx){
                                //mdw 11-21-2005 modified if statements to address particle in 3-walled cells
                                dsigwdni.at(id)=(sigwi.at(idip1)-sigwi.at(id))/dx;
                                dsigvdni.at(id)=(sigvi.at(idip1)-sigvi.at(id))/dx;
                                dsigudni.at(id)=(sigui.at(idip1)-sigui.at(id))/dx;
                                dupwpdni.at(id)=(upwpi.at(idip1)-upwpi.at(id))/dx;
                                dsigwdni.at(id)=0.;
                                dsigvdni.at(id)=0.;
                                dsigudni.at(id)=0.;
                                dupwpdni.at(id)=0.;
                                sigwi.at(id)=std::max(sigwi.at(idip1),sigwi.at(id));
                                sigvi.at(id)=std::max(sigvi.at(idip1),sigvi.at(id));
                                sigui.at(id)=std::max(sigui.at(idip1),sigui.at(id));
                                ustarij.at(id)=std::max(ustarij.at(idip1),ustarij.at(id));
                                if(fabs(upwpi.at(id))<fabs(upwpi.at(idip1))){
                                    upwpi.at(id)=upwpi.at(idip1);
                                }
                                else{
                                    upwpi.at(idip1)=upwpi.at(id);
                                }
                            }
                            if(dxp.at(id)<dx && dxm.at(id)>dx){
                                //mdw 11-21-2005 modified if statements to address particle in 3-walled cells
                                if(idim1<0){
                                    dsigwdni.at(id)=0.;//(sigwi.at(id)-sigwi.at(id))/dx;
                                    dsigvdni.at(id)=0.;//(sigvi.at(id)-sigvi.at(id))/dx;
                                    dsigudni.at(id)=0.;//(sigui.at(id)-sigui.at(id))/dx;
                                    dupwpdni.at(id)=0.;//(upwpi.at(id)-upwpi.at(id))/dx;
                                    sigwi.at(id)=std::max(0.f,sigwi.at(id));
                                    sigvi.at(id)=std::max(0.f,sigvi.at(idim1));
                                    sigui.at(id)=std::max(0.f,sigui.at(id));
                                    ustarij.at(id)=std::max(0.f,ustarij.at(id));
                                }
                                else{
                                    dsigwdni.at(id)=(sigwi.at(idim1)-sigwi.at(id))/dx;
                                    dsigvdni.at(id)=(sigvi.at(idim1)-sigvi.at(id))/dx;
                                    dsigudni.at(id)=(sigui.at(idim1)-sigui.at(id))/dx;
                                    dupwpdni.at(id)=(upwpi.at(idim1)-upwpi.at(id))/dx;
                                    sigwi.at(id)=std::max(sigwi.at(idim1),sigwi.at(id));
                                    sigvi.at(id)=std::max(sigvi.at(idim1),sigvi.at(idim1));
                                    sigui.at(id)=std::max(sigui.at(idim1),sigui.at(id));
                                    ustarij.at(id)=std::max(ustarij.at(idim1),ustarij.at(id));
                                }
                                dsigwdni.at(id)=0.f;
                                dsigvdni.at(id)=0.f;
                                dsigudni.at(id)=0.f;
                                dupwpdni.at(id)=0.f;
                                
                                if(idim1<0){
                                    if(fabs(upwpi.at(id))<fabs(0.f)){
                                        upwpi.at(id)=upwpi.at(id);
                                    }
                                }
                                else{
                                    if(fabs(upwpi.at(id))<fabs(upwpi.at(idim1))){
                                        upwpi.at(id)=upwpi.at(idim1);
                                    }
                                    else{
                                        upwpi.at(idim1)=upwpi.at(id);
                                    }
                                }
                            }  
                        }
                    }
                    
                    if(dym.at(id)>=dy && dyp.at(id)>=dy && j!=0 && j!=ny-1){
                        if(idjm1<0){
                            dsigwdy=.5*(sigwi.at(idjp1)-sigwi.at(id))/dy;
                            dsigvdy=.5*(sigvi.at(idjp1)-sigvi.at(id))/dy;
                            dsigudy=.5*(sigui.at(idjp1)-sigui.at(id))/dy;
                            dupwpdy=.5*(upwpi.at(idjp1)-upwpi.at(id))/dy;
                        }
                        else{
                            dsigwdy=.5*(sigwi.at(idjp1)-sigwi.at(idjm1))/dy;
                            dsigvdy=.5*(sigvi.at(idjp1)-sigvi.at(idjm1))/dy;
                            dsigudy=.5*(sigui.at(idjp1)-sigui.at(idjm1))/dy;
                            dupwpdy=.5*(upwpi.at(idjp1)-upwpi.at(idjm1))/dy;
                        }
                    }
                    else{
                        if(j==1||j==ny-1){
                            dsigwdy=0.;
                        }
                        else{
                            if(dym.at(id)<dy && dyp.at(id)>dy){
                                //mdw 11-21-2006 modified if statements to address particle in 3-walled cells
                                dsigwdni.at(id)=(sigwi.at(idjp1)-sigwi.at(id))/dy;
                                dsigvdni.at(id)=(sigvi.at(idjp1)-sigvi.at(id))/dy;
                                dsigudni.at(id)=(sigui.at(idjp1)-sigui.at(id))/dy;
                                dupwpdni.at(id)=(upwpi.at(idjp1)-upwpi.at(id))/dy;
                                dsigwdni.at(id)=0.;
                                dsigvdni.at(id)=0.;
                                dsigudni.at(id)=0.;
                                dupwpdni.at(id)=0.;
                                sigwi.at(id)=std::max(sigwi.at(idjp1),sigwi.at(id));
                                sigvi.at(id)=std::max(sigvi.at(idjp1),sigvi.at(id));
                                sigui.at(id)=std::max(sigui.at(idjp1),sigui.at(id));
                                ustarij.at(id)=std::max(ustarij.at(idjp1),ustarij.at(id));
                                if(fabs(upwpi.at(id))<fabs(upwpi.at(idjp1))){
                                    upwpi.at(id)=upwpi.at(idjp1);
                                }
                                else{
                                    upwpi.at(idjp1)=upwpi.at(id);
                                }
                            }
                            if(dyp.at(id)<dy && dym.at(id)>dy){
                                //mdw 11-21-2005 modified if statements to address particle in 3-walled cells
                                if(idjm1<0){
                                    dsigwdni.at(id)=0.;//(sigwi.at(id)-sigwi.at(id))/dy;
                                    dsigvdni.at(id)=0.;//(sigvi.at(id)-sigvi.at(id))/dy;
                                    dsigudni.at(id)=0.;//(sigui.at(id)-sigui.at(id))/dy;
                                    dupwpdni.at(id)=0.;//(upwpi.at(id)-upwpi.at(id))/dy;
                                    
                                    sigwi.at(id)=std::max(0.f,sigwi.at(id));
                                    sigvi.at(id)=std::max(0.f,sigvi.at(id));
                                    sigui.at(id)=std::max(0.f,sigui.at(id));
                                    ustarij.at(id)=std::max(0.f,ustarij.at(id));
                                }
                                else{
                                    dsigwdni.at(id)=(sigwi.at(idjm1)-sigwi.at(id))/dy;
                                    dsigvdni.at(id)=(sigvi.at(idjm1)-sigvi.at(id))/dy;
                                    dsigudni.at(id)=(sigui.at(idjm1)-sigui.at(id))/dy;
                                    dupwpdni.at(id)=(upwpi.at(idjm1)-upwpi.at(id))/dy;
                                    
                                    sigwi.at(id)=std::max(sigwi.at(idjm1),sigwi.at(id));
                                    sigvi.at(id)=std::max(sigvi.at(idjm1),sigvi.at(id));
                                    sigui.at(id)=std::max(sigui.at(idjm1),sigui.at(id));
                                    ustarij.at(id)=std::max(ustarij.at(idjm1),ustarij.at(id));
                                }
                                
                                dsigwdni.at(id)=0.f;
                                dsigvdni.at(id)=0.f;
                                dsigudni.at(id)=0.f;
                                dupwpdni.at(id)=0.f;
                                if(idjm1<0){
                                    if(fabs(upwpi.at(id))<fabs(0.f)){
                                        upwpi.at(id)=upwpi.at(idjm1);
                                    }
                                }
                                else{
                                    if(fabs(upwpi.at(id))<fabs(upwpi.at(idjm1))){
                                        upwpi.at(id)=upwpi.at(idjm1);
                                    }
                                    else{
                                        upwpi.at(idjm1)=upwpi.at(id);
                                    }
                                }
				
                            }
                        }
                    }
                    if(dzm.at(id)>dz && k!=nzdz-1 && dzp.at(id)>dz){
                        if(idkm1<0){
                            dsigwdz=(sigwi.at(idkp1))/(dz_array.at(k)+0.5+dz_array.at(k+1));
                            dsigvdz=(sigvi.at(idkp1))/(dz_array.at(k)+0.5+dz_array.at(k+1));
                            dsigudz=(sigui.at(idkp1))/(dz_array.at(k)+0.5+dz_array.at(k+1));
                            dupwpdz=(upwpi.at(idkp1))/(dz_array.at(k)+0.5+dz_array.at(k+1));
                        }
                        else{
                            dsigwdz=(sigwi.at(idkp1)-sigwi.at(idkm1))/(dz_array.at(k)+.5*dz_array.at(k-1)+dz_array.at(k+1));
                            dsigvdz=(sigvi.at(idkp1)-sigvi.at(idkm1))/(dz_array.at(k)+.5*dz_array.at(k-1)+dz_array.at(k+1));
                            dsigudz=(sigui.at(idkp1)-sigui.at(idkm1))/(dz_array.at(k)+.5*dz_array.at(k-1)+dz_array.at(k+1));
                            dupwpdz=(upwpi.at(idkp1)-upwpi.at(idkm1))/(dz_array.at(k)+.5*dz_array.at(k-1)+dz_array.at(k+1));
                        }
                    }
                    if(dzm.at(id)<=dz  && k!=nzdz-1 && dzp.at(id)>dz){
                        dsigwdn=(sigwi.at(idkp1)-sigwi.at(id))/(dz_array.at(k)+.5+dz_array.at(k+1));
                        dsigvdn=(sigvi.at(idkp1)-sigvi.at(id))/(dz_array.at(k)+.5+dz_array.at(k+1));
                        dsigudn=(sigui.at(idkp1)-sigui.at(id))/(dz_array.at(k)+.5+dz_array.at(k+1));
                        dupwpdn=(upwpi.at(idkp1)-upwpi.at(id))/(dz_array.at(k)+.5+dz_array.at(k+1));
                        //mdw 9-26-2005 force dsigwdn to be zero near surface
                        dsigwdn=0.;
                        dsigudn=0.;
                        dsigvdn=0.;
                        dupwpdn=0.;
                        sigui.at(id)=std::max(sigui.at(idkp1),sigui.at(id));
                        sigvi.at(id)=std::max(sigvi.at(idkp1),sigvi.at(id));
                        sigwi.at(id)=std::max(sigwi.at(idkp1),sigwi.at(id));
                        ustarij.at(id)=std::max(ustarij.at(idkp1),ustarij.at(id));
                        if(fabs(upwpi.at(id))<fabs(upwpi.at(idkp1))){
                            upwpi.at(id)=upwpi.at(idkp1);
                        }
                        else{
                            upwpi.at(idkp1)=upwpi.at(id);
                        }
                    }
                    if(dzp.at(id)<=dz && k!=nzdz-1 && dzm.at(id)>dz){
                        //mdw 9-26-2005 force dsigwdn to be zero near surface
                        dsigwdn=0.;
                        dsigudn=0.;
                        dsigvdn=0.;
                        dupwpdn=0.;
                        if(idkm1<0){
                            sigui.at(id)=std::max(0.f,sigui.at(id));
                            sigvi.at(id)=std::max(0.f,sigvi.at(id));
                            sigwi.at(id)=std::max(0.f,sigwi.at(id));
                            ustarij.at(id)=std::max(0.f,ustarij.at(id));
                        }
                        else{
                            sigui.at(id)=std::max(sigui.at(idkm1),sigui.at(id));
                            sigvi.at(id)=std::max(sigvi.at(idkm1),sigvi.at(id));
                            sigwi.at(id)=std::max(sigwi.at(idkm1),sigwi.at(id));
                            ustarij.at(id)=std::max(ustarij.at(idkm1),ustarij.at(id));
                        }
                        if(idkm1<0){
                            if(fabs(upwpi.at(id))<0.f){
                                upwpi.at(id)=0.;
                            }
                        }
                        else{
                            if(fabs(upwpi.at(id))<fabs(upwpi.at(idkm1))){
                                upwpi.at(id)=upwpi.at(idkm1);
                            }
                            else{
                                upwpi.at(idkm1)=upwpi.at(id);
                            }
                        }
                        
                    }
                    // erp June 18, 2009: only need diswdni used to calculate tau gradients
                    // export this value to be used, send to shaders
                    if((dxm.at(id)>=dx)&&(dxp.at(id)>=dx)&&(dym.at(id)>=dy)&& 
                       (dyp.at(id)>=dy)&&(dzm.at(id)>=dz)&&(dzp.at(id)>=dz)){
                        dsigwdn=ani.at(id)*dsigwdx+bni.at(id)*dsigwdy+cni.at(id)*dsigwdz;
                        dsigvdn=ani.at(id)*dsigvdx+bni.at(id)*dsigvdy+cni.at(id)*dsigvdz;
                        dsigudn=ani.at(id)*dsigudx+bni.at(id)*dsigudy+cni.at(id)*dsigudz;
                        dupwpdn=ani.at(id)*dupwpdx+bni.at(id)*dupwpdy+cni.at(id)*dupwpdz;
                    }
                    dsigwdni.at(id)=dsigwdn;
                    dsigvdni.at(id)=dsigvdn;
                    dsigudni.at(id)=dsigudn;
                    dupwpdni.at(id)=dupwpdn;
                    
                    // limiting form for near wall circumstances
                    data3[texidx+3] = dsigwdni.at(id);
                }
                dsigwdnOut<<xi.at(i)<<" "<<yi.at(j)<<" "<<zi.at(k)<<" "<<dsigwdni.at(id)<<std::endl;
            }//   lp030
        }//   lp031
    } //  lp032
    //500 line code ends
    createTexture(duvw_dz, GL_RGBA32F_ARB, width,height, data3);
    createTexture(dxyz_wall, GL_RGBA32F_ARB, width,height, data4);
    delete[] data3;
    delete[] data4;
        
}//end turbinit

//rotation subroutines

void ParticleControl:: detang(const int iomega,const int &id, float &dutotds,float &dutotdn,int &i,int &j,int &k){
    float  e11,e12,e13,e21,e22,e23,e31,e32,e33;
    float cospsi,sinpsi,sinphiw,cosphiw,omenum,omeden,cosomeg,sinomeg,tanomeg,omeg1,cosomeg1,sinomeg1;
    float omeg2,cosomeg2,sinomeg2,dutotdn2,dutotdn1,pi,omeg;

	pi=4.*atan(1.0);
    
    if(sqrt(wind_vel[id].u*wind_vel[id].u+wind_vel[id].v*wind_vel[id].v)>1.e-05){
        cospsi=wind_vel[id].u/(sqrt(wind_vel[id].u*wind_vel[id].u+wind_vel[id].v*wind_vel[id].v));
        sinpsi=wind_vel[id].v/(sqrt(wind_vel[id].u*wind_vel[id].u+wind_vel[id].v*wind_vel[id].v));
        sinphiw=wind_vel[id].w/(sqrt(wind_vel[id].u*wind_vel[id].u+wind_vel[id].v*wind_vel[id].v+wind_vel[id].w*wind_vel[id].w)+1.e-10);
        cosphiw=sqrt(wind_vel[id].u*wind_vel[id].u+wind_vel[id].v*wind_vel[id].v)/(sqrt(wind_vel[id].u*wind_vel[id].u+wind_vel[id].v*wind_vel[id].v+wind_vel[id].w*wind_vel[id].w)+1.e-10);
        omenum=-dutotdxi.at(id)*sinpsi+dutotdyi.at(id)*cospsi;
        omeden=dutotdxi.at(id)*cospsi*sinphiw+dutotdyi.at(id)*sinpsi*sinphiw-dutotdzi.at(id)*cosphiw;
        if(iomega==0){
            if(fabs(omeden)<1.e-10){
                cosomeg=0.;
                sinomeg=1.;
                  dutotdn=dutotdxi.at(id)*(sinpsi*sinomeg-cospsi*sinphiw*cosomeg) 
                       -dutotdyi.at(id)*(cospsi*sinomeg+sinpsi*sinphiw*cosomeg) 
                      +dutotdzi.at(id)*cosomeg*cosphiw;
                  if(dutotdn<0.)sinomeg=-1.;
            }
            else{
                tanomeg=omenum/omeden;
				omeg1=atan(tanomeg);
                cosomeg1=cos(omeg1);
                sinomeg1=sin(omeg1);
                dutotdn1=dutotdxi.at(id)*(sinpsi*sinomeg1-cospsi*sinphiw*cosomeg1) 
                    -dutotdyi.at(id)*(cospsi*sinomeg1+sinpsi*sinphiw*cosomeg1) 
                    +dutotdzi.at(id)*cosomeg1*cosphiw;
                omeg2=omeg1+pi;
                cosomeg2=cos(omeg2);
                sinomeg2=sin(omeg2);
                dutotdn2=dutotdxi.at(id)*(sinpsi*sinomeg2-cospsi*sinphiw*cosomeg2) 
                    -dutotdyi.at(id)*(cospsi*sinomeg2+sinpsi*sinphiw*cosomeg2) 
                    +dutotdzi.at(id)*cosomeg2*cosphiw;
                if(dutotdn2>dutotdn1){
                    dutotdn=dutotdn2;
                    omeg=omeg2;
                    cosomeg=cosomeg2;
                    sinomeg=sinomeg2;
                }
                else{
                    dutotdn=dutotdn1;
                    omeg=omeg1;
                    cosomeg=cosomeg1;
                    sinomeg=sinomeg1;
                }
                

            }
        }
        else{
            if(iomega==1){
                omeg=0.;
                cosomeg=1.;
                sinomeg=0.;
            }
            else{
                if(fabs(cospsi)>0.5){
                    omeg=pi/2.;
                    cosomeg=0.;
					sinomeg=-sign(cospsi,1.0);
                    if(iomega==3){
                        sinomeg=sign(cospsi,1.0);
                        omeg=-omeg;
                    }
                }
                else{
                    return;
                }
            }
        }
        alph1ij.at(id)=cospsi*cosphiw;
        alph2ij.at(id)=-sinpsi*cosomeg-cospsi*sinphiw*sinomeg;
        alph3ij.at(id)=sinpsi*sinomeg-cospsi*sinphiw*cosomeg;
        bet1ij.at(id)=sinpsi*cosphiw;
        bet2ij.at(id)=cospsi*cosomeg-sinpsi*sinphiw*sinomeg;
        bet3ij.at(id)=-cospsi*sinomeg-sinpsi*sinphiw*cosomeg;
        gam1ij.at(id)=sinphiw;
        gam2ij.at(id)=cosphiw*sinomeg;
        gam3ij.at(id)=cosphiw*cosomeg;
        alphn1ij.at(id)=cospsi*cosphiw;
        alphn2ij.at(id)=sinpsi*cosphiw;
        alphn3ij.at(id)=sinphiw;
        betn1ij.at(id)=-sinpsi*cosomeg-cospsi*sinphiw*sinomeg;
        betn2ij.at(id)=cospsi*cosomeg-sinpsi*sinphiw*sinomeg;
        betn3ij.at(id)=cosphiw*sinomeg;
        gamn1ij.at(id)=sinpsi*sinomeg-cospsi*sinphiw*cosomeg;
        gamn2ij.at(id)=-cospsi*sinomeg-sinpsi*sinphiw*cosomeg;
        gamn3ij.at(id)=cosphiw*cosomeg;
        dutotdn=dutotdxi.at(id)*(sinpsi*sinomeg-cospsi*sinphiw*cosomeg) 
            -dutotdyi.at(id)*(cospsi*sinomeg+sinpsi*sinphiw*cosomeg) 
            +dutotdzi.at(id)*cosomeg*cosphiw;
        
        dutotds=dutotdxi.at(id)*cospsi*cosphiw+dutotdyi.at(id)*
            sinpsi*cosphiw+dutotdzi.at(id)*sinphiw;
        dutotdni.at(id)=dutotdn;
        dutotdsi.at(id)=dutotds;
        ani.at(id)=(sinpsi*sinomeg-cospsi*sinphiw*cosomeg);
        bni.at(id)=-(cospsi*sinomeg+sinpsi*sinphiw*cosomeg);
        cni.at(id)=cosomeg*cosphiw;
    }
    else{
        if(fabs(wind_vel[id].w)<1.e-05){
            if(sqrt(dutotdxi.at(id)*dutotdxi.at(id)+dutotdyi.at(id)*dutotdyi.at(id))>0.){
                cospsi=dutotdxi.at(id)/(sqrt(dutotdxi.at(id)*dutotdxi.at(id)+dutotdyi.at(id)*dutotdyi.at(id)));
                sinpsi=dutotdyi.at(id)/(sqrt(dutotdxi.at(id)*dutotdxi.at(id)+dutotdyi.at(id)*dutotdyi.at(id)));
            }
            else{
                cospsi=1.;
                sinpsi=0.;
            }
            if(sqrt(dutotdxi.at(id)*dutotdxi.at(id)+dutotdyi.at(id)*dutotdyi.at(id)+dutotdzi.at(id)*dutotdzi.at(id))>0){
                cosphiw=dutotdzi.at(id)/(sqrt(dutotdxi.at(id)*dutotdxi.at(id)
                                              +dutotdyi.at(id)*dutotdyi.at(id)+dutotdzi.at(id)*dutotdzi.at(id)));
                  sinphiw=sqrt(dutotdxi.at(id)*dutotdxi.at(id)+dutotdyi.at(id)*dutotdyi.at(id))
                      /(sqrt(dutotdxi.at(id)*dutotdxi.at(id)+dutotdyi.at(id)*dutotdyi.at(id)+dutotdzi.at(id)*dutotdzi.at(id)));
            }
            else{
                cosphiw=1.;
                sinphiw=0.;
            }
            alphn1ij.at(id)=cospsi*cosphiw;
            alphn2ij.at(id)=sinpsi*cosphiw;
            alphn3ij.at(id)=-sinphiw;
            betn1ij.at(id)=-sinpsi;
            betn2ij.at(id)=cospsi;
            betn3ij.at(id)=0.;
            gamn1ij.at(id)=sinphiw*cospsi;
            gamn2ij.at(id)=sinphiw*sinpsi;
            gamn3ij.at(id)=cosphiw;
            alph1ij.at(id)=cospsi*cosphiw;
            alph2ij.at(id)=-sinpsi;
            alph3ij.at(id)=sinphiw*cospsi;
            bet1ij.at(id)=sinpsi*cosphiw;
            bet2ij.at(id)=cospsi;
            bet3ij.at(id)=sinphiw*sinpsi;
            gam1ij.at(id)=-sinphiw;
            gam2ij.at(id)=0.;
            gam3ij.at(id)=cosphiw;
            dutotdni.at(id)=sqrt(dutotdxi.at(id)*dutotdxi.at(id)+dutotdyi.at(id)*dutotdyi.at(id)+dutotdzi.at(id)*dutotdzi.at(id));
            dutotdsi.at(id)=0.;
            ani.at(id)=sinphiw*cospsi;
            bni.at(id)=sinphiw*sinpsi;
            cni.at(id)=cosphiw;
        }
        else{
            if(wind_vel[id].w>0.){
                if(sqrt(dutotdxi.at(id)*dutotdxi.at(id)+dutotdyi.at(id)*dutotdyi.at(id))>0.){
                    cospsi=dutotdxi.at(id)/sqrt(dutotdxi.at(id)*dutotdxi.at(id)+dutotdyi.at(id)*dutotdyi.at(id));
                    sinpsi=dutotdyi.at(id)/sqrt(dutotdxi.at(id)*dutotdxi.at(id)+dutotdyi.at(id)*dutotdyi.at(id));
                }
                else{
                    cospsi=1.;
                    sinpsi=0.;
                }
                alph1ij.at(id)=0.;
                alph2ij.at(id)=sinpsi;
                alph3ij.at(id)=cospsi;
                bet1ij.at(id)=0.;
                bet2ij.at(id)=-cospsi;
                bet3ij.at(id)=sinpsi;
                gam1ij.at(id)=1.;
                gam2ij.at(id)=0.;
                gam3ij.at(id)=0.;
                alphn1ij.at(id)=0.;
                alphn2ij.at(id)=0.;
                alphn3ij.at(id)=1.;
                betn1ij.at(id)=sinpsi;
                betn2ij.at(id)=-cospsi;
                betn3ij.at(id)=0.;
                gamn1ij.at(id)=cospsi;
                gamn2ij.at(id)=sinpsi;
                gamn3ij.at(id)=0.;
                dutotdni.at(id)=sqrt(dutotdxi.at(id)*dutotdxi.at(id)+dutotdyi.at(id)*dutotdyi.at(id));
                dutotdni.at(id)=std::max(1.e-12f,dutotdni.at(id));
                
                dutotdsi.at(id)=dutotdzi.at(id);
                ani.at(id)=cospsi;
                bni.at(id)=sinpsi;
                cni.at(id)=0.;
            }
            else{
                if(sqrt(dutotdxi.at(id)*dutotdxi.at(id)+dutotdyi.at(id)*dutotdyi.at(id))>0.){
                    cospsi=dutotdxi.at(id)/sqrt(dutotdxi.at(id)*dutotdxi.at(id)+dutotdyi.at(id)*dutotdyi.at(id));
                    sinpsi=dutotdyi.at(id)/sqrt(dutotdxi.at(id)*dutotdxi.at(id)+dutotdyi.at(id)*dutotdyi.at(id));
                }
                else{
                    cospsi=1.;
                    sinpsi=0.;
                  }
                alphn1ij.at(id)=0.;
                alphn2ij.at(id)=0.;
                alphn3ij.at(id)=-1.;
                betn1ij.at(id)=-sinpsi;
                betn2ij.at(id)=cospsi;
                betn3ij.at(id)=0.;
                gamn1ij.at(id)=cospsi;
                gamn2ij.at(id)=sinpsi;
                gamn3ij.at(id)=0.;
                alph1ij.at(id)=0.;
                alph2ij.at(id)=-sinpsi;
                alph3ij.at(id)=cospsi;
                bet1ij.at(id)=0.;
                bet2ij.at(id)=cospsi;
                bet3ij.at(id)=sinpsi;
                gam1ij.at(id)=-1.;
                gam2ij.at(id)=0.;
                gam3ij.at(id)=0.;
                dutotdni.at(id)=sqrt(dutotdxi.at(id)*dutotdxi.at(id)+dutotdyi.at(id)*dutotdyi.at(id));
                
                dutotdsi.at(id)=-dutotdzi.at(id);
            }
        }
    }
    e11=alph1ij.at(id)*alphn1ij.at(id)+alph2ij.at(id)*betn1ij.at(id)+alph3ij.at(id)*gamn1ij.at(id);
    e12=alph1ij.at(id)*alphn2ij.at(id)+alph2ij.at(id)*betn2ij.at(id)+alph3ij.at(id)*gamn2ij.at(id);
    e13=alph1ij.at(id)*alphn3ij.at(id)+alph2ij.at(id)*betn3ij.at(id)+alph3ij.at(id)*gamn3ij.at(id);
    e21=bet1ij.at(id)*alphn1ij.at(id)+bet2ij.at(id)*betn1ij.at(id)+bet3ij.at(id)*gamn1ij.at(id);
    e22=bet1ij.at(id)*alphn2ij.at(id)+bet2ij.at(id)*betn2ij.at(id)+bet3ij.at(id)*gamn2ij.at(id);
    e23=bet1ij.at(id)*alphn3ij.at(id)+bet2ij.at(id)*betn3ij.at(id)+bet3ij.at(id)*gamn3ij.at(id);
    e31=gam1ij.at(id)*alphn1ij.at(id)+gam2ij.at(id)*betn1ij.at(id)+gam3ij.at(id)*gamn1ij.at(id);
    e32=gam1ij.at(id)*alphn2ij.at(id)+gam2ij.at(id)*betn2ij.at(id)+gam3ij.at(id)*gamn2ij.at(id);
    e33=gam1ij.at(id)*alphn3ij.at(id)+gam2ij.at(id)*betn3ij.at(id)+gam3ij.at(id)*gamn3ij.at(id);
    
}
float ParticleControl::sign(const float &A,const float &B){
    float R;
    if(B>=0.f)
        R=fabs(A);
    else
        R=-fabs(A);
    return R;	

}
void ParticleControl::rotu3psq(const int &id,const float &u3psq,const float &utot,const float &upvp,const float &upwp
							   ,const float &vpwp,const float &v3psq,const float &w3psq,float &ufsqb,float &wfsqb,float &vfsqb
							   ,float &ufvf,float &ufwf,float &vfwf){
    //this subroutine rotates the fluctuating quanitities back into the normal
    // coordinate sytem
								   
    ufsqb=u3psq*alph1ij.at(id)*alph1ij.at(id)+utot*utot*alph1ij.at(id)*alph1ij.at(id)+2.*upvp*alph1ij.at(id)*alph2ij.at(id) 
        +2.*upwp*alph1ij.at(id)*alph3ij.at(id)+v3psq*alph2ij.at(id)*alph2ij.at(id) 
        +2.*vpwp*alph2ij.at(id)*alph3ij.at(id)+w3psq*alph3ij.at(id)*alph3ij.at(id)-wind_vel[id].u*wind_vel[id].u;
    wfsqb=u3psq*gam1ij.at(id)*gam1ij.at(id)+utot*utot*gam1ij.at(id)*gam1ij.at(id)+2.*upvp*gam1ij.at(id)*gam2ij.at(id) 
        +2.*upwp*gam1ij.at(id)*gam3ij.at(id)+v3psq*gam2ij.at(id)*gam2ij.at(id)+ 
        2.*vpwp*gam2ij.at(id)*gam3ij.at(id)+w3psq*gam3ij.at(id)*gam3ij.at(id)-wind_vel[id].w*wind_vel[id].w;
    vfsqb=u3psq*bet1ij.at(id)*bet1ij.at(id)+utot*utot*bet1ij.at(id)*bet1ij.at(id)+2.*upvp*bet1ij.at(id)*bet2ij.at(id) 
        +2.*upwp*bet1ij.at(id)*bet3ij.at(id)+v3psq*bet2ij.at(id)*bet2ij.at(id) 
        +2.*vpwp*bet2ij.at(id)*bet3ij.at(id)+w3psq*bet3ij.at(id)*bet3ij.at(id)-wind_vel[id].v*wind_vel[id].v;
    ufvf=u3psq*alph1ij.at(id)*bet1ij.at(id)+utot*utot*alph1ij.at(id)*bet1ij.at(id) 
        +upvp*(alph1ij.at(id)*bet2ij.at(id)+alph2ij.at(id)*bet1ij.at(id)) 
        +upwp*(alph1ij.at(id)*bet3ij.at(id)+alph3ij.at(id)*bet1ij.at(id)) 
        +v3psq*alph2ij.at(id)*bet2ij.at(id)+vpwp*(alph2ij.at(id)*bet3ij.at(id) 
        +alph3ij.at(id)*bet2ij.at(id))+w3psq*alph3ij.at(id)*bet3ij.at(id) 
        -wind_vel[id].u*wind_vel[id].v;
    ufwf=u3psq*alph1ij.at(id)*gam1ij.at(id)+utot*utot*alph1ij.at(id)*gam1ij.at(id) 
        +upvp*(alph1ij.at(id)*gam2ij.at(id)+alph2ij.at(id)*gam1ij.at(id)) 
        +upwp*(alph1ij.at(id)*gam3ij.at(id)+alph3ij.at(id)*gam1ij.at(id)) 
        +v3psq*alph2ij.at(id)*gam2ij.at(id)+vpwp*(alph2ij.at(id)*gam3ij.at(id) 
        +alph3ij.at(id)*gam2ij.at(id))+w3psq*alph3ij.at(id)*gam3ij.at(id) 
        -wind_vel[id].u*wind_vel[id].w;
    vfwf=u3psq*bet1ij.at(id)*gam1ij.at(id)+utot*utot*bet1ij.at(id)*gam1ij.at(id)+upvp* 
        (bet1ij.at(id)*gam2ij.at(id)+bet2ij.at(id)*gam1ij.at(id))+upwp*(bet1ij.at(id)*gam3ij.at(id) 
        +bet3ij.at(id)*gam1ij.at(id))+v3psq*bet2ij.at(id)*gam2ij.at(id) 
        +vpwp*(bet2ij.at(id)*gam3ij.at(id)+bet3ij.at(id)*gam2ij.at(id)) 
        +w3psq*bet3ij.at(id)*gam3ij.at(id)-wind_vel[id].v*wind_vel[id].w;
}

void ParticleControl:: rotufsq(const int &id,float &u3psq,float &upwp
							   ,float &v3psq,float &w3psq,const float &ufsq,const float &ufvf
							   ,const float &ufwf,const float &vfsq,const float &vfwf,const float &wfsq){
    float  e11,e12,e13,e21,e22,e23,e31,e32,e33;    
    
    u3psq=ufsq*alphn1ij.at(id)*alphn1ij.at(id)+wind_vel[id].u*wind_vel[id].u*alphn1ij.at(id)*alphn1ij.at(id) 
        +2.*ufvf*alphn1ij.at(id)*alphn2ij.at(id)+2.*wind_vel[id].u*wind_vel[id].v*alphn1ij.at(id) 
        *alphn2ij.at(id)+2.*ufwf*alphn1ij.at(id)*alphn3ij.at(id)+2.*wind_vel[id].u*wind_vel[id].w 
        *alphn1ij.at(id)*alphn3ij.at(id)+vfsq*alphn2ij.at(id)*alphn2ij.at(id)+wind_vel[id].v*wind_vel[id].v 
        *alphn2ij.at(id)*alphn2ij.at(id)+2.*vfwf*alphn2ij.at(id)*alphn3ij.at(id)+ 
        2.*wind_vel[id].v*wind_vel[id].w*alphn2ij.at(id)*alphn3ij.at(id)+wfsq 
        *alphn3ij.at(id)*alphn3ij.at(id)+wind_vel[id].w*wind_vel[id].w*alphn3ij.at(id)*alphn3ij.at(id)
        -(wind_vel[id].u*wind_vel[id].u+wind_vel[id].v*wind_vel[id].v +wind_vel[id].w*wind_vel[id].w);
    
    v3psq=ufsq*betn1ij.at(id)*betn1ij.at(id)+wind_vel[id].u*wind_vel[id].u*betn1ij.at(id)*betn1ij.at(id)+2.*ufvf 
        *betn1ij.at(id)*betn2ij.at(id)+2.*wind_vel[id].u*wind_vel[id].v*betn1ij.at(id) 
        *betn2ij.at(id)+2.*ufwf*betn1ij.at(id)*betn3ij.at(id)+2.*wind_vel[id].u 
        *wind_vel[id].w*betn1ij.at(id)*betn3ij.at(id)+vfsq*betn2ij.at(id)*betn2ij.at(id) 
        +wind_vel[id].v*wind_vel[id].v*betn2ij.at(id)*betn2ij.at(id)+2.*vfwf*betn2ij.at(id)*betn3ij.at(id) 
        +2.*wind_vel[id].v*wind_vel[id].w*betn2ij.at(id)*betn3ij.at(id)+wfsq*betn3ij.at(id)*betn3ij.at(id) 
        +wind_vel[id].w*wind_vel[id].w*betn3ij.at(id)*betn3ij.at(id);
    w3psq=ufsq*gamn1ij.at(id)*gamn1ij.at(id)+wind_vel[id].u*wind_vel[id].u*gamn1ij.at(id)*gamn1ij.at(id)+2.*ufvf 
        *gamn1ij.at(id)*gamn2ij.at(id)+2.*wind_vel[id].u*wind_vel[id].v*gamn1ij.at(id) 
        *gamn2ij.at(id)+2.*ufwf*gamn1ij.at(id)*gamn3ij.at(id)+2.*wind_vel[id].u 
        *wind_vel[id].w*gamn1ij.at(id)*gamn3ij.at(id)+vfsq*gamn2ij.at(id)*gamn2ij.at(id) 
        +wind_vel[id].v*wind_vel[id].v*gamn2ij.at(id)*gamn2ij.at(id)+2.*vfwf*gamn2ij.at(id)*gamn3ij.at(id) 
        +2.*wind_vel[id].v*wind_vel[id].w*gamn2ij.at(id)*gamn3ij.at(id)+wfsq*gamn3ij.at(id)*gamn3ij.at(id) 
        +wind_vel[id].w*wind_vel[id].w*gamn3ij.at(id)*gamn3ij.at(id);
    upwp=ufsq*alphn1ij.at(id)*gamn1ij.at(id)+wind_vel[id].u*wind_vel[id].u*alphn1ij.at(id)* 
        gamn1ij.at(id)+ufvf*(alphn1ij.at(id)*gamn2ij.at(id)+alphn2ij.at(id)*gamn1ij.at(id))
        +wind_vel[id].u*wind_vel[id].v*(alphn1ij.at(id)*gamn2ij.at(id)+ alphn2ij.at(id)*gamn1ij.at(id))
        +ufwf*(alphn1ij.at(id)*gamn3ij.at(id)+alphn3ij.at(id)*gamn1ij.at(id))
        +wind_vel[id].u*wind_vel[id].w*(alphn1ij.at(id)*gamn3ij.at(id)+alphn3ij.at(id)*gamn1ij.at(id))+vfsq*alphn2ij.at(id) 
        *gamn2ij.at(id)+wind_vel[id].v*wind_vel[id].v*alphn2ij.at(id)*gamn2ij.at(id) 
        +vfwf*(alphn2ij.at(id)*gamn3ij.at(id)+alphn3ij.at(id)*gamn2ij.at(id)) 
        +wind_vel[id].v*wind_vel[id].w*(alphn2ij.at(id)*gamn3ij.at(id)+alphn3ij.at(id)*gamn2ij.at(id))+wfsq*alphn3ij.at(id)
        *gamn3ij.at(id)+wind_vel[id].w*wind_vel[id].w *alphn3ij.at(id)*gamn3ij.at(id);

    e11=alph1ij.at(id)*alphn1ij.at(id)+alph2ij.at(id)*betn1ij.at(id)+alph3ij.at(id)*gamn1ij.at(id);
    e12=alph1ij.at(id)*alphn2ij.at(id)+alph2ij.at(id)*betn2ij.at(id)+alph3ij.at(id)*gamn2ij.at(id);
    e13=alph1ij.at(id)*alphn3ij.at(id)+alph2ij.at(id)*betn3ij.at(id)+alph3ij.at(id)*gamn3ij.at(id);
    e21=bet1ij.at(id)*alphn1ij.at(id)+bet2ij.at(id)*betn1ij.at(id)+bet3ij.at(id)*gamn1ij.at(id);
    e22=bet1ij.at(id)*alphn2ij.at(id)+bet2ij.at(id)*betn2ij.at(id)+bet3ij.at(id)*gamn2ij.at(id);
    e23=bet1ij.at(id)*alphn3ij.at(id)+bet2ij.at(id)*betn3ij.at(id)+bet3ij.at(id)*gamn3ij.at(id);
    e31=gam1ij.at(id)*alphn1ij.at(id)+gam2ij.at(id)*betn1ij.at(id)+gam3ij.at(id)*gamn1ij.at(id);
    e32=gam1ij.at(id)*alphn2ij.at(id)+gam2ij.at(id)*betn2ij.at(id)+gam3ij.at(id)*gamn2ij.at(id);
    e33=gam1ij.at(id)*alphn3ij.at(id)+gam2ij.at(id)*betn3ij.at(id)+gam3ij.at(id)*gamn3ij.at(id);
    
}



void ParticleControl::rotate2d(const int &id,const float &cosphi,const float &sinphi,const float &upsqg,
                                   const float &upvpg,const float &vpsqg,const float &wpsqg,const float &upwpg,const float &vpwpg ){
    // this subroutine rotates variables from primed system aligned with the overall wind
    // into the regular grid system
    float ub,vb,wb,upb,vpb;
    ub=wind_vel[id].u;
    vb=wind_vel[id].v;
    wb=wind_vel[id].w;
    upb=ub*cosphi+vb*sinphi;
    vpb=-ub*sinphi+vb*cosphi;
    ufsqgi.at(id)=upsqg*cosphi*cosphi+upb*upb*cosphi*cosphi-2.f*cosphi*sinphi*upvpg 
        -2.f*cosphi*sinphi*upb*vpb+vpsqg*sinphi*sinphi+sinphi*sinphi*vpb*vpb-ub*ub;

    vfsqgi.at(id)=upsqg*sinphi*sinphi+upb*upb*sinphi*sinphi+2.f*upvpg*sinphi*cosphi 
        +2.f*upb*vpb*sinphi*cosphi+vpsqg*cosphi*cosphi+vpb*vpb*cosphi*cosphi-vb*vb;
    wfsqgi.at(id)=wpsqg;
    ufvfgi.at(id)=upsqg*cosphi*sinphi+upb*upb*cosphi*sinphi+upvpg
        *(cosphi*cosphi-sinphi*sinphi)+upb*vpb*(cosphi*cosphi-sinphi*sinphi)-vpsqg*sinphi*cosphi
        -vpb*vpb*sinphi*cosphi-ub*vb;
    ufwfgi.at(id)=upwpg*cosphi+upb*wb*cosphi-vpwpg*sinphi-vpb*wb*sinphi-ub*wb;
    vfwfgi.at(id)=upwpg*sinphi+upb*wb*sinphi+vpwpg*cosphi+vpb*wb*cosphi-vb*wb;
}

int ParticleControl::nint(const float t1){
    return int(t1 + .5f);
}
int ParticleControl::nint(const int t1){
    return int(t1 + .5f);
}

