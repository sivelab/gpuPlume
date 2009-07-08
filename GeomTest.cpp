#include <math.h>
#include "GeomTest.h"
#include "glErrorUtil.h"
#include <iostream>
#include <fstream>
#include <sstream>


#ifdef WIN32
#include <windows.h>
#include <stdio.h>
#include <conio.h>
#include <tchar.h>

#endif

//GLenum texType2;
//GLuint tex3d[2];
//GLuint case_to_numpoly[1];
//GLuint edge_connect_list[1];
int oneTime;
//GLuint query;
//int numP;
//int mesh;
bool doOnce;

GeomTest::GeomTest(Util* u){
  doOnce = true;

  util = u;

  //from util
  twidth = util->twidth;
  theight = util->theight;
  nx = util->nx;
  ny = util->ny;
  nz = util->nz;
  time_step = util->time_step;
  nxdx = (int)(nx*(1.0/util->dx));
  nydy = (int)(ny*(1.0/util->dy));
  nzdz = (int)(nz*(1.0/util->dz));
 
  texType = GL_TEXTURE_RECTANGLE_ARB;
  //texType2 = GL_TEXTURE_3D;//GL_TEXTURE_RECTANGLE_ARB;
  int_format = GL_RGBA32F_ARB;

  odd = true;
  dump_contents = true;
  oneTime = 0;
  //paused = false;

  //mesh = 4;
  firstTime = true;
  
}
GeomTest::~GeomTest(){}

void GeomTest::init(bool OSG){

 
  //glFrontFace(GL_CW);
 
  
  pc = new ParticleControl(texType, twidth,theight,nx,ny,nz,util->dx,util->dy,util->dz);
  pc->setUstarAndSigmas(util->ustar);
  pc->setBuildingParameters(util->numBuild,util->numSides,util->xfo,util->yfo,util->zfo,util->ht,util->wti,util->lti,util->gamma);
  pc->setQuicFilesPath(util->quicFilesPath);

  dc = new DisplayControl(nx,ny,nz, texType, util->dx,util->dy,util->dz);  
  dc->initVars(util->numBuild,util->numSides,util->xfo,util->yfo,util->zfo,util->ht,util->wti,util->lti,util->gamma);

  if(util->numBuild == 0){
    dc->draw_buildings = false;  
  }
  else{   
    dc->draw_buildings = true;
  }
  /*
  int num_vertices;
  if(GL_EXT_geometry_shader4)
    {
      std::cout << "Ready for geom shader!" << std::endl;
      //int numV = 0;
      //glGetIntegerv(GL_MAX_GEOMETRY_OUTPUT_VERTICES_EXT, &numV);
      
      glGetIntegerv(GL_MAX_GEOMETRY_OUTPUT_VERTICES_EXT, &num_vertices);
      std::cout << num_vertices << " number of output vertices." << std::endl;
    }
  
  geomShader.addShader("Shaders/isoSurface_vp.glsl", GLSLObject::VERTEX_SHADER);
  geomShader.addShader("Shaders/isoSurface_gp.glsl", GLSLObject::GEOMETRY_SHADER);
  //geomShader.addShader("Shaders/streamLines_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  geomShader.setInputandOutput(GL_POINTS,GL_TRIANGLE_STRIP,num_vertices);
  geomShader.createProgram();
  u_tau3D = geomShader.createUniform("tau");
  u_case = geomShader.createUniform("case_to_numpoly");
  u_edge = geomShader.createUniform("edge_connect_list");
  u_dx = geomShader.createUniform("dx");
  u_dy = geomShader.createUniform("dy");
  u_dz = geomShader.createUniform("dz");
  u_mesh = geomShader.createUniform("mesh");

  geomShader.activate();
  
  float scale = 1.0/(float)mesh;
  
  glUniform1fARB(u_dx,scale/(float)(nx));
  glUniform1fARB(u_dy,scale/(float)(ny));
  glUniform1fARB(u_dz,scale/(float)(nz));
  glUniform1fARB(u_mesh,scale);

  geomShader.deactivate();


  testShader.addShader("Shaders/render3D_vp.glsl", GLSLObject::VERTEX_SHADER);
  testShader.addShader("Shaders/render3D_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  testShader.createProgram();
  u_slice = testShader.createUniform("slice");
  uniform_tau = testShader.createUniform("tau");
  */
  
  //isoShader.addShader("Shaders/isoColor_vp.glsl", GLSLObject::VERTEX_SHADER);
  //isoShader.addShader("Shaders/isoColor_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  //isoShader.createProgram();
  


  glEnable(texType);
  glGenTextures(18, texid);
  
  //Textures used:
  positions0 = texid[0];
  positions1 = texid[1];
  windField = texid[3];
  randomValues = texid[4];
  prime0 = texid[5];
  prime1 = texid[6];
  lambda = texid[7];
  tau_dz = texid[8];
  duvw_dz = texid[9];

  //Texture for mean velocities
  meanVel0 = texid[10];
  meanVel1 = texid[11];
  //Texture used to hold current velocity
  currDirection = texid[12];
  //Texture used for building information
  buildings = texid[13];
  //Texture used for cell type information
  cellType = texid[14];
 
  //Texture to store path lines
  pathTex = texid[15];
  
  //Texture for tau11,tau22,tau33,and tau13
  tau = texid[16];

  /////////////////////////////
  setupTextures(); 
  /////////////////////////////
  

  //glDisable(texType);

  //isoSurface = new IsoSurface(util,pc);


  /*
  glEnable(texType2);
  glGenTextures(2,tex3d);


  //Texture to render to
  GLfloat* data = new GLfloat[(nx+1)*(ny+1)*(nz+1)*4];
  for(int k=0;k<=nz;k++){
    for(int i=0; i <= ny; i++){
      for(int j=0; j <= nx; j++){
	//int idx = k*ny*nx + i*nx + j;

	int tidx = k*ny*nx*4 + i*nx*4 + j*4;
	
	data[tidx] = 0.0;//pc->tau[idx].t11;
	data[tidx+1] = 0.0;//pc->tau[idx].t22;
	data[tidx+2] = 0.0;//pc->tau[idx].t33;
	data[tidx+3] = 1.0;//pc->tau[idx].t13;
      }
    }
  }
  


  glBindTexture(texType2, tex3d[0]);
  //glTexParameteri(texType2, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  //glTexParameteri(texType2, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(texType2, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(texType2, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(texType2, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(texType2, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(texType2, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
  glTexImage3D(texType2, 0, int_format, nx+1, ny+1, nz+1,0,GL_RGBA, GL_FLOAT,data); 
  
  delete [] data;

  GLfloat* data2 = new GLfloat[nx*ny*nz*4];

  //Texture to read from
  for(int k=0;k<nz;k++){
    for(int i=0; i < ny; i++){
      for(int j=0; j < nx; j++){
	int idx = k*ny*nx + i*nx + j;

	int tidx = k*ny*nx*4 + i*nx*4 + j*4;
	
	data2[tidx] = pc->tau[idx].t11;
	data2[tidx+1] = pc->tau[idx].t22;
	data2[tidx+2] = pc->tau[idx].t33;
	data2[tidx+3] = pc->tau[idx].t13;
      }
    }
  }
  
  glBindTexture(texType2, tex3d[1]);
  glTexParameteri(texType2, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(texType2, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  //glTexParameteri(texType2, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  //glTexParameteri(texType2, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(texType2, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(texType2, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(texType2, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
  glTexImage3D(texType2, 0, int_format, nx, ny, nz,0,GL_RGBA, GL_FLOAT,data2); 

  delete[] data2;

  glBindTexture(texType2,0);
  glDisable(texType2);
  */

  /*
  glGenBuffersARB(1, iso_surface_buffer);

 
  //for(int i=0; i < 1; i++){
    glBindBufferARB(GL_ARRAY_BUFFER_ARB,iso_surface_buffer[0]);
    //Each buffer is for the size of one layer
    //make a '4' = mesh
    glBufferDataARB(GL_ARRAY_BUFFER_ARB, (nz)*(nx)*(ny)*15*4*sizeof(GLfloat),0, GL_STREAM_COPY);
    //glBufferDataARB(GL_ARRAY_BUFFER_ARB, 100000*4*sizeof(GLfloat),0, GL_STREAM_COPY);
    
    //}
  glBindBuffer(GL_ARRAY_BUFFER,0);
  */
  /*
  int* attr = new int[1];
  attr[0] = GL_POSITION;

 
  
  geomShader.setVaryingOutput(1,attr,GL_INTERLEAVED_ATTRIBS_NV);
  //glTransformFeedbackAttribsNV(1,attr,GL_INTERLEAVED_ATTRIBS_NV);
  */

  initFBO();

  //Read in lookup tables and store in textures
  //readInTables();
  
  //glGenQueries(1,&query);

}
/*
void GeomTest::readInTables(){
  glEnable(GL_DEPTH_TEST);

  std::ifstream fp;
  std::string path;
  path = "case_to_numpolys";
  std::string path2 = "edge_connect_list";

  fp.open(path.c_str());
  if(!fp.is_open()){
    std::cerr << "Couldn't open file: case_to_numpolys";
    exit(1);
  }

  GLfloat* data = new GLfloat[256*4];
  std::string comma;

  for(int i = 0; i < 256; i++){
    int idx = i*4;
    fp >> data[idx];
    fp >> data[idx+1];
    fp >> data[idx+2];
    fp >> data[idx+3];

    //A comma follows each set 
    fp >> comma;
  }
  
  //for(int i = 0; i < 256; i++){
  //std::cout << data[i*4] << " " << 0 << " " << 0 << " " << 0 << " , ";
  //}

  fp.close();

  fp.open(path2.c_str());
  if(!fp.is_open()){
    std::cerr << "Couldn't open file: edge_connect_list";
    exit(1);
  }

  //read in edge list
  GLfloat* data2 = new GLfloat[256*4*5];
  for(int i=0; i < 256; i++){
    for(int j=0; j < 5; j++){
      int idx = i*5*4 + j*4;
      
      fp >> data2[idx];
      fp >> data2[idx+1];
      fp >> data2[idx+2];
      fp >> data2[idx+3];

      //A comma follows each set 
      fp >> comma;

    }
  }
  //for(int i=0; i < 256; i++){
    //for(int j=0; j < 5; j++){
      //int idx = i*5*4 + j*4;
      //std::cout << data2[idx] << " " << data2[idx+1] << " " << data2[idx+2] << " " << data2[idx+3] << " , ";
    //}
   // }


  fp.close();

  //Create case to numpoly texture
  glEnable(GL_TEXTURE_1D);
  glGenTextures(1,case_to_numpoly);
  
  glBindTexture(GL_TEXTURE_1D,case_to_numpoly[0]);

  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexImage1D(GL_TEXTURE_1D, 0, int_format, 256, 0,GL_RGBA, GL_FLOAT,data); 

  glBindTexture(GL_TEXTURE_1D,0);

  glDisable(GL_TEXTURE_1D);

  //Create edge connect list texture
  glEnable(GL_TEXTURE_2D);
  glGenTextures(1,edge_connect_list);
  
  glBindTexture(GL_TEXTURE_2D,edge_connect_list[0]);

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexImage2D(GL_TEXTURE_2D, 0, int_format, 5,256,0,GL_RGBA, GL_FLOAT,data2); 

  glDisable(GL_TEXTURE_2D);

  delete [] data;
  delete [] data2;

}
*/
int GeomTest::display(){
    

  
   
    glGetIntegerv(GL_DRAW_BUFFER, &draw_buffer);

    if(oneTime < 2){

      //render the 3D density function texture
      //isoSurface->render3DTexture(fbo);
    
    
      glDrawBuffer(draw_buffer);
      ///////////////////////////////////////////////////
    
      /////////////////////////////////////////////////////
      //Render geometry shader outputs to the vertex buffer
      /////////////////////////////////////////////////////
      // isoSurface->createIsoSurface(iso_surface_buffer[0]);
    
    

      oneTime++;
    
    }
  
  ////////////////////////////////////////////////////////

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  

  glViewport(0, 0, glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60.0, glutGet(GLUT_WINDOW_WIDTH)/float(glutGet(GLUT_WINDOW_HEIGHT)), 1.0, 250.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0.0,0.0,0.0,1.0);
	 
  dc->drawVisuals(vertex_buffer, duvw_dz, color_buffer, numInRow, twidth, theight, 0, 0);
     
  
  ///////////////////////////////////
  //Draw the vertex buffer
  /////////////////////////////////// 

  //isoSurface->draw(iso_surface_buffer[0]);
  
  

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(-1, 1, -1, 1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  

  //glDisable(texType2);
  CheckErrorsGL("END : visualization");
  firstTime = false;

  glutSwapBuffers();
  return 1;

}
void GeomTest::initFBO(void){
  // Get the Framebuffer Object ready
  fbo = new FramebufferObject();
  fbo->Bind();
   
  //Attach textures to framebuffer object
  fbo->AttachTexture(GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_3D, isoSurface->tex3d[0]);
 
  fbo->IsValid();
  FramebufferObject::Disable();


}
void GeomTest::setupTextures(){
  std::ofstream out;
  out.open(texFileName.c_str());

  CheckErrorsGL("BEGIN : Creating textures");

  int sz = 4;
  GLfloat *data = new GLfloat[ twidth * theight * sz];
		
  // Creates wind field data texture
  //The variable numInRow gets set here and should not be changed after being set!
  pc->initWindTex(windField, &numInRow, util->windFieldData);
  CheckErrorsGL("\tcreated texid[3], the wind field texture...");

  //Creates lambda, tau/dz, and duvw/dz textures
  if(util->windFieldData < 5)
    pc->initLambda_and_TauTex(lambda, tau_dz, duvw_dz);
  else if(util->windFieldData == 5)
    pc->initLambda_and_TauTex_fromQUICFILES(windField, lambda, tau_dz, duvw_dz, tau);
  else
    pc->initLambda_and_Taus_withCalculations(windField, lambda, tau_dz, duvw_dz, tau);
  CheckErrorsGL("\tcreated texid[7], the lambda texture...");

  //Print out max and min taus
  /*std::cout << "Max and Min Taus" << std::endl;
  std::cout << pc->tauMax[0] << " " << pc->tauMax[1] << " " << pc->tauMax[2] << " " << 
    pc->tauMax[3] << std::endl;
  
  std::cout << pc->tauMin[0] << " " << pc->tauMin[1] << " " << pc->tauMin[2] << " " << 
  pc->tauMin[3] << std::endl;*/

  pc->addBuildingsInWindField(cellType);
  CheckErrorsGL("\tcreated texid[14], the cell type texture...");

  for (int j=0; j<theight; j++)
    for (int i=0; i<twidth; i++)
      {
	int idx = j*twidth*sz + i*sz;
	data[idx] = 100.0;
	data[idx+1] = 100.0;
	data[idx+2] = 100.0;
	data[idx+3] = lifeTime+1;
      }
  
  // create the base texture with inital vertex positions
  pc->createTexture(texid[0], int_format, twidth, theight, data);
  CheckErrorsGL("\tcreated texid[0], the position texture...");

  // create a second texture to double buffer the vertex positions
  pc->createTexture(texid[1], int_format, twidth, theight, data);
  CheckErrorsGL("\tcreated texid[1], the position texture (double buffer)...");
  

  CheckErrorsGL("END : Creating textures");
 
}
