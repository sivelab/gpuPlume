#include <math.h>
#include "GeomTest.h"
#include "glErrorUtil.h"


#ifdef WIN32
#include <windows.h>
#include <stdio.h>
#include <conio.h>
#include <tchar.h>

#endif

GeomTest::GeomTest(Util* u){
  
  util = u;

  pwidth = 100;
  pheight = 10;

  //from util
  twidth = util->twidth;
  theight = util->theight;
  nx = util->nx;
  ny = util->ny;
  nz = util->nz;
  time_step = util->time_step;

  odd = true;
  dump_contents = true;
  //paused = false;

  texType = GL_TEXTURE_RECTANGLE_ARB;
  int_format = GL_RGBA32F_ARB;
  int_format_init = GL_RGBA;

}
GeomTest::~GeomTest(){}

void GeomTest::init(bool OSG){

  //pc = new ParticleControl(texType, twidth,theight,nx,ny,nz);
 

  glEnable(texType);
  glGenTextures(3, texid);
  positions0 = texid[0];
  positions1 = texid[1];
  paths = texid[2];

  setupTextures();

  glGenBuffersARB(4, vbo_buffer);
  vertex_buffer = vbo_buffer[0];
  path_buffer = vbo_buffer[1];

  glBindBufferARB(GL_ARRAY_BUFFER_ARB, vertex_buffer);
  glBufferDataARB(GL_ARRAY_BUFFER_ARB, twidth*theight*4*sizeof(GLfloat), 0, GL_STREAM_COPY);

  for(int i=1; i < 4; i++){

    glBindBufferARB(GL_ARRAY_BUFFER_ARB, vbo_buffer[i]);
    glBufferDataARB(GL_ARRAY_BUFFER_ARB, pwidth*4*sizeof(GLfloat), 0, GL_STREAM_COPY);
  }

  glBindBuffer(GL_ARRAY_BUFFER, 0);

  CheckErrorsGL("before fbo");

  initFBO();

  testAdvectShader.addShader("Shaders/test_vp.glsl", GLSLObject::VERTEX_SHADER);
  testAdvectShader.addShader("Shaders/test_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  testAdvectShader.createProgram();
  uniform_postPP = testAdvectShader.createUniform("positions");

  pathLineShader.addShader("Shaders/pathLine_vp.glsl", GLSLObject::VERTEX_SHADER);
  pathLineShader.addShader("Shaders/pathLine_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  pathLineShader.createProgram();
  uniform_posPP = pathLineShader.createUniform("positions");
  uniform_x = pathLineShader.createUniform("x");
  uniform_y = pathLineShader.createUniform("y");

  /*
  int num_vertices;
  if(GL_EXT_geometry_shader4)
    {
      glGetIntegerv(GL_MAX_GEOMETRY_OUTPUT_VERTICES_EXT, &num_vertices);
    }
  geomShader.addShader("Shaders/streamLines_vp.glsl", GLSLObject::VERTEX_SHADER);
  geomShader.addShader("Shaders/streamLines_gp.glsl", GLSLObject::GEOMETRY_SHADER);
  //geomShader.addShader("Shaders/streamLines_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  geomShader.setInputandOutput(GL_POINTS,GL_POINTS,num_vertices);
  geomShader.createProgram();
  */

  pathIdx pIdx;
  pathIdx p2Idx;
  pathIdx p3Idx;
  pIdx.x = 0;
  pIdx.y = 0;
  pIdx.s = 0;
  pIdx.t = 0;

  p2Idx.x = 1;
  p2Idx.y = 1;
  p2Idx.s = 0;
  p2Idx.t = 1;

  p3Idx.x = 1;
  p3Idx.y = 0;
  p3Idx.s = 0;
  p3Idx.t = 2;

  pathList.push_back(pIdx);
  pathList.push_back(p2Idx);
  pathList.push_back(p3Idx);


}

int GeomTest::display(){
  glGetIntegerv(GL_DRAW_BUFFER, &draw_buffer);
  glEnable(texType);
  // bind the framebuffer object so we can render to the 2nd texture
  fbo->Bind();

  if(!paused || !inPauseMode){

     //////////////////////////////////////////////////////////////////
    //Update path texture
    //////////////////////////////////////////////////////////////////
    pathFbo->Bind();
   
    glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, pwidth, 0, pheight);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity(); 
    
    pathLineShader.activate();
    glActiveTexture(GL_TEXTURE0);
    if (odd)
      glBindTexture(texType, positions0);  // read from texture 0
    else 
      glBindTexture(texType, positions1);  // read from texture 1
 
    glUniform1iARB(uniform_posPP, 0);
    
    
    pathIdx pIndex;

    pathIter = pathList.begin();
    while(pathIter != pathList.end()){

      pIndex = *pathIter;
      pathIdx &pIdx = *pathIter;


      glUniform1fARB(uniform_x,pIndex.x);
      glUniform1fARB(uniform_y,pIndex.y);

      //set viewport with s and t
      glViewport(pIndex.s,pIndex.t,1,1);

     
      //Punch Hole new point in path line into path line texture
      glBegin(GL_POINTS);
      {
	//pass texture coordinates into position texture to shader
	//using the color attribute variable
	glColor4f(pIndex.x,pIndex.y,0.0,1.0);
	glVertex2f(0.5,0.5);
      }
      glEnd();
      

      //update s coordinate into path line texture
      pIdx.s++;

      pathIter++;
    }
    pathLineShader.deactivate();

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(-1, 1, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    //////////////////////////////////////////////////////////////////
    // Now do advection
    //////////////////////////////////////////////////////////////////
    fbo->Bind();

    if (odd)
      glDrawBuffer(GL_COLOR_ATTACHMENT1_EXT);
    else 
      glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glViewport(0, 0, twidth, theight);

    testAdvectShader.activate();
    
    glActiveTexture(GL_TEXTURE0);
    if (odd)
      glBindTexture(texType, positions0);  // read from texture 0
    else 
      glBindTexture(texType, positions1);  // read from texture 1
 
    glUniform1iARB(uniform_postPP, 0);

    glBegin(GL_QUADS);
    {
      glTexCoord2f(0, 0);				glVertex3f(-1, -1, -0.5f);
      glTexCoord2f(float(twidth), 0);			glVertex3f( 1, -1, -0.5f);
      glTexCoord2f(float(twidth), float(theight));	glVertex3f( 1,  1, -0.5f);
      glTexCoord2f(0, float(theight));			glVertex3f(-1,  1, -0.5f);
    }
    glEnd();
  
    testAdvectShader.deactivate();
 
    glBindTexture(texType, 0);
    //////////////////////////////////////////////////////////////////////
    if (dump_contents)
    {
      if(odd)
	glReadBuffer(GL_COLOR_ATTACHMENT1_EXT);
      else
	glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);

      GLfloat* buffer_mem = new GLfloat[ twidth * theight * 4 ];  
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
      //FramebufferObject::Disable();
      
      pathFbo->Bind();
      glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);

      GLfloat* buffer = new GLfloat[ pwidth * pheight * 4 ];  
      glReadPixels(0, 0, pwidth, pheight, GL_RGBA, GL_FLOAT, buffer);
      pn = 0;
      for (int j=0; j<3; j++){
	std::cout << pn << " : ";
	for (int i=0; i<10; i++){
	  
	  int idx = j*pwidth*4 + i*4;
	  std::cout << "( " << buffer[idx] << "  ,";
	  std::cout << buffer[idx+1] << " ,";
	  std::cout << buffer[idx+2] << " ,";
	  std::cout << buffer[idx+3] << " )";
	}
	pn++;
	std::cout << "\n\n";
      }
      delete [] buffer;

      //FramebufferObject::Disable();
      fbo->Bind();

      dump_contents = false;
    }

    odd = !odd;
    paused = true;

  
  }
  
  glGetIntegerv(GL_READ_BUFFER, &read_buffer);
  if(odd){
    glReadBuffer(GL_COLOR_ATTACHMENT1_EXT);
  }
  else 
    glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);

  glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB, vertex_buffer);
  glReadPixels(0, 0, twidth, theight, GL_RGBA, GL_FLOAT, 0);
 
  ////////////////////////////////////////////////////////////
  //Get path lines into vbos
  ////////////////////////////////////////////////////////////
  
  pathFbo->Bind();
  
  glReadBuffer(GL_COLOR_ATTACHMENT0_EXT); 

  pathIdx pIndex;
  pathIter = pathList.begin();
  int i=1;
  while(pathIter != pathList.end()){

    pIndex = *pathIter;
    glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB, vbo_buffer[i]);
    glReadPixels(0, pIndex.t, pIndex.s, 1, GL_RGBA, GL_FLOAT, 0);

    pathIter++;
    i++;
  }
  
  //////////////////////////////////////////////////////////////

  glBindBuffer(GL_PIXEL_PACK_BUFFER_ARB, 0);
  CheckErrorsGL("after glReadPixels");
  glReadBuffer(read_buffer);
  
  
  // Disable the framebuffer object
  FramebufferObject::Disable();
  glDrawBuffer(draw_buffer); // send it to the original buffer
  CheckErrorsGL("END : after 2nd pass");

  

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glViewport(0, 0, glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60.0, glutGet(GLUT_WINDOW_WIDTH)/float(glutGet(GLUT_WINDOW_HEIGHT)), 1.0, 250.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();	
  //glClearColor(util->bcolor[0],util->bcolor[1],util->bcolor[2],1.0);	
  glClearColor(0.0,0.0,0.0,1.0);	

  gluLookAt(0,10,0,0,0,0, 0, 0, 1 );
  
  
  
  glEnableClientState(GL_VERTEX_ARRAY);
  i=1;
  pathIter = pathList.begin();
  while(pathIter != pathList.end()){

    glBindBufferARB(GL_ARRAY_BUFFER, vbo_buffer[i]);
    glVertexPointer(4,GL_FLOAT,0,0);
     
    glColor4f(1.0,1.0,1.0,1.0);
 
    pIndex = *pathIter;
    glDrawArrays(GL_LINE_STRIP,0,pIndex.s);


    pathIter++;
    i++;
  }
  //geomShader.activate();

  glPointSize(4.0);
  glBegin(GL_POINTS);

  glColor4f(0.0,0.0,1.0,1.0);
  glVertex3f(0.0,0.0,0.0);

  glEnd();

  //geomShader.deactivate();

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(-1, 1, -1, 1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
     
  glDisable(texType);
  CheckErrorsGL("END : visualization");


  glutSwapBuffers();
  return 1;

}
void GeomTest::initFBO(void){
  // Get the Framebuffer Object ready
  fbo = new FramebufferObject();
  pathFbo = new FramebufferObject();
  fbo->Bind();
      
  glGetIntegerv(GL_MAX_COLOR_ATTACHMENTS_EXT, (GLint*)&maxColorAttachments);
  std::cout << "Max color attachments: " << maxColorAttachments << std::endl;
  
  //Attach textures to framebuffer object
  fbo->AttachTexture(GL_COLOR_ATTACHMENT0_EXT, texType, texid[0]); 
  fbo->AttachTexture(GL_COLOR_ATTACHMENT1_EXT, texType, texid[1]);
  
  fbo->IsValid();
  //FramebufferObject::Disable();
  
  pathFbo->Bind();

  pathFbo->AttachTexture(GL_COLOR_ATTACHMENT0_EXT, texType, texid[2]);

  pathFbo->IsValid();
  FramebufferObject::Disable();

}

void GeomTest::setupTextures(){
  
  CheckErrorsGL("BEGIN : Creating textures");

  int sz = 4;
  GLfloat *data = new GLfloat[ twidth * theight * sz];
		
  for (int j=0; j<theight; j++)
    for (int i=0; i<twidth; i++)
      {
	int idx = j*twidth*sz + i*sz;
	data[idx] = 0.0;
	data[idx+1] = 0.0;
	data[idx+2] = j*twidth + i;
	data[idx+3] = 1.0;
      }
  
  glBindTexture(texType, texid[0]);
  glTexParameteri(texType, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(texType, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(texType, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(texType, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  glTexImage2D(texType, 0, int_format, twidth, theight, 0, GL_RGBA, GL_FLOAT, data);
  
  // create the base texture with inital vertex positions
  //pc->createTexture(texid[0], int_format, twidth, theight, data);
  CheckErrorsGL("\tcreated texid[0], the position texture...");

  glBindTexture(texType, texid[1]);
  glTexParameteri(texType, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(texType, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(texType, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(texType, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  glTexImage2D(texType, 0, int_format, twidth, theight, 0, GL_RGBA, GL_FLOAT, data);
  
  // create a second texture to double buffer the vertex positions
  //pc->createTexture(texid[1], int_format, twidth, theight, data);
  CheckErrorsGL("\tcreated texid[1], the position texture (double buffer)...");


  GLfloat *pdata = new GLfloat[ pwidth * pheight * sz];
  for (int j=0; j<pheight; j++)
    for (int i=0; i<pwidth; i++)
      {
	int idx = j*pwidth*sz + i*sz;
	pdata[idx] = 0.0;
	pdata[idx+1] = 0.0;
	pdata[idx+2] = 0.0;
	pdata[idx+3] = 1.0;
      }
  //Create path Line texture
  glBindTexture(texType, texid[2]);
  glTexParameteri(texType, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(texType, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(texType, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(texType, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  glTexImage2D(texType, 0, int_format, pwidth, pheight, 0, GL_RGBA, GL_FLOAT, pdata);
  
  delete [] data;
  delete [] pdata;

  CheckErrorsGL("END : Creating textures");
 
}

