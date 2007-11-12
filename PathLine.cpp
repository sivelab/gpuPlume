#include "PathLine.h"
#include <math.h>

PathLine::PathLine(int w, int h, GLenum type){
  pwidth = w;
  pheight = h;
  texType = type;

  pathNum = 0;
  startStream = false;
  update = true;

  path_buffer = new GLuint[pheight];

  //Set up vbos for the paths
  glGenBuffersARB(pheight, path_buffer);
 
  for(int i=0; i < pheight; i++){

    glBindBufferARB(GL_ARRAY_BUFFER_ARB, path_buffer[i]);
    glBufferDataARB(GL_ARRAY_BUFFER_ARB, pwidth*4*sizeof(GLfloat), 0, GL_STREAM_COPY);
  }

  glBindBuffer(GL_ARRAY_BUFFER, 0);

  //Set up shader to update path lines
  pathLineShader.addShader("Shaders/pathLine_vp.glsl", GLSLObject::VERTEX_SHADER);
  pathLineShader.addShader("Shaders/pathLine_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  pathLineShader.createProgram();
  uniform_pos = pathLineShader.createUniform("positions");

  /*
  //Setup path line texture
  glEnable(texType);
  glGenTextures(1,texId);
  GLfloat *pdata = new GLfloat[ pwidth * pheight * 4];
  for (int j=0; j<pheight; j++)
    for (int i=0; i<pwidth; i++)
      {
	int idx = j*pwidth*4 + i*4;
	pdata[idx] = 0.0;
	pdata[idx+1] = 0.0;
	pdata[idx+2] = 0.0;
	pdata[idx+3] = 1.0;
      }
  
  glBindTexture(texType, texId[0]);
  glTexParameteri(texType, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(texType, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(texType, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(texType, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  glTexImage2D(texType, 0, GL_RGBA32F_ARB, pwidth, pheight, 0, GL_RGBA, GL_FLOAT, pdata);

  glBindTexture(texType, 0);
  delete [] pdata;

  //Setup fbo
  pathFbo = new FramebufferObject();
  pathFbo->Bind();
  pathFbo->AttachTexture(GL_COLOR_ATTACHMENT0_EXT, texType, texId[0]);
  pathFbo->IsValid();
  FramebufferObject::Disable();
  */
}
void PathLine::setupGeomShader(){
  /*streamLineShader.addShader("Shaders/streamLines_vp.glsl", GLSLObject::VERTEX_SHADER);
  streamLineShader.addShader("Shaders/streamLines_gp.glsl", GLSLObject::GEOMETRY_SHADER);
  streamLineShader.addShader("Shaders/streamLines_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  streamLineShader.setInputandOutput(GL_LINES,GL_LINES);
  streamLineShader.createProgram();*/

}
void PathLine::addNewPath(ParticleEmitter* pe){
  float x,y,z;
  int xindex,yindex; 
  pathIndex pIndex;

  pe->getReleasedPosition(&x,&y,&z);

  //pheight is how many path lines that can be generated
  if(pathNum < pheight){
    pe->getIndex(&xindex,&yindex);
    //Need to Punch Hole method to put initial starting position into 
    //path line texture

    //std::cout << "Index emitted: " << xindex << " "<< yindex << std::endl;

    pIndex.x = xindex;
    pIndex.y = yindex;
    pIndex.s = 0;
    pIndex.t = pathNum;

    pathList.push_back(pIndex);
     
    pathNum++;
    startStream = true;
    update = true;
  }

}

void PathLine::updatePathLines(GLuint positions0, GLuint positions1, bool odd){
  //pathFbo->Bind();
  //Need to bind the previous fbo after this call is made!!!

  glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0, pwidth, 0, pheight);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity(); 
    

  glActiveTexture(GL_TEXTURE0);
  if (odd)
    glBindTexture(texType, positions0);  // read from texture 0
  else 
    glBindTexture(texType, positions1);  // read from texture 1
 
  glUniform1iARB(uniform_pos, 0);
    
  pathLineShader.activate();
  pathIndex pIndex;

  pathIter = pathList.begin();
  while(pathIter != pathList.end()){

    pIndex = *pathIter;
    pathIndex &pIdx = *pathIter;

    //Length of path line is bound to pwidth
    if(pIndex.s < pwidth){
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
    }

    pathIter++;
  }
  pathLineShader.deactivate();

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(-1, 1, -1, 1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
      
  update = true;
  
  //updateVBOS();
}
void PathLine::updateVBOS(){
  //pathFbo->Bind();
  glReadBuffer(GL_COLOR_ATTACHMENT0_EXT); 

  pathIndex pIndex;
  pathIter = pathList.begin();
  int i=0;
  while(pathIter != pathList.end()){

    pIndex = *pathIter;
    glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB, path_buffer[i]);
    glReadPixels(0, pIndex.t, pIndex.s, 1, GL_RGBA, GL_FLOAT, 0);

    pathIter++;
    i++;
  }

}

void PathLine::draw(){

  glDisable(GL_COLOR_ARRAY);
  glEnableClientState(GL_VERTEX_ARRAY);
  int i=0;

  //glLineWidth(2.0);

  pathIndex pIndex;
  pathIter = pathList.begin();
  while(pathIter != pathList.end()){

    glBindBufferARB(GL_ARRAY_BUFFER, path_buffer[i]);
    glVertexPointer(4,GL_FLOAT,0,0);
     
    glColor4f(1.0,1.0,1.0,1.0);
 
    pIndex = *pathIter;
    glDrawArrays(GL_LINE_STRIP,0,pIndex.s);


    pathIter++;
    i++;
  }
  

}
bool PathLine::doUpdate(){
  if(pathList.empty())
    return false;
  /*
  streamIndex sIndex;
  if(update){
    indexIter = indexList.begin();
    while(indexIter != indexList.end()){
      sIndex = *indexIter;
      if(!sIndex.done){
	return update;
      }
    
      indexIter++;
    }

    }*/
  //update = false;
  return update;

}
void PathLine::clear(){
  pathList.clear();
  startStream = true;
  update = true;
  pathNum = 0;

}
void PathLine::printPathLineTexture(){
  //pathFbo->Bind();
  glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);

  buffer = new GLfloat[ pwidth * pheight * 4 ];  
  glReadPixels(0, 0, pwidth, pheight, GL_RGBA, GL_FLOAT, buffer);
  
  int pn =0;
  for (int j=0; j<pheight; j++){
    std::cout << pn << " ";
    for (int i=0; i<4; i++){
	  
       int idx = j*pwidth*4 + i*4;
       
       std::cout << "( " << buffer[idx] << " ,";
       std::cout << buffer[idx+1] << " ,";
       std::cout << buffer[idx+2] << " ,";
       std::cout << buffer[idx+3] << " )";
      
    }
    pn++;
    std::cout << "\n\n";
  }
  delete [] buffer;

}
