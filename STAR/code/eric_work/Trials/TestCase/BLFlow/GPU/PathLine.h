#ifndef __PATHLINE_H__
#define __PATHLINE_H__


#include <iostream>
#include <list>
#include <vector>
#include "framebufferObject.h"
#include "particleEmitter.h"

typedef struct{
  //Index into position texture
  int x;
  int y;
  //Index into pathLine texture
  int s;
  int t;
}pathIndex;

class PathLine{
 public:

  PathLine(int,int,GLenum);

  void setupGeomShader();

  void addNewPath(ParticleEmitter*);

  void updatePathLines(GLuint,GLuint,bool);

  void updateVBOS();

  void draw();

  void printPathLineTexture();

  bool doUpdate();

  void clear();

  int pathNum;
  bool startStream;

 private:
  //FramebufferObject* pathFbo;

  std::list<pathIndex> pathList;
  
  std::list<pathIndex>::iterator pathIter;

  int pwidth;
  int pheight;

  GLfloat* buffer;

  bool update;
  GLuint* path_buffer;
  GLint uniform_pos;

  GLSLObject pathLineShader;

  GLenum texType;
  //GLuint texId[1];

};

#endif //__PATHLINE_H__
