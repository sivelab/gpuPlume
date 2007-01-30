#include <GL/glew.h>
#include <GL/glut.h>
#include <iostream>
#include "framebufferObject.h"
#include "GLSL.h"

class ParticleControl{

 public:
 
  ParticleControl(GLenum);

  void createTexture(GLuint texId, GLenum format,  int w, int h, GLfloat* data); 
 
  void initWindTex(GLuint texId, int* numInRow);

  void initParticlePositions(FramebufferObject*, int, int, GLSLObject, GLuint);

  void dumpContents(int w, int h);

  void getDomain(int* , int*, int*);
  void test1();
  void test2();
  void test3();  
  
 private:
  
  typedef struct{
    float u;
    float v;
    float w;
  }wind;
  wind* data3d;

  int nx;
  int ny;
  int nz;

  GLenum texType;
  GLfloat* buffer_mem;
};
