//////////////////////////////////////////////////
// This class takes care of the represenation of
// the particle positions and wind field on the GPU.
// It also loads the Quic PLume Fortran References.
//////////////////////////////////////////////////

#include <GL/glew.h>
#include <GL/glut.h>
#include <iostream>
#include "framebufferObject.h"
#include "GLSL.h"

class ParticleControl{

 public:
 
  ParticleControl(GLenum,int,int);

  void setupAdvectShader(float*,int*);

  //This function puts the values held in the variable, data, into a 2D texture 
  //on the GPU. 
  void createTexture(GLuint texId, GLenum format,  int w, int h, GLfloat* data); 
  void createWrappedTexture(GLuint texId, GLenum format,  int w, int h, GLfloat* data); 
 
  //This function maps the 3D wind field into a 2D texture that can be
  //stored on the GPU.
  void initWindTex(GLuint texId, int* numInRow);

  //This function is used to initialize the particle positions.
  //It was needed when we weren't able to directly put 32-bit floating point
  //values directly into a texture. The help of a shader program was needed
  //to do that.  Thanks to a driver update, we don't have to do this anymore.
  void initParticlePositions(FramebufferObject*, GLuint);

  //This will output the values of the current texture being read.
  void dumpContents();

  //Advects particle positions using the windfield.
  //First GLuint is the windfield texture.
  //Second and third GLuint are the two position textures. 
  void advect(FramebufferObject*,bool,GLuint,GLuint,GLuint,GLuint);

  void getDomain(int* , int*, int*);
  
  
 private:

  void test1();
  void test2();
  void test3();  
  
  typedef struct{
    float u;
    float v;
    float w;
  }wind;
  wind* data3d;

  int nx;
  int ny;
  int nz;

  int twidth, theight;

  GLSLObject init_shader, pass1_shader;
  GLint uniform_postex, uniform_wind, uniform_randomTexture, uniform_timeStep;

  GLenum texType;
  GLfloat* buffer_mem;
};
