//////////////////////////////////////////////////
// This class takes care of the represenation of
// the particle positions and wind field on the GPU.
// It also loads the Quic PLume Fortran References.
//////////////////////////////////////////////////

#include <GL/glew.h>
#include "framebufferObject.h"
#include "GLSL.h"

class ParticleControl{

 public:
 
  ParticleControl(GLenum,int,int,int,int,int,double*,double*,double*);

  void setupAdvectShader(float*,int*, float);

  void setupPrimeShader();

  //This function puts the values held in the variable, data, into a 2D texture 
  //on the GPU. 
  void createTexture(GLuint texId, GLenum format,  int w, int h, GLfloat* data); 
  void createWrappedTexture(GLuint texId, GLenum format,  int w, int h, GLfloat* data); 
 
  //This function maps the 3D wind field into a 2D texture that can be
  //stored on the GPU.
  void initWindTex(GLuint,GLuint, int* numInRow, int dataSet);

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
  void advect(FramebufferObject*,bool,GLuint,GLuint,GLuint,GLuint,GLuint,GLuint,float);

  void updatePrime(FramebufferObject*,bool,GLuint,GLuint,GLuint,GLuint,GLuint,float);

  void getDomain(int* , int*, int*);

  //Sets the ustar and sigma values which can now be used to create
  //the lambda texture and CoE/2 values. 
  void setUstarAndSigmas(float);
  
 private:

  void test1();
  void randomWindField();
  void quicPlumeWindField();  
  void uniformUWindField();
  void gravity();
  
  typedef struct{
    float u;
    float v;
    float w;
  }wind;
  wind* data3d;

  double* u_quicPlumeData;
  double* v_quicPlumeData;
  double* w_quicPlumeData;

  int nx;
  int ny;
  int nz;

  int twidth, theight;

  float ustar,sigU,sigV,sigW;

  GLSLObject init_shader, pass1_shader, prime_shader;
  //Variables for prime shader
  GLint uniform_prime, uniform_windTex, uniform_random;
  GLint uniform_dt,uniform_lambda;
  //Variables for advect shader
  GLint uniform_postex, uniform_wind, uniform_randomTexture;
  GLint uniform_primePrev, uniform_primeCurr, uniform_timeStep;

  GLenum texType;
  GLfloat* buffer_mem;
};
