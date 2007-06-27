#ifndef  __PC_H__
#define  __PC_H__

//////////////////////////////////////////////////
// This class takes care of the represenation of
// the particle positions and wind field on the GPU.
// It also loads the Quic PLume Fortran References.
//////////////////////////////////////////////////

#include <GL/glew.h>
#include "framebufferObject.h"
#include "GLSL.h"
#include <string>

class ParticleControl{

 public:
 
  ParticleControl(GLenum,int,int,int,int,int,double*,double*,double*);

  void setupAdvectShader(int, float);

  void setupPrimeShader(int); //Included argument -- Balli(04/12/07)

  void setupPrime_and_AdvectShader(int,float);

  void setupNonGaussianShader(int,float);

  void nonGaussianAdvect(FramebufferObject*,bool,GLuint,GLuint,GLuint,
			 GLuint,GLuint,GLuint,GLuint,GLuint,GLuint,float);

  //Performs update prime and advect with one shader using multiple render targets.
  void updatePrimeAndAdvect(FramebufferObject*,bool,GLuint,GLuint,GLuint,
			    GLuint,GLuint,GLuint,GLuint,float);

  //This function puts the values held in the variable, data, into a 2D texture 
  //on the GPU. 
  void createTexture(GLuint texId, GLenum format,  int w, int h, GLfloat* data); 
  void createWrappedTexture(GLuint texId, GLenum format,  int w, int h, GLfloat* data); 
 
  //This function maps the 3D wind field into a 2D texture that can be
  //stored on the GPU.
  void initWindTex(GLuint, int* numInRow, int dataSet);

  //Call this function for uniform wind field.
  void initLambdaTex(GLuint, int);
  
  void initLambda_and_TauTex(GLuint,GLuint,GLuint,int);

  //This function is used to initialize the particle positions.
  //It was needed when we weren't able to directly put 32-bit floating point
  //values directly into a texture. The help of a shader program was needed
  //to do that.  Thanks to a driver update, we don't have to do this anymore.
  void initParticlePositions(FramebufferObject*, GLuint);

  //This will output the values of the current texture being read.
  void dumpContents();

  //Writes the texture to a ppm image.
  void writePPM(const std::string&);
  short c2Short(float);

  //Advects particle positions using the windfield.
  //First GLuint is the windfield texture.
  //Second and third GLuint are the two position textures. 
  void advect(FramebufferObject*,bool,GLuint,GLuint,GLuint,GLuint,GLuint,float);

  void updatePrime(FramebufferObject*,bool,GLuint,GLuint,GLuint,GLuint,GLuint,GLuint,GLuint,float);
  // included two more arguments in the above function for position textures. --Balli(04/12/07)
  void getDomain(int* , int*, int*);

  //Sets the ustar and sigma values which can now be used to create
  //the lambda texture and CoE/2 values. 
  void setUstarAndSigmas(float);
  
  void createPosImages(bool);
  void createPrimeImages(bool);

  bool outputPrime;
   
  typedef struct{
    float u;
    float v;
    float w;
  }wind;

  wind* sig;

 private:
  float min,max;

  void test1();
  void randomWindField();
  void quicPlumeWindField();  
  void uniformUWindField();
  void variedUWindField();

  void printPrime(bool,bool);
 
  wind* wind_vel;

  typedef struct{
    float t11;
    float t22;
    float t33;
    float t13;
  }Matrix;
  Matrix* tau;

  double* u_quicPlumeData;
  double* v_quicPlumeData;
  double* w_quicPlumeData;

  int nx;
  int ny;
  int nz;

  //Number of Particles
  int twidth, theight;

  //Width and height of wind,lambda textures
  int width, height;

  float ustar,sigU,sigV,sigW;

  GLSLObject init_shader, pass1_shader, prime_shader, mrt_shader;
  GLSLObject nonGaussian_shader;

  //Variables for prime shader
  GLint uniform_prime, uniform_windTex, uniform_random,uniform_pos;
  GLint uniform_dt,uniform_lambda;
  GLint uniform_randomTexCoordOffset, uniform_randomTexWidth, uniform_randomTexHeight;

  //Variables for advect shader
  GLint uniform_postex, uniform_wind, uniform_randomTexture;
  GLint uniform_primePrev, uniform_primeCurr, uniform_timeStep;

  GLint uniform_tau_dz, uniform_duvw_dz;

  GLenum texType;
  GLfloat* buffer_mem;

  GLint currentbuffer;

};
#endif //__PC_H__
