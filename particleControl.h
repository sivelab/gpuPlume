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

typedef struct{
  float t11;
  float t22;
  float t33;
  float t13;
}Matrix;
  

class ParticleControl{

 public:
 
  ParticleControl(GLenum,int,int,int,int,int);

  void setupAdvectShader(int, float);

  void setupPrimeShader(int); //Included argument -- Balli(04/12/07)

  void setupPrime_and_AdvectShader(int,float);

  void setupNonGaussianShader(int,float);
  
  void setupReflectionShader(int,float);

  void setupMultipleBuildingsShader(int,float,int);

  void multipleBuildingsAdvect(bool,GLuint,GLuint,GLuint,
			 GLuint,GLuint,GLuint,GLuint,GLuint,GLuint,float,GLuint,GLuint);

  void nonGaussianAdvect(bool,GLuint,GLuint,GLuint,
			 GLuint,GLuint,GLuint,GLuint,GLuint,GLuint,float);

  //Performs update prime and advect with one shader using multiple render targets.
  void updatePrimeAndAdvect(bool,GLuint,GLuint,GLuint,
			    GLuint,GLuint,GLuint,GLuint,float);

  void reflectionAdvect(bool,GLuint,GLuint,GLuint,
			 GLuint,GLuint,GLuint,GLuint,GLuint,GLuint,float,float*);

  //This function puts the values held in the variable, data, into a 2D texture 
  //on the GPU. 
  void createTexture(GLuint texId, GLenum format,  int w, int h, GLfloat* data); 
  void createIntTexture(GLuint texId, GLenum format,  int w, int h, GLint* data); 
  void createWrappedTexture(GLuint texId, GLenum format,  int w, int h, GLfloat* data); 
 
  //This function maps the 3D wind field into a 2D texture that can be
  //stored on the GPU.
  void initWindTex(GLuint, int* numInRow, int dataSet);

  //Call this function for uniform wind field.
  void initLambdaTex(GLuint, int);
  
  void initLambda_and_TauTex(GLuint,GLuint,GLuint,int);
  void initLambda_and_TauTex_fromQUICFILES(GLuint,GLuint,GLuint,GLuint,GLuint,int);
  void initLambda_and_Taus_withCalculations(GLuint,GLuint,GLuint,GLuint,GLuint,int);


  //This function is used to initialize the particle positions.
  //It was needed when we weren't able to directly put 32-bit floating point
  //values directly into a texture. The help of a shader program was needed
  //to do that.  Thanks to a driver update, we don't have to do this anymore.
  void initParticlePositions(FramebufferObject*, GLuint);

  //This will output the values of the current texture being read.
  void dumpContents();

  //Prints out the previous and updated position values
  void printPositions(bool);

  //Prints out the mean velocity values
  void printMeanVelocities(bool);

  //Writes the texture to a ppm image.
  void writePPM(const std::string&);
  short c2Short(float);

  //Advects particle positions using the windfield.
  //First GLuint is the windfield texture.
  //Second and third GLuint are the two position textures. 
  void advect(bool,GLuint,GLuint,GLuint,GLuint,GLuint,float);

  void updatePrime(bool,GLuint,GLuint,GLuint,GLuint,GLuint,GLuint,GLuint,float);
  // included two more arguments in the above function for position textures. --Balli(04/12/07)

  //Sets the ustar and sigma values which can now be used to create
  //the lambda texture and CoE/2 values. 
  void setUstarAndSigmas(float);
  
  void createPosImages(bool);
  void createPrimeImages(bool);

  void setupMeanVel_shader(int);
  void findMeanVel(bool,GLuint,GLuint,GLuint,GLuint,GLuint,GLuint,GLuint);

  void setupCurrVel_shader(int);
  void updateCurrVel(bool,GLuint,GLuint,GLuint,GLuint,GLuint);

  void setRandomTexCoords();

  void setBuildingParameters(int,float*,float*,float*,float*,float*,float*);
  void addBuildingsInWindField(GLuint,int);
  void setQuicFilesPath(std::string);
  
  bool outputPrime;

  bool osgPlume;
   
  typedef struct{
    float u;
    float v;
    float w;
    float id;
  }wind;

  
  typedef struct{
    int u;
    int v;
    int w;
    int id;
  }intCells;

  //!!!ASK BALLI ABOUT THIS!!!!
  typedef struct{
    float c;
  }cellType;
  
  float getMinDistance(int,int,int); //Added function for obtaining minimum distance to a surface --Balli

  wind* sig;

  GLenum meanVelBuffer0, meanVelBuffer1, currVelBuffer;

  std::string randomFile;
  bool alreadyOpen;

  //RandomTexCoords;
  float t1,t2;

  //Max and Min Tau values 
  float tauMax[4];// tau11Max,tau22Max,tau33Max,tau13Max;  
  float tauMin[4];// tau11Min,tau22Min,tau33Min,tau13Min;

  //Max and Min Tau values per height value(local max and min)
  float* tauLocalMax;
  float* tauLocalMin;

  Matrix* tau;

  //Max and Min Wind Field values
  float windMax[4];
  float windMin[4];


  int nx;
  int ny;
  int nz;

  GLenum texType;

 private:
  float min,max; 

  void test1();
  void randomWindField();
  void uniformUWindField();
  void variedUWindField();
  void QUICWindField();
  
  void initCellType();

  void printPrime(bool,bool);

  void updateMaxandMinTaus(float,float,float,float);
  void find_tauLocalMax();

  void updateMaxandMinWindVel(float,float,float,float);

  wind* wind_vel;
  cellType* cellQuic;

  int numBuild;
  float* xfo;
  float* yfo;
  float* zfo;
  float* ht;
  float* wti;
  float* lti;

  //Number of Particles
  int twidth, theight;

  //Width and height of wind,lambda textures
  int width, height;

  float ustar,sigU,sigV,sigW;

  GLSLObject init_shader, pass1_shader, prime_shader, mrt_shader;
  GLSLObject nonGaussian_shader, reflection_shader, meanVel_shader;
  GLSLObject currVel_shader, multipleBuildings_shader;

  //Variables for prime shader
  GLint uniform_prime, uniform_windTex, uniform_random,uniform_pos;
  GLint uniform_dt,uniform_lambda;
  GLint uniform_randomTexCoordOffset, uniform_randomTexWidth, uniform_randomTexHeight;

  //Variables for advect shader
  GLint uniform_postex, uniform_wind, uniform_randomTexture;
  GLint uniform_primePrev, uniform_primeCurr, uniform_timeStep;
  GLint uniform_tau_dz, uniform_duvw_dz;

  //Uniform variables for building
  GLint uniform_numBuild;
  GLint uniform_xfo, uniform_yfo, uniform_zfo, uniform_ht, uniform_wti, uniform_lti;
  GLint uniform_xfo2, uniform_yfo2, uniform_zfo2, uniform_ht2, uniform_wti2, uniform_lti2;

  GLint uniform_buildings, uniform_cellType;

  //Uniform variables for Mean Velocity shader
  GLint uniform_prevMean, uniform_currVel, uniform_position, uniform_windVel;
  GLint unir;
  
  //Uniform variables for the Current Velocity shader
  GLint uniform_currentPrime, uniform_windVelocity,uniform_partPos, uniform_prevPartPos;


  
  GLfloat* buffer_mem;

  GLint currentbuffer;

  std::string quicFilesPath;

};
#endif //__PC_H__
