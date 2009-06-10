
#include <iostream>
#include <GL/glut.h>
#include <math.h>

// //////////////////////////////////////
// BEGIN -----> QUIC PLUME FORTRAN REFERENCES
// //////////////////////////////////////

#define USE_PLUME_DATA

#ifdef USE_PLUME_DATA

extern "C"
{
  void readfiles_();
}

// Domain size stored in nx, ny, and nz
extern "C" int __datamodule__nx;
extern "C" int __datamodule__ny;
extern "C" int __datamodule__nz;

// UVW contains the wind field
extern "C" double* __datamodule__u;
extern "C" double* __datamodule__v;
extern "C" double* __datamodule__w;

#endif
// //////////////////////////////////////
// END ----> QUIC PLUME FORTRAN REFERENCES
// //////////////////////////////////////

int width, height;
static GLuint texName[2];
static GLuint list;

void init(void);
void display(void);

typedef struct{
  float u;
  float v;
  float w;

}wind;

void init(void)
{
  glClearColor(1.0, 0.0, 0.0, 1.0);

  glShadeModel(GL_SMOOTH);
  glEnable(GL_DEPTH_TEST);

  int nx = 60;
  int ny = 20;
  int nz = 60;

  /*wind data3d[nx][nz][ny];

  for(int k = 0; k < ny; k++){   
    for(int i = 0; i < nx; i++){
      for(int j = 0; j < nz; j++){
	//data3d[i][j][k] = 100*k + 10*i + j;
	if(i == 0){
	  	data3d[i][j][k].u = 1.0;
		data3d[i][j][k].v = 1.0;
		data3d[i][j][k].w = 0;
	}
	else{
	data3d[i][j][k].u = 0;
	data3d[i][j][k].v = 0;
	data3d[i][j][k].w = 1.0;
	}
      }
    }
  }
  */
  
  int num;
  int total = nx*ny*nz;
  width = (int)sqrt(total);
  width = width - (width%nz);
  bool done = false;
  while(!done){  
    num = width/nz;
    if((num*num) >= ny){
      done = true;
    }
    else{
      width = width+nz;
    }
  }
  if(width % 2 != 0) width++;
  height = width;
  /*
std:: cout << height << " " << width << std::endl;
  //wind data2d[height][width];
  wind **data2d = new wind*[height];
  for(int i = 0;  i < height; i++) data2d[i] = new wind[width];
  
  for(int i = 0; i < height; i++){
    for(int j = 0; j < width; j++){
      data2d[i][j].u = 0.0;
      data2d[i][j].v = 0.0;
      data2d[i][j].w = 0.0;
    }
  }
  
  //Convert this to 2D Structure
  int numInRow = (width-(width%nz))/nz;
  int s = 0;
  int t = 0;
  
  for(int k = 0; k < ny; k++){
    if((k%numInRow == 0) && (k != 0)){
      s += nx;
      t = 0;
    }
    
    for(int i = s; i < (s+nx); i++){
      for(int j = t; j < (t+nz); j++){
	data2d[i][j].u = data3d[i%nx][j%nz][k].u;
	data2d[i][j].v = data3d[i%nx][j%nz][k].v;
	data2d[i][j].w = data3d[i%nx][j%nz][k].w;
      }
    }
    t += nz;
  }

  GLfloat *data = new GLfloat[ width * height * 4 ];
  for(int i = 0; i < height; i++){
    for(int j = 0; j < width; j++){
      int idx = i*width*4 + j*4;
      data[idx] = data2d[i][j].u;
      data[idx+1] = data2d[i][j].v;
      data[idx+2] = data2d[i][j].w;
      data[idx+3] = 1.0;
    }
    }*/
  /*
  int qi, qj, qk;
    int p2idx = 0, texidx = 0;
    int ncols = width / nz;
    int nrows = height / nx;
    int row = 0;
    int sz = 4;
    GLfloat *data = new GLfloat[ width * height * sz ];
    
    std::cout << "Creating tiled, flat 3D texture: ncols = " << ncols << ", nrows = " << nrows << ", nz=" << nz << std::endl;

    for (qk=0; qk<ny; qk++) 
      for (qi=0; qi<nx; qi++)
	for (qj=0; qj<nz; qj++)
	  {
	    // Correct way to index
	    p2idx = qk*nz*nx + qi*nz + qj;
	    
	    row = qk / ncols;
	    texidx = row * width * nx * sz +
	      qi * width * sz +
	      qk % ncols * nx * sz +
	      qj * sz;

	    //assert(texidx < (width*height*sz));
	    
	    if ( qi == 0 || qj == 0 || qi == nz-1 || qj == nx-1)
	    
	      {
		// set boundary conditions...
		data[texidx] = 1.0;
		data[texidx+1] = 1.0;
		data[texidx+2] = 1.0;
		data[texidx+3] = 1.0;
	      }		
	    else 
	      {
		data[texidx] = qk*.1;
		data[texidx+1] = 0;
		data[texidx+2] = qk*.2;
		data[texidx+3] = 1.0;
		/*
		data[texidx] = __datamodule__e[p2idx];  // store e in the first component
		data[texidx+1] = __datamodule__f[p2idx];  // store f in the first component
		data[texidx+2] = __datamodule__g[p2idx];  // store g in the first component
		data[texidx+3] = __datamodule__h[p2idx];  // store h in the first component
		
	      }
	  }*/
      double val1, val2, val3;
      int qi, qj, qk, count = 0;
      int p2idx = 0, texidx = 0;
      int ncols = width / nz;
      int nrows = height / nx;
      int row = 0;
      int sz = 4;
      GLbyte *data = new GLbyte[ width * height * sz ];

      std::cout << "Creating tiled, flat 3D texture: ncols = " << ncols << ", nrows = " << nrows << ", nz=" << nz << std::endl;
      double maxmag = 0, mag;
      double maxu = 0, maxv = 0, maxw = 0;
      for (qk=0; qk<ny; qk++) 
	for (qj=0; qj<nx; qj++)
	  for (qi=0; qi<nz; qi++)
	    {
	      // Correct way to index
	      p2idx = qk*nx*nz + qj*nz + qi;
	      
	      if (__datamodule__u[p2idx] > maxu)
		maxu = __datamodule__u[p2idx];

	      if (__datamodule__v[p2idx] > maxv)
		maxv = __datamodule__v[p2idx];

	      if (__datamodule__w[p2idx] > maxw)
		maxw = __datamodule__w[p2idx];

	      mag = sqrt(__datamodule__u[p2idx]*__datamodule__u[p2idx] +
			 __datamodule__v[p2idx]*__datamodule__v[p2idx] +
			 __datamodule__w[p2idx]*__datamodule__w[p2idx]);
	      if (mag > maxmag) maxmag = mag;
	    }
      std::cout << "max: " << maxu << ", " << maxv << ", " << maxw << ", mag=" << mag << std::endl;

      for (qk=0; qk<ny; qk++) 
	for (qj=0; qj<nx; qj++)
	  for (qi=0; qi<nz; qi++)
	    {
	      // Correct way to index
	      p2idx = qk*nx*nz + qj*nz + qi;
	      
	      row = qk / ncols;
	      texidx = row * width * nx * sz +
		qj * width * sz +
		qk % ncols * nz * sz +
		qi * sz;

	      //assert(texidx < (width*height*sz));

	      if (qk == 0 || qi == 0 || qj == 0 ||
		  qk == ny-1 || qi == nz-1 || qj == nx-1)
		{
		  // set boundary conditions...
		  data[texidx] = 255;
		  data[texidx+1] = 255;
		  data[texidx+2] = 255;
		  data[texidx+3] = 255;
		}		
	      else 
		{
		  double len;
#if 0
		  val1 = __datamodule__pwtx[p2idx] * 255;
		  data[texidx] = (int)val1;
		  data[texidx+1] = (int)val1;
		  data[texidx+2] = (int)val1;
		  data[texidx+3] = 255;
#endif
		  
		  val1 = __datamodule__u[p2idx];
		  val2 = __datamodule__v[p2idx];
		  val3 = __datamodule__w[p2idx];
		  
		  double eps = 1.0e-2;
		  if ((val1 < eps && val1 > -eps) && 
		      (val2 < eps && val2 > -eps) && 
		      (val3 < eps && val3 > -eps))
		    {
		      data[texidx] = 128;
		      data[texidx+1] = 128;
		      data[texidx+2] = 128;
		      data[texidx+3] = 255;
		    }
		  else 
		    {
		      len = sqrt(val1*val1 + val2*val2 + val3*val3);
		      val1 /= maxu;
		      val2 /= maxv;
		      val3 /= maxw;
		  
		  if (len/maxmag > 0.92)
		    {
		      data[texidx] = 255;
		      data[texidx+1] = 0;
		      data[texidx+2] = 0;
		      data[texidx+3] = 255;
		    }
		  else if (len/maxmag > 0.85)
		    {
		      data[texidx] = 255;
		      data[texidx+1] = 128;
		      data[texidx+2] = 0;
		      data[texidx+3] = 255;
		    }
		  else if (len/maxmag > 0.77)
		    {
		      data[texidx] = 255;
		      data[texidx+1] = 191;
		      data[texidx+2] = 0;
		      data[texidx+3] = 255;
		    }
		  else if (len/maxmag > 0.69)
		    {
		      data[texidx] = 255;
		      data[texidx+1] = 255;
		      data[texidx+2] = 0;
		      data[texidx+3] = 255;
		    }
		  else if (len/maxmag > 0.62)
		    {
		      data[texidx] = 191;
		      data[texidx+1] = 255;
		      data[texidx+2] = 0;
		      data[texidx+3] = 255;
		    }
		  else if (len/maxmag > 0.54)
		    {
		      data[texidx] = 128;
		      data[texidx+1] = 255;
		      data[texidx+2] = 0;
		      data[texidx+3] = 255;
		    }
		  else if (len/maxmag > 0.46)
		    {
		      data[texidx] = 0;
		      data[texidx+1] = 255;
		      data[texidx+2] = 0;
		      data[texidx+3] = 255;
		    }
		  else if (len/maxmag > 0.38)
		    {
		      data[texidx] = 0;
		      data[texidx+1] = 255;
		      data[texidx+2] = 128;
		      data[texidx+3] = 255;
		    }
		  else if (len/maxmag > 0.31)
		    {
		      data[texidx] = 0;
		      data[texidx+1] = 255;
		      data[texidx+2] = 191;
		      data[texidx+3] = 255;
		    }
		  else if (len/maxmag > 0.23)
		    {
		      data[texidx] = 0;
		      data[texidx+1] = 255;
		      data[texidx+2] = 255;
		      data[texidx+3] = 255;
		    }
		  else if (len/maxmag > 0.15)
		    {
		      data[texidx] = 0;
		      data[texidx+1] = 191;
		      data[texidx+2] = 255;
		      data[texidx+3] = 255;
		    }
		  else if (len/maxmag > 0.08)
		    {
		      data[texidx] = 0;
		      data[texidx+1] = 128;
		      data[texidx+2] = 255;
		      data[texidx+3] = 255;
		    }
		  else 
		    {
		      data[texidx] = 0;
		      data[texidx+1] = 0;
		      data[texidx+2] = 255;
		      data[texidx+3] = 255;
		    }
		    }
		}
	    }




  glEnable(GL_TEXTURE_RECTANGLE_ARB);
  glBindTexture(GL_TEXTURE_RECTANGLE_ARB, texName[1]);
  glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);

  list = glGenLists(1);
  glNewList(list, GL_COMPILE);
  glBegin(GL_QUADS);
  
  glTexCoord2f(0.0, 0.0); glVertex3f(-3.0, -3.0, 0.0);
  glTexCoord2f(width, 0.0); glVertex3f(3.0, -3.0, 0.0);
  glTexCoord2f(width, height); glVertex3f(3.0, 3.0, 0.0);
  glTexCoord2f(0.0, height); glVertex3f(-3.0, 3.0, 0.0);
 
  glEnd();
  glEndList();

}
void display(void)
{
  
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
 
  glEnable(GL_TEXTURE_RECTANGLE_ARB);
  
  glBindTexture(GL_TEXTURE_RECTANGLE_ARB, texName[1]);
  
  glCallList(list);
 
  glutSwapBuffers();
  //glDisable(GL_TEXTURE_RECTANGLE_ARB);

}

void reshape(int w, int h)
{
  glViewport(0,0,(GLsizei) w, (GLsizei) h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(65.0,(GLfloat) w / (GLfloat) h, 1.0, 50.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0.0, 0.0,5.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0); 
  
}


int main(int argc, char** argv)
{
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  glutInitWindowSize(600, 600);
  glutInitWindowPosition(100, 100);
  glutCreateWindow("test");
  
#ifdef USE_PLUME_DATA
  // Call the PLUME code to read in the data files.
  std::cout << "Reading data using PLUME code..." << std::endl;
  readfiles_();
#endif

  init();
  
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  //glutMouseFunc(mouse);
  //glutKeyboardFunc(keyboard);
  glutMainLoop();
  return 0;

}

