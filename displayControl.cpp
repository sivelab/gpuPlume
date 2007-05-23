#include <iostream>
#include <math.h>
#include "displayControl.h"
#include "glm.h"

static char text_buffer[128];

DisplayControl::DisplayControl(int x, int y, int z, GLenum type)
{
  nx = x;
  ny = y;
  nz = z;
  texType = type;

  eye_pos[0] = nx+50;
  eye_pos[1] = 0;
  eye_pos[2] = 0;
  
  eye_gaze[0] = eye_pos[0] - 5;
  eye_gaze[1] = 0;
  eye_gaze[2] = 0;

  rotate_sphere = false;
  rotate_object = false;
  translate_view = false;
  azimuth = 0.0;
  elevation = 0.0;

  frame_rate = true;
  visual_layer = -1;
  osgPlume = false;
  
  //This shader is used to make final changes before rendering to the screen
  render_shader.addShader("Shaders/particleVisualize_vp.glsl", GLSLObject::VERTEX_SHADER);
  render_shader.addShader("Shaders/particleVisualize_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  render_shader.createProgram();
  
  // Create a high resolution clock timer - only works on Linux, x86
  // systems.  The basic timer works on Windows.  Setting the argument
  // to true will have no affect on windows implementations.
  clock_timer = new Timer(true);
  
}

void DisplayControl::drawVisuals(GLuint vertex_buffer,GLuint texid3, int numInRow, 
				 int twidth, int theight)
{
  if(!osgPlume){
    gluLookAt( eye_pos[0], eye_pos[1], eye_pos[2],
	     eye_gaze[0], eye_gaze[1], eye_gaze[2],
	     0, 0, 1 );

    if (!rotate_sphere) 
      {
	// allow rotation of this object
	glRotatef(elevation, 0,1,0);
	glRotatef(azimuth, 0,0,1);
	// glTranslatef(0,0,5.0);
      }
  }
  
  // render the vertices in the VBO (the particle positions) as points in the domain
  glBindBufferARB(GL_ARRAY_BUFFER, vertex_buffer);
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(4, GL_FLOAT, 0, 0);
  glPointSize(3.0);
  render_shader.activate();
  glDrawArrays(GL_POINTS, 0, twidth*theight);
  render_shader.deactivate();

  drawAxes();
  if(!osgPlume){
    if(draw_buildings){
      drawFeatures();
    }
  }
  drawGround();

  drawLayers(texid3, numInRow);  

  // spit out frame rate
  if(!osgPlume){
    if (frame_rate){
      drawFrameRate(twidth, theight);
    }
  }
  
}
void DisplayControl::increaseVisualLayer(){
  visual_layer++;
  if(visual_layer > nz) visual_layer = nz;

}
void DisplayControl::decreaseVisualLayer(){
  visual_layer--;
  if(visual_layer < -1) visual_layer = -1;
}
void DisplayControl::setEyeValues(float change){
  eye_pos[0] = eye_pos[0] + change;
  eye_gaze[0] = eye_pos[0] - 5.0;
}
void DisplayControl::setAzimuth(float change, float rate){
  azimuth = azimuth + change*rate;
}
void DisplayControl::setElevation(float change, float rate){
  elevation = elevation + change*rate;
}
void DisplayControl::drawGround(){

  glDisable(texType);
  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D, displayTex[0]);

  glBegin(GL_QUADS);
  {
    glTexCoord2f(0,0);       glVertex3f(0.0,0.0,0.0);
    glTexCoord2f(1,0);       glVertex3f(nx,0.0,0.0);
    glTexCoord2f(1,1);       glVertex3f(nx,ny,0.0);
    glTexCoord2f(0,1);       glVertex3f(0.0,ny,0.0);
  }
  glEnd();

  glBindTexture(GL_TEXTURE_2D,0);
  glDisable(GL_TEXTURE_2D);
  glEnable(texType);

}
void DisplayControl::drawAxes(){
  // query the current line width so we can set it back at the end of
  // the function
  GLint lwidth;
  glGetIntegerv(GL_LINE_WIDTH, &lwidth);
  glLineWidth(3.0);

  glBegin(GL_LINES);
  glColor3f(1.0, 0.0, 0.0);
  glVertex3f(0.0, 0.0, 0.0);
  glVertex3f(nx, 0.0, 0.0);

  glColor3f(0.0, 1.0, 0.0);
  glVertex3f(0.0, 0.0, 0.0);
  glVertex3f(0.0, ny, 0.0);

  glColor3f(0.0, 0.0, 1.0);
  glVertex3f(0.0, 0.0, 0.0);
  glVertex3f(0.0, 0.0, nz);
  glEnd();
  //This draws a label for the axis by drawing a textured quad.
  /*
  glDisable(texType);
  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D, axisLabel[0]);
  glPushMatrix();
  glRotatef(180.0, 0.0, 0.0, 1.0);

  glBegin(GL_QUADS);
  {
    glTexCoord2f(0, 0);      glVertex3f(-0.5, -0.5, nz+1.0);
    glTexCoord2f(1, 0);      glVertex3f(0.5, -0.5, nz+1.0);
    glTexCoord2f(1, 1);      glVertex3f(0.5, 0.5, nz+1.0);
    glTexCoord2f(0, 1);      glVertex3f(-0.5, 0.5, nz+1.0);
  }
  

  glEnd();

  glPopMatrix();
  glDisable(GL_TEXTURE_2D);
  glEnable(texType);*/

  // set the line width back to what it was
  glLineWidth(lwidth);

}
void DisplayControl::drawLayers(GLuint texId, int numInRow){
  if (visual_layer >= 0 && visual_layer < nz)
    {
      //glEnable(GL_LIGHTING);
      //glEnable(GL_LIGHT0);
      glEnable(GL_COLOR_MATERIAL);
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);     

      glEnable(texType);
      glBindTexture(texType, texId);
      glTexParameteri(texType, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(texType, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      
      //Since the alpha value is the epsilon value, we need
      //to make sure alpha value of displayed layer is 1.0;
      static GLfloat col[4] = {0.0,0.0,0.0,1.0};
      glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_COLOR, col);

      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_COMBINE);  
      glTexEnvf(GL_TEXTURE_ENV, GL_COMBINE_ALPHA, GL_ADD);
      glTexEnvf(GL_TEXTURE_ENV, GL_SRC0_ALPHA, GL_TEXTURE);
      
      //glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);

      // The s and t parameters reference pixels on the UVW texture
      // map.  Because we are using texture rectangles, the s and t
      // parameters do not need to be in normalized device coordinates
      // so the values represent actual pixels positions in the
      // texture map.
      int s = 0;
      int t = 0;

      // The texture we use here is the packed texture containing all
      // cells of the 3D uvw wind field data.  It is packed because we
      // flatten the 3D structure into a series of discrete 2D
      // elements. The variable numInRow is the number of these
      // discrete 2D layers that can fit in each row the texture. 

      // The coordinate frame we use is with Y up, so X and Z (at Y=0)
      // forms the ground plane.
      //

      // s (or the value in the x dimension of the texture) can be
      // determined with a mod of the layer by the number of layers
      // that can be held in each row of the texutre. [ COMPLETE DESCRIPTION ]
      s = (int)(visual_layer % numInRow) * nx;

      // t (or the value in the y dimension of the texture) can be 
      // calculated by the floor of the layer to be visualized divided
      // by the number of layers that can be held in each row of
      // the texture. 
      t = (int)(floor(visual_layer/(float)numInRow) * ny);

      // Create a quad at this layer with 50% transparency
      glColor4f(1.0, 1.0, 1.0, 0.8);
      glBegin(GL_QUADS);
      {
	glNormal3f(0.0, 0.0, 1.0);
	glTexCoord2f(s, t);         glVertex3f(0, 0, visual_layer);
	glTexCoord2f(s+nx, t);      glVertex3f(nx, 0,visual_layer);
	glTexCoord2f(s+nx, t+ny);   glVertex3f(nx, ny,visual_layer);
	glTexCoord2f(s, t+ny);      glVertex3f(0, ny,visual_layer);
      }
      glEnd();
      glTexParameteri(texType, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(texType, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      //glBindTexture(texType, 0);
      glDisable(texType);
      glDisable(GL_BLEND);
      glDisable(GL_COLOR_MATERIAL);
      //glDisable(GL_LIGHT0);
      //glDisable(GL_LIGHTING);
    }
}

void instanceCube()
{
}
void DisplayControl::initVars(int nb,double* x, double* y, double* z,
				    double* h, double* w, double* l)
{
  numBuild = nb;
  xfo = x;
  yfo = y;
  zfo = z;
  ht = h;
  wti = w;
  lti = l;
  
  glDisable(texType);
  glEnable(GL_TEXTURE_2D);
  glGenTextures(3,displayTex);
  
  createImageTex(displayTex[0], "concrete.ppm");
 
  glDisable(GL_TEXTURE_2D);
  glEnable(texType);

}

void DisplayControl::createImageTex(GLuint texture, char* filename){
  GLubyte* testImage;
  int w, h;

  testImage = glmReadPPM(filename, &w, &h);
  if(testImage == 0){
    std::cout << "Didn't Load Texture File" << std::endl;
  }
  glBindTexture(GL_TEXTURE_2D, texture); 
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, testImage);
}

void DisplayControl::drawFeatures(void)
{
  float grid_scale = 1.0;  // currently, just 1 but likely needs to come from QUICPLUME

  //#ifdef USE_PLUME_DATA
  // Draw the building
  for (int qi=0; qi<numBuild; qi++) 
    {
      glPushMatrix();
      glColor3f(0.5, 0.5, 0.5);

      glTranslatef(xfo[qi]*grid_scale+ (lti[qi]*grid_scale)/2.0,
		   yfo[qi]*grid_scale,
		   zfo[qi]*grid_scale + (ht[qi]*grid_scale)/2.0);

      glScalef(lti[qi]*grid_scale,
	       wti[qi]*grid_scale,
	       ht[qi]*grid_scale);

      glutSolidCube(1.0);
      glPopMatrix();
    }
  //#endif
}

//text: draws a string of text with an 18 point helvetica bitmap font
// at position (x,y) in window space(bottom left corner is (0,0).

void DisplayControl::drawFrameRate(int twidth, int theight)
{
  // record end clock time
  graphics_time[1] = clock_timer->tic();

  float avg_frame_rate = 1.0/( clock_timer->deltas( graphics_time[0], graphics_time[1] ) );
  sprintf(text_buffer, "%d particles, %0d fps", twidth*theight, (int)(avg_frame_rate+0.5));
  OpenGLText(5, 5, text_buffer);

  // record start clock time
  graphics_time[0] = clock_timer->tic();
}

void DisplayControl::OpenGLText(int x, int y, char* s)
{
  int lines;
  char* p;
  GLint* vp = new GLint[4];
  glGetIntegerv(GL_VIEWPORT,vp);

  // glDisable(GL_LIGHTING);
  glDisable(texType);
  glDisable(GL_DEPTH_TEST);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0, vp[2], 
	  0, vp[3], -1, 1);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glColor3ub(0, 0, 0);
  glRasterPos2i(x+1, y-1);
  for (p=s, lines=0; *p; p++) {
    if (*p == '\n') {
      lines++;
      glRasterPos2i(x+1, y-1-(lines*18));
    }
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *p);
  }
  glColor3ub(255, 255, 0);
  glRasterPos2i(x, y);
  for (p=s, lines=0; *p; p++) {
    if (*p == '\n') {
      lines++;
      glRasterPos2i(x, y-(lines*18));
    }
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *p);
  }
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glEnable(GL_DEPTH_TEST);
  glEnable(texType);
  // glEnable(GL_LIGHTING);

}
