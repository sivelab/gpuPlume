#include "displayControl.h"
#include <iostream>
#include <math.h>
#include <GL/glut.h>

using namespace std;

static char text_buffer[128];

DisplayControl::DisplayControl(int x, int y, int z, GLenum type)
{
  nx = x;
  ny = y;
  nz = z;
  texType = type;

  
  // Create a high resolution clock timer - only works on Linux, x86
  // systems.  The basic timer works on Windows.  Setting the argument
  // to true will have no affect on windows implementations.
  clock_timer = new Timer(true);
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

  // set the line width back to what it was
  glLineWidth(lwidth);
}
void DisplayControl::drawLayers(int layer, GLuint texId, int numInRow){
  if (layer >= 0 && layer < ny)
    {
      glEnable(GL_LIGHTING);
      glEnable(GL_LIGHT0);
      glEnable(GL_COLOR_MATERIAL);
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

      glEnable(texType);
      glBindTexture(texType, texId);
      glTexParameteri(texType, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(texType, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

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
      // elements.

      // The coordinate frame we use is with Y up, so X and Z (at Y=0)
      // forms the ground plane.
      //

      // s (or the value in the x dimension of the texture) can be
      // determined with a mod of the layer by the number of elements
      // in a row. [ COMPLETE DESCRIPTION ]
      s = (layer % numInRow) * nz;

      // t (or ...) can be calculated by 
      t = (int)(floor(layer/(float)numInRow) * nx);

      // Create a quad at this layer with 50% transparency
      glColor4f(1.0, 1.0, 1.0, 0.8);
      glBegin(GL_QUADS);
      {
	glNormal3f(0.0, 1.0, 0.0);
	glTexCoord2f(s, t);           glVertex3f(0, layer, 0);
	glTexCoord2f(s+nz-1, t);      glVertex3f(nz, layer, 0);
	glTexCoord2f(s+nz-1, t+nx-1); glVertex3f(nz, layer, nx);
	glTexCoord2f(s, t+nx-1);      glVertex3f(0, layer, nx);
      }
      glEnd();
      glTexParameteri(texType, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(texType, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glDisable(texType);
      glDisable(GL_BLEND);
      glDisable(GL_COLOR_MATERIAL);
      glDisable(GL_LIGHT0);
      glDisable(GL_LIGHTING);
    }
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

  // glDisable(GL_LIGHTING);
  glDisable(texType);
  //glDisable(GL_TEXTURE_2D);
  glDisable(GL_DEPTH_TEST);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0, glutGet(GLUT_WINDOW_WIDTH), 
	  0, glutGet(GLUT_WINDOW_HEIGHT), -1, 1);
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
  //glEnable(GL_TEXTURE_2D);
  glEnable(texType);
  // glEnable(GL_LIGHTING);

}
