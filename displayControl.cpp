#include <iostream>
#include <math.h>
#include <assert.h>
#include "displayControl.h"
#include "glErrorUtil.h"

static char text_buffer[128];
#ifdef WIN32
#define M_PI 3.14159
#endif

DisplayControl::DisplayControl(int x, int y, int z, GLenum type)
{
  nx = x;
  ny = y;
  nz = z;
  texType = type;

  eye_pos[0] = nx;
  eye_pos[1] = 0;
  eye_pos[2] = 5;
  
  eye_gaze[0] = -1.0;
  eye_gaze[1] = 0;
  eye_gaze[2] = 0;

  angle = M_PI;
  yangle = 0.0;

  xslide = 0.0;
  yslide = 1.0;

  change_height = false;
  change_look = false;
  rotate_around = false;
  azimuth = 0.0;
  elevation = 0.0;

  frame_rate = true;
  visual_layer = -1;
  osgPlume = false;

  // Set particle visual state to point based particles initially
  particle_visual_state = PARTICLE_SPRITE;  // POINT or SPRITE;

  //This shader is used to make final changes before rendering to the screen
  render_shader.addShader("Shaders/particleVisualize_vp.glsl", GLSLObject::VERTEX_SHADER);
  render_shader.addShader("Shaders/particleVisualize_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  render_shader.createProgram();

  // for point sprites, we need the uniform variables for the texture
  // units that hold the point sprite and the normal map
  uniform_pointsprite_tex = render_shader.createUniform("pointsprite_texunit");
  uniform_normalmap_tex = render_shader.createUniform("pointspritenormal_texunit");

  // determine which visual to use
  uniform_pointsprite_visuals = render_shader.createUniform("point_visuals");

  // POINT_SPRITE
  //
  // Initialize the texture data used for sprite based particle
  // visuals.  This currently means creating two textures: one
  // representing the base shape/color of the particle (white circle),
  // and a second to represent a normal map for a sphere to provide
  // Phong based lighting to the particle representations in a shader.
  createPointSpriteTextures();

  // float quadratic[] =  { 0.0f, 0.0f, 0.1f };
  // glPointParameterfvARB( GL_POINT_DISTANCE_ATTENUATION_ARB, quadratic );
  float maxSize = 0.0f, currSize = 1.0f;
  glGetFloatv(GL_POINT_SIZE_MAX_ARB, &maxSize);
  glGetFloatv(GL_POINT_SIZE, &currSize);
  glPointParameterfARB(GL_POINT_SIZE_MAX_ARB, currSize + 3.0);
  glPointParameterfARB(GL_POINT_SIZE_MIN_ARB, 1.0);
  
  // Create a high resolution clock timer - only works on Linux, x86
  // systems.  The basic timer works on Windows.  Setting the argument
  // to true will have no affect on windows implementations.
  clock_timer = new Timer(true);

  graphics_time[0] = clock_timer->tic();
  
}

void DisplayControl::drawVisuals(GLuint vertex_buffer,GLuint texid3, GLuint color_buffer, 
				 int numInRow, int twidth, int theight)
{
  drawSky();

  if(!osgPlume){
    gluLookAt( eye_pos[0], eye_pos[1], eye_pos[2],
	     eye_gaze[0]+eye_pos[0], eye_gaze[1]+eye_pos[1], 
	       eye_gaze[2]+eye_pos[2], 0, 0, 1 );

    // allow rotation of this object
    //glRotatef(elevation, 0,1,0);
    //glRotatef(azimuth, 0,0,1);  
  }
  
  drawAxes();
  
  if(!osgPlume)
    drawGrid();

  if(draw_buildings){
    drawFeatures();
  }
  
  drawGround();

  // drawLayers(texid3, numInRow);  

  // render the vertices in the VBO (the particle positions) as points in the domain
  
  if(color_buffer != 0){
    glEnableClientState(GL_COLOR_ARRAY); 
    glBindBufferARB(GL_ARRAY_BUFFER, color_buffer);
    glColorPointer(4, GL_FLOAT, 0, 0);
  }
  glBindBufferARB(GL_ARRAY_BUFFER, vertex_buffer);
  glVertexPointer(4, GL_FLOAT, 0, 0);
  
  glEnableClientState(GL_VERTEX_ARRAY); 
  
  render_shader.activate();  
  if(color_buffer == 0)
    glColor4f(1.0,1.0,1.0,1.0);

  if (particle_visual_state == PARTICLE_SPRITE)
    {
      glPointSize(6.0);
      glUniform1iARB(uniform_pointsprite_visuals, 1);

      glEnable(GL_TEXTURE_2D);

      glActiveTextureARB(GL_TEXTURE0_ARB);
      glUniform1iARB(uniform_pointsprite_tex, 0);
      glBindTexture(GL_TEXTURE_2D, displayTex[0]); // point_sprite_textures[0]);

      glActiveTextureARB(GL_TEXTURE1_ARB);
      glUniform1iARB(uniform_normalmap_tex, 1);
      glBindTexture(GL_TEXTURE_2D, point_sprite_textures[1]);

      glActiveTextureARB(GL_TEXTURE0_ARB);

      glTexEnvf(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);
      glEnable(GL_POINT_SPRITE_ARB);
      glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
    }
  else 
    {
      glPointSize(3.0);
      glUniform1iARB(uniform_pointsprite_visuals, 0);
    }

  glDrawArrays(GL_POINTS, 0, twidth*theight);

  if (particle_visual_state == PARTICLE_SPRITE)
    {
      glDisable(GL_POINT_SPRITE_ARB);

      // glBindTexture(GL_TEXTURE_2D, 0);
      glDisable(GL_TEXTURE_2D);
    }

  render_shader.deactivate();


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
void DisplayControl::moveForwardorBack(float change){
  //eye_pos[0] = eye_pos[0] + change;
  //eye_gaze[0] = eye_pos[0] - 5.0;
  if(change < 0.0){
    eye_pos[0] -= eye_gaze[0];
    eye_pos[1] -= eye_gaze[1];
    eye_pos[2] -= eye_gaze[2];
  }
  else{
    eye_pos[0] += eye_gaze[0];
    eye_pos[1] += eye_gaze[1];
    eye_pos[2] += eye_gaze[2];
  }
}
void DisplayControl::slideLeftorRight(float direction){

  if(direction < 0.0){
    eye_pos[0] -= xslide;
    eye_pos[1] -= yslide;
  }
  else{
    eye_pos[0] += xslide;
    eye_pos[1] += yslide;
  }

}
void DisplayControl::lookUporDown(float change){
  if(change < 0.0)
    yangle = yangle + (M_PI/90.0);
  else
    yangle = yangle - (M_PI/90.0);

  if(yangle > M_PI/2.0)
    yangle = yangle - (M_PI/90.0);
  else if(yangle < -M_PI/2.0)
    yangle = yangle + (M_PI/90.0);

  eye_gaze[2] = sin(yangle);

}
void DisplayControl::setAzimuth(float change, float rate){
  azimuth = azimuth + change*rate;
}
void DisplayControl::setElevation(float change, float rate){
  //elevation = elevation + change*rate;
  eye_pos[2] = eye_pos[2] + change*rate;
  //eye_gaze[2] = eye_gaze[2] + change*rate;
}
void DisplayControl::setRotateAround(float change){

  if(change < 0)
    angle = angle + (M_PI/90.0);
  else 
    angle = angle - (M_PI/90.0);

  if(angle > 2*M_PI)
    angle = 0.0;
  else if(angle < 0.0)
    angle = (2*M_PI - M_PI/90.0);

  eye_gaze[0] = cos(angle);
  eye_gaze[1] = sin(angle);

  xslide = cos(angle-M_PI/2.0);
  yslide = sin(angle-M_PI/2.0);


}
void DisplayControl::drawSky(){
  // double skyr = 0.4, skyg = 0.5, skyb = 0.8;
  double skyr = 0.5, skyg = 0.5, skyb = 0.6;
   
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glDisable ( GL_TEXTURE_2D );
  glDepthMask(GL_FALSE);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
 
  //Put sky in background
  glBegin(GL_POLYGON);
  glColor4f(0.9,0.9,1.0, 1.0); //mid color
  glVertex3f(-1,-0.5,-1);
  glColor4f(0.9,0.9,1.0, 1.0); //white
  glVertex3f(-1,-1,-1);
  glVertex3f(1,-1,-1);
  glColor4f(0.9,0.9,1.0, 1.0); //mid color
  glVertex3f(1,-.5,-1);
  glColor4f(skyr, skyg, skyb, 1.0);
  glVertex3f(1,1,-1);
  glVertex3f(-1,1,-1);
  glEnd();
 
  glDepthMask(GL_TRUE);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glPopAttrib();

}
void DisplayControl::drawGround(){

  glDisable(texType);
  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D, displayTex[0]);
  
  glColor4f(1.0,1.0,1.0,1.0);
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
void DisplayControl::drawGrid(){
  GLint lwidth;
  glGetIntegerv(GL_LINE_WIDTH, &lwidth);

  glLineWidth(0.5);
  glBegin(GL_LINES);
  {
    glColor3f(0.4,0.4,0.4);
    for(int i=-ny; i <= ny*2; i+=5){
      glVertex3f(-nx,i,-0.1);
      glVertex3f(nx*2,i,-0.1);
    }
    for(int i=-nx; i <= nx*2; i+=5){
      glVertex3f(i,-ny,-0.1);
      glVertex3f(i,ny*2,-0.1);
    }
  
  }
  glEnd();

  glLineWidth(lwidth);
}

void DisplayControl::drawAxes(){
  // query the current line width so we can set it back at the end of
  // the function
  //int numLines = ny + 10;

  GLint lwidth;
  glGetIntegerv(GL_LINE_WIDTH, &lwidth);
  glLineWidth(3.0);

  // Modified colors of axes to represent the manner in which
  // direction is being visualized.
  // 
  // Which makes me think we might want to make the axes' colors
  // dependent on what is being visualized.

  glBegin(GL_LINES);
  glColor3f(0.0, 0.0, 1.0);
  glVertex3f(0.0, 0.0, 0.0);
  glColor3f(1.0, 1.0, 0.0);
  glVertex3f(nx, 0.0, 0.0);

  glColor3f(0.0, 1.0, 0.0);
  glVertex3f(0.0, 0.0, 0.0);
  glColor3f(1.0, 0.0, 0.0);
  glVertex3f(0.0, ny, 0.0);

  glColor3f(0.0, 0.0, 0.0);
  glVertex3f(0.0, 0.0, 0.0);
  glColor3f(1.0, 1.0, 1.0);
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
      glPushMatrix();
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

      glPopMatrix();
    }
  
}

void instanceCube()
{
}
void DisplayControl::initVars(int nb,float* x, float* y, float* z,
				    float* h, float* w, float* l)
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
  createImageTex(displayTex[1], "building.ppm");
  createImageTex(displayTex[2], "buildingRoof.ppm");
 
  glDisable(GL_TEXTURE_2D);
  glEnable(texType);

}

void DisplayControl::createImageTex(GLuint texture, char* filename){
  GLubyte* testImage;
  int w, h;

  testImage = readPPM(filename, &w, &h);
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
  //float grid_scale = 1.0;  // currently, just 1 but likely needs to come from QUICPLUME
  glDisable(texType);
  glEnable(GL_TEXTURE_2D);
  
  // Draw the building
  for (int qi=0; qi<numBuild; qi++) 
    {
      //glPushMatrix();
      glColor4f(1.0, 1.0, 1.0,1.0);

      /*glTranslatef(xfo[qi]*grid_scale+ (lti[qi]*grid_scale)/2.0,
		   yfo[qi]*grid_scale,
		   zfo[qi]*grid_scale + (ht[qi]*grid_scale)/2.0);

      glScalef(lti[qi]*grid_scale,
	       wti[qi]*grid_scale,
	       ht[qi]*grid_scale);


	       glutSolidCube(1.0);*/
      //glPopMatrix();
      glBindTexture(GL_TEXTURE_2D, displayTex[1]);

      glBegin(GL_QUADS);
      {
	glTexCoord2f(0,0);       glVertex3f(xfo[qi],yfo[qi]-(wti[qi]/2.0),zfo[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]+lti[qi],yfo[qi]-(wti[qi]/2.0),zfo[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]+lti[qi],yfo[qi]-(wti[qi]/2.0),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi],yfo[qi]-(wti[qi]/2.0),zfo[qi]+ht[qi]);
	
	glTexCoord2f(0,0);       glVertex3f(xfo[qi],yfo[qi]+(wti[qi]/2.0),zfo[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]+lti[qi],yfo[qi]+(wti[qi]/2.0),zfo[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]+lti[qi],yfo[qi]+(wti[qi]/2.0),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi],yfo[qi]+(wti[qi]/2.0),zfo[qi]+ht[qi]);
	
	glTexCoord2f(0,0);       glVertex3f(xfo[qi],yfo[qi]+(wti[qi]/2.0),zfo[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi],yfo[qi]-(wti[qi]/2.0),zfo[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi],yfo[qi]-(wti[qi]/2.0),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi],yfo[qi]+(wti[qi]/2.0),zfo[qi]+ht[qi]);

	glTexCoord2f(0,0);       glVertex3f(xfo[qi]+lti[qi],yfo[qi]+(wti[qi]/2.0),zfo[qi]);
	glTexCoord2f(1,0);       glVertex3f(xfo[qi]+lti[qi],yfo[qi]-(wti[qi]/2.0),zfo[qi]);
	glTexCoord2f(1,1);       glVertex3f(xfo[qi]+lti[qi],yfo[qi]-(wti[qi]/2.0),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);       glVertex3f(xfo[qi]+lti[qi],yfo[qi]+(wti[qi]/2.0),zfo[qi]+ht[qi]);
      }
      glEnd();

      glBindTexture(GL_TEXTURE_2D,displayTex[2]);
      glBegin(GL_QUADS);
      {
	glTexCoord2f(0,0);	glVertex3f(xfo[qi],yfo[qi]-(wti[qi]/2.0),zfo[qi]+ht[qi]);
	glTexCoord2f(1,0);	glVertex3f(xfo[qi]+lti[qi],yfo[qi]-(wti[qi]/2.0),zfo[qi]+ht[qi]);
	glTexCoord2f(1,1);	glVertex3f(xfo[qi]+lti[qi],yfo[qi]+(wti[qi]/2.0),zfo[qi]+ht[qi]);
	glTexCoord2f(0,1);	glVertex3f(xfo[qi],yfo[qi]+(wti[qi]/2.0),zfo[qi]+ht[qi]);
		
      }
      glEnd();
    }
	
  glDisable(GL_TEXTURE_2D);
	
  glEnable(texType);
	
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

GLubyte* DisplayControl::readPPM(char* filename, int* width, int* height)
{
    FILE* fp;
    int i, w, h, d;
    unsigned char* image;
    char head[70];          /* max line <= 70 in PPM (per spec). */
    
    fp = fopen(filename, "rb");
    if (!fp) {
        perror(filename);
        return NULL;
    }
    
    /* grab first two chars of the file and make sure that it has the
       correct magic cookie for a raw PPM file. */
    fgets(head, 70, fp);
    if (strncmp(head, "P6", 2)) {
        fprintf(stderr, "%s: Not a raw PPM file\n", filename);
        return NULL;
    }
    
    /* grab the three elements in the header (width, height, maxval). */
    i = 0;
    while(i < 3) {
        fgets(head, 70, fp);
        if (head[0] == '#')     /* skip comments. */
            continue;
        if (i == 0)
            i += sscanf(head, "%d %d %d", &w, &h, &d);
        else if (i == 1)
            i += sscanf(head, "%d %d", &h, &d);
        else if (i == 2)
            i += sscanf(head, "%d", &d);
    }
    
    /* grab all the image data in one fell swoop. */
    image = (unsigned char*)malloc(sizeof(unsigned char)*w*h*3);
    fread(image, sizeof(unsigned char), w*h*3, fp);
    fclose(fp);
    
    *width = w;
    *height = h;
    return image;
}

void DisplayControl::createPointSpriteTextures()
{
  unsigned int width = 128;
  unsigned int height = 128;

  // for this, I want width & height to be the same
  assert(width == height);
  GLfloat *data = new GLfloat[ width * height * 4 ];

  glEnable(GL_TEXTURE_2D);
  glGenTextures(2, point_sprite_textures);

  glBindTexture(GL_TEXTURE_2D, point_sprite_textures[0]);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  // make the shape/intensity map
  for (unsigned int y=0; y<height; y++)
    for (unsigned int x=0; x<width; x++)
      {
	int idx = y*width*4 + x*4;

	// convert to sphere coordinates with radius half of width/height
	double r = (double)width/2.0;
	double xr = x - r;
	double yr = y - r;

	// if coordinate is on or inside the circle, keep going
	double mag = sqrt(xr*xr + yr*yr);
	if (mag <= r)
	  {
	    // White for now... may use later
	    data[idx] = 1.0; data[idx+1] = 1.0; data[idx+2] = 1.0; data[idx+3] = 1.0;
	  }
	else 
	  {
	    // store transparent white in image space
	    data[idx] = data[idx+1] = data[idx+2] = 1.0;
	    data[idx+3] = 0.0;
	  }
      }

  gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, width, height, GL_RGBA, GL_FLOAT, data);
  // glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_FLOAT, data);
  
  // 
  // make the normal map
  //
  glBindTexture(GL_TEXTURE_2D, point_sprite_textures[1]);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  for (unsigned int y=0; y<height; y++)
    for (unsigned int x=0; x<width; x++)
      {
	int idx = y*width*4 + x*4;

	// convert to sphere coordinates with radius half of width/height
	double r = (double)width/2.0;
	double xr = x - r;
	double yr = y - r;

	// if coordinate is on or inside the circle, keep going
	double mag = sqrt(xr*xr + yr*yr);
	if (mag <= r)
	  {
	    // calculate what z would be for this x, y, radius combination
	    double z = sqrt( r*r - xr*xr - yr*yr );

	    // convert to positive half space since some x and y will
	    // be negative (bad for RGBA values)
	    xr /= r;
	    xr = (xr + 1.0) / 2.0;

	    yr /= r;
	    yr = (yr + 1.0) / 2.0;

	    z /= r;
	    z = (z + 1.0) / 2.0;

	    // store in image space
	    data[idx] = xr;
	    data[idx+1] = yr;
	    data[idx+2] = z;
	    data[idx+3] = 1.0;
	  }
	else 
	  {
	    // store transparent white in image space
	    data[idx] = data[idx+1] = data[idx+2] = 1.0;
	    data[idx+3] = 0.0;
	  }
      }
	
  // glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_FLOAT, data);
  gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, width, height, GL_RGBA, GL_FLOAT, data);

  glDisable(GL_TEXTURE_2D);

  delete [] data;
}
