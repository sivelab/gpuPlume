#include "VisualPlane.h"
#include "glErrorUtil.h"

static char number[128];
VisualPlane::VisualPlane(Matrix* tauData,int x, int y, int z, float* TMax,
			 float * TMin, float* max, float* min){

  texType = GL_TEXTURE_3D;
  format = GL_RGBA32F_ARB;

  tau = tauData;

  nx = x;
  ny = y;
  nz = z;

  tauMax = TMax[0];
  if(TMax[1]>tauMax)
    tauMax = TMax[1];
  if(TMax[2]>tauMax)
    tauMax = TMax[2];
  if(TMax[3]>tauMax)
    tauMax = TMax[3];
  tauMin = TMin[0];
  if(TMin[1]<tauMin)
    tauMin = TMin[1];
  if(TMin[2]<tauMin)
    tauMin = TMin[2];
  if(TMin[3]<tauMin)
    tauMin = TMin[3];

  Taus[0] = "t11";
  Taus[1] = "t22";
  Taus[2] = "t33";
  Taus[3] = "t13";

  TauMax = TMax;
  TauMin = TMin;
  tauPos_x = new int[5];
  //Color scale position on screen 
  scale_xstart = 10.0;
  scale_xend = 310.0;
  scale_ystart = 40.0;
  scale_yend = 55.0;
 
  rangeButton_x = (int)scale_xstart + 10;
  rangeButton_y = (int)scale_ystart - 13;

  localValues = true;

  visual_field = 0;

  slider = 0.2;

  //tauLocalMax = max;
  //tauLocalMin = min;
  horiz = false;
  if(horiz)
    max_layer = nz;
  else
    max_layer = nx;

  
  getLocalTauValues();
	
  plane_layer = -1;

  //GLfloat data2[nz][ny][nx][4];

  GLfloat* data = new GLfloat[nx*ny*nz*4];
  
  for(int k=0; k < nz; k++){
    for(int i=0; i < ny; i++){
      for(int j = 0; j < nx; j++){
	
	  int idx = k*ny*nx + i*nx + j;
	
	  int texidx = k*ny*nx*4 + i*nx*4 + j*4;
	  /*
	  data2[k][i][j][0] =  tauData[idx].t11;
	  data2[k][i][j][1] =  tauData[idx].t22;
	  data2[k][i][j][2] =  tauData[idx].t33;
	  data2[k][i][j][3] =  tauData[idx].t13;*/

	  
	  data[texidx] = tauData[idx].t11;
	  data[texidx+1] = tauData[idx].t22;
	  data[texidx+2] = tauData[idx].t33;
	  data[texidx+3] = tauData[idx].t13;

      }
    }
  }

  glDisable(GL_TEXTURE_RECTANGLE_ARB);
  glEnable(texType);
  glGenTextures(1,tex_id);

  glBindTexture(texType, tex_id[0]);
  glTexParameteri(texType, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(texType, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(texType, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(texType, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(texType, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
	
  glTexImage3D(texType, 0, format, nx, ny, nz,0,GL_RGBA, GL_FLOAT,data); 
  CheckErrorsGL("3D texture\n");

  glBindTexture(texType, 0);

  glDisable(texType);

  delete [] data;

  plane_shader.addShader("Shaders/turbulencePlane_vp.glsl", GLSLObject::VERTEX_SHADER);
  plane_shader.addShader("Shaders/turbulencePlane_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  plane_shader.createProgram();
  u_tauTex = plane_shader.createUniform("Tau");
  u_max11 = plane_shader.createUniform("max11");
  u_max22 = plane_shader.createUniform("max22");
  u_max33 = plane_shader.createUniform("max33");
  u_max13 = plane_shader.createUniform("max13");
  u_min11 = plane_shader.createUniform("min11");
  u_min22 = plane_shader.createUniform("min22");
  u_min33 = plane_shader.createUniform("min33");
  u_min13 = plane_shader.createUniform("min13");
  u_controlTau = plane_shader.createUniform("controlTau");
  u_sliderTurb = plane_shader.createUniform("slider");
  
  glEnable(GL_TEXTURE_RECTANGLE_ARB);

  //Turbulence Color Scale
  scale_shader.addShader("Shaders/scale_vp.glsl", GLSLObject::VERTEX_SHADER);
  scale_shader.addShader("Shaders/scale_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  scale_shader.createProgram();
  uniform_xmax = scale_shader.createUniform("xmax");
  uniform_xmin = scale_shader.createUniform("xmin");
  uniform_tauMin = scale_shader.createUniform("tauMin");
  uniform_tauMax = scale_shader.createUniform("tauMax");
  uniform_sliderScale = scale_shader.createUniform("slider");
  

}
void VisualPlane::drawPlane(){
  
  glDisable(GL_TEXTURE_RECTANGLE_ARB);
  glDisable(GL_TEXTURE_2D);
  
  if(horiz)
    max_layer = nz;
  else
    max_layer = nx;
    
  if(visual_field > 0){

    if (plane_layer >= 0 && plane_layer < max_layer)
      {
	
	glPushMatrix();
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);     

	glEnable(texType);
	glBindTexture(texType, tex_id[0]);
	//glTexParameteri(texType, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	//glTexParameteri(texType, GL_TEXTURE_MIN_FILTER, GL_NEAREST);    

	//Horizontal Texture coordinates
	float s,t,r;
	
	if(horiz){
	  s = 0.0;
	  t = 0.0;
	  r = (float)plane_layer/((float)nz - 1.0);
	}
	//Vertical Texture coordinates
	else{
	  s = (float)plane_layer/((float)nx - 1.0);
	  t = 0.0;
	  r = 0.0;
	}
	plane_shader.activate();
      
	glUniform1iARB(u_controlTau, visual_field);
	glUniform1fARB(u_sliderTurb, slider);

	int tidx = plane_layer*4;
	if(localValues){
	  	  
	  glUniform1fARB(u_max11, tauLocalMax[tidx]);
	  glUniform1fARB(u_max22, tauLocalMax[tidx+1]);
	  glUniform1fARB(u_max33, tauLocalMax[tidx+2]);
	  glUniform1fARB(u_max13, tauLocalMax[tidx+3]);
  
	  glUniform1fARB(u_min11, tauLocalMin[tidx]);
	  glUniform1fARB(u_min22, tauLocalMin[tidx+1]);
	  glUniform1fARB(u_min33, tauLocalMin[tidx+2]);
	  glUniform1fARB(u_min13, tauLocalMin[tidx+3]);
	}
	else{
	  glUniform1fARB(u_max11, TauMax[0]);
	  glUniform1fARB(u_max22, TauMax[1]);
	  glUniform1fARB(u_max33, TauMax[2]);
	  glUniform1fARB(u_max13, TauMax[3]);
  
	  glUniform1fARB(u_min11, TauMin[0]);
	  glUniform1fARB(u_min22, TauMin[1]);
	  glUniform1fARB(u_min33, TauMin[2]);
	  glUniform1fARB(u_min13, TauMin[3]);
	}

	glUniform1iARB(u_tauTex, 0); 
	glColor4f(1.0,1.0,1.0,0.8);

	if(horiz){
	  glBegin(GL_QUADS);
	  {
	    glNormal3f(0.0, 0.0, 1.0);
	    glTexCoord3f( s, t, r);         glVertex3f(0, 0, plane_layer);
	    glTexCoord3f( s+1, t, r);       glVertex3f(nx, 0,plane_layer);
	    glTexCoord3f( s+1, t+1, r);     glVertex3f(nx, ny,plane_layer);
	    glTexCoord3f( s,  t+1, r);      glVertex3f(0, ny,plane_layer);
	  }
	  glEnd();
	}
	else{
	  glBegin(GL_QUADS);
	  glNormal3f(1.0,0.0,0.0);
	  glTexCoord3f( s, t, r);         glVertex3f(plane_layer, 0, 0);
	  glTexCoord3f( s, t, r+1);       glVertex3f(plane_layer, 0, nz);
	  glTexCoord3f( s, t+1, r+1);     glVertex3f(plane_layer, ny, nz);
	  glTexCoord3f( s, t+1, r);       glVertex3f(plane_layer, ny, 0);

	  glEnd();
	}

	plane_shader.deactivate();

	//glTexParameteri(texType, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	//glTexParameteri(texType, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glDisable(texType);
	glDisable(GL_BLEND);
	glDisable(GL_COLOR_MATERIAL);

	glPopMatrix();

      }
  }

  glEnable(GL_TEXTURE_RECTANGLE_ARB);

}
void VisualPlane::getLocalTauValues(){

  tauLocalMax = new float[4*max_layer];
  tauLocalMin = new float[4*max_layer];
 
  //Initialize max and min 
  for(int k=0; k<max_layer; k++){
    
    int idx;
    if(horiz){
      idx = k*ny*nx;
    }
    else{
      idx = k;
    }
    int tidx = k*4;

    tauLocalMax[tidx] = tau[idx].t11;
    tauLocalMax[tidx+1] = tau[idx].t22;
    tauLocalMax[tidx+2] = tau[idx].t33;
    tauLocalMax[tidx+3] = tau[idx].t13;

    
    tauLocalMin[tidx] = tau[idx].t11;
    tauLocalMin[tidx+1] = tau[idx].t22;
    tauLocalMin[tidx+2] = tau[idx].t33;
    tauLocalMin[tidx+3] = tau[idx].t13;
    
  }

  //Find max and min tau values for each height value

  for(int k=0; k<nz; k++){
    for(int i=0; i<ny; i++){
      for(int j=0; j<nx; j++){

	int idx = k*ny*nx + i*nx + j;
	int tidx;
	if(horiz){
	  //idx = k*ny*nx + i*nx + j;
	  tidx = k*4;
	}
	else{
	  tidx = j*4;	    
	}
	

	if(tau[idx].t11 > tauLocalMax[tidx])
	  tauLocalMax[tidx] = tau[idx].t11;
	if(tau[idx].t22 > tauLocalMax[tidx+1])
	  tauLocalMax[tidx+1] = tau[idx].t22;
	if(tau[idx].t33 > tauLocalMax[tidx+2])
	  tauLocalMax[tidx+2] = tau[idx].t33;
	if(tau[idx].t13 > tauLocalMax[tidx+3])
	  tauLocalMax[tidx+3] = tau[idx].t13;

	if(tau[idx].t11 < tauLocalMin[tidx])
	  tauLocalMin[tidx] = tau[idx].t11;
	if(tau[idx].t22 < tauLocalMin[tidx+1])
	  tauLocalMin[tidx+1] = tau[idx].t22;
	if(tau[idx].t33 < tauLocalMin[tidx+2])
	  tauLocalMin[tidx+2] = tau[idx].t33;
	if(tau[idx].t13 < tauLocalMin[tidx+3])
	  tauLocalMin[tidx+3] = tau[idx].t13;



      }
    }
  }
  

}
void VisualPlane::switchPlane(){
  horiz = !horiz;
  if(horiz){
    max_layer = nz;
  }
  else{
    max_layer = nx;
  }
  delete [] tauLocalMax;
  delete [] tauLocalMin;

  getLocalTauValues();

  if(plane_layer > max_layer)
    plane_layer = max_layer;

}

void VisualPlane::increasePlaneLayer(){
  plane_layer++;
  if(plane_layer > max_layer) plane_layer = max_layer;

}
void VisualPlane::decreasePlaneLayer(){
  plane_layer--;
  if(plane_layer < -1) plane_layer = -1;
}
void VisualPlane::drawScale(){
  GLint* vp = new GLint[4];
  glGetIntegerv(GL_VIEWPORT,vp);

  glDisable(GL_TEXTURE_RECTANGLE_ARB);
  glDisable(GL_DEPTH_TEST);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0, vp[2], 
	  0, vp[3], -1, 1);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
 

  //Draw Colored Scale
  char* disp;
  float max[4];
  float min[4];
  float globalMin = 0;
  float globalMax = 0;
  char* allDisp[4];

  /*if(visual_field == 0){
    for(int i=0; i < 4; i++){
      max[i] = WindMax[i];
      min[i] = WindMin[i];
      allDisp[i] = windTitle[i];
    }
    globalMin = windMin;
    globalMax = windMax;
    }*/
  if(visual_field != 0){
    int tidx = plane_layer*4;
    for(int i=0; i < 4; i++){
      
      if(localValues){
	max[i] = tauLocalMax[tidx+i];
	min[i] = tauLocalMin[tidx+i];
      }
      else{
	max[i] = TauMax[i];
	min[i] = TauMin[i];
      }
      allDisp[i] = Taus[i];
    }

    if(localValues){
      globalMax = tauLocalMax[tidx];
      if(tauLocalMax[tidx+1]>globalMax)
        globalMax = tauLocalMax[tidx+1];
      if(tauLocalMax[tidx+2]>globalMax)
	globalMax = tauLocalMax[tidx+2];
      if(tauLocalMax[tidx+3]>globalMax)
	globalMax = tauLocalMax[tidx+3];

      globalMin = tauLocalMin[tidx];
      if(tauLocalMin[tidx+1]<globalMin)
	globalMin = tauLocalMin[tidx+1];
      if(tauLocalMin[tidx+2]<globalMin)
	globalMin = tauLocalMin[tidx+2];
      if(tauLocalMin[tidx+3]<globalMin)
	globalMin = tauLocalMin[tidx+3];
    }
    else{
      globalMin = tauMin;
      globalMax = tauMax;
    }
  }

  switch(visual_field){
  case 0:
    disp = "Displaying: Wind Field";
    break;
  case 1:
    disp = "Displaying: Tau11";
    break;
  case 2:
    disp = "Displaying: Tau22";
    break;
  case 3:
    disp = "Displaying: Tau33";
    break;
  case 4:
    disp = "Displaying: Tau13";
    break;
  default:
    disp = "";
  }

  //Grey Background
  ////////////////////////////////////////////////////
  glColor3f(0.7,0.7,0.7);
  glBegin(GL_QUADS);
  {
    glVertex3f(scale_xstart-10,scale_ystart-18, 0.0);
    glVertex3f(scale_xend+35,scale_ystart-18, 0.0);
    glVertex3f(scale_xend+35,scale_yend+28, 0.0);
    glVertex3f(scale_xstart-10,scale_yend+28, 0.0);
  }
  glEnd();
  ////////////////////////////////////////////////////


  //Draw Colored scale
  ////////////////////////////////////////////////////
  scale_shader.activate();
  glUniform1fARB(uniform_xmin, scale_xstart);
  glUniform1fARB(uniform_xmax, scale_xend);
  glUniform1fARB(uniform_sliderScale, slider);
  //glUniform1fARB(uniform_tauMin, tauMin);
  //glUniform1fARB(uniform_tauMax, tauMax);
  glColor3f(1.0,1.0,1.0);

  glBegin(GL_QUADS);
  {
    glVertex3f(scale_xstart,scale_ystart, 0.0);
    glVertex3f(scale_xend,scale_ystart, 0.0);
    glVertex3f(scale_xend,scale_yend, 0.0);
    glVertex3f(scale_xstart,scale_yend, 0.0);
  }
  glEnd();

  scale_shader.deactivate();
  ///////////////////////////////////////////////////

  int y = (int)scale_yend + 5;

  //Display Text
  glColor3ub(255,255,0);
  glRasterPos2i((int)((scale_xstart+scale_xend)/2)- 50, (int)scale_ystart-13);
  for(int i=0; i < (int)strlen(disp); i++){
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, disp[i]);
  }

  if(visual_field != 0){

    if(localValues)
      disp = "(Local Range)";
    else
      disp = "(Global Range)";
    //Local or Global Display
    glColor3ub(255,255,0);
    glRasterPos2i(rangeButton_x, rangeButton_y);
    for(int i=0; i < (int)strlen(disp); i++){
      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, disp[i]);
    }

    sprintf(number, "%.2f", globalMin);
    glColor3ub(255, 255, 0);
    glRasterPos2i((int)scale_xstart-5, y);
    for(int i=0; i < (int)strlen(number); i++){
      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, number[i]);
    }
 
    //Display Max Values on Scale
    /////////////////////////////////////////////////////////////////
    for(int j=0; j <= 3; j++){
      tauPos_x[j] = (int)(((max[j]-globalMin)/(globalMax-globalMin))*(scale_xend-scale_xstart)+scale_xstart);
    }
    tauPos_x[4] = (int)scale_xstart;

    for(int j=0; j <= 3; j++){
      //map tau value onto x position of scale
      int x = (int)(((max[j]-globalMin)/(globalMax-globalMin))*(scale_xend-scale_xstart)+scale_xstart);

      int yPos = y;
        
      for(int i=j; i <=4; i++){
	if(abs(x-tauPos_x[i]) < 10 && abs(x-tauPos_x[i]) != 0)
	  yPos += 20;
      }

      sprintf(number, "%.2f", max[j]);
      glColor3ub(255,255,0);
      glRasterPos2i(x,yPos);
      for(int i=0; i < (int)strlen(number); i++)
	glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, number[i]); 

      glRasterPos2i(x,yPos+12);
      disp = allDisp[j];
      for(int i=0; i < (int)strlen(disp); i++)
	glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, disp[i]); 

    
      glColor3f(0.0,0.0,0.0);
      glBegin(GL_LINES);
      {
	glVertex3f(x,y,0.0);
	glVertex3f(x,scale_yend,0.0);
      }
      glEnd();

    }
    ////////////////////////////////////////////////////////////////

  }
  //Find position of slider bar and draw it on scale
  ////////////////////////////////////////////////////////////////
  slider_x = (int)((scale_xend-scale_xstart)*slider+scale_xstart);
  glColor3f(0.0,0.0,0.0);
  glBegin(GL_LINES);
  {
    glVertex3f(slider_x,scale_ystart,0.0);
    glVertex3f(slider_x,scale_yend,0.0);
  }
  glEnd();
  ////////////////////////////////////////////////////////////////

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_TEXTURE_RECTANGLE_ARB);

}
bool VisualPlane::clickedSlider(int x,int y){
  y = glutGet(GLUT_WINDOW_HEIGHT)-y;

  //std::cout << x << " " << y << std::endl;
  if((x < (slider_x + 5)) && (x > (slider_x - 5)) &&
     (y <= (int)scale_yend) && ( y >= (int)scale_ystart)){
    return true;
  }
  else return false;
}
void VisualPlane::moveSlider(int x){
  //slider_x = (int)((scale_xend-scale_xstart)*slider+scale_xstart);
  //slider_x = x;
  if((x <= scale_xend) && (x >= scale_xstart)) 
    slider = (float)(x - scale_xstart)/(float)(scale_xend-scale_xstart);
  
}
void VisualPlane::clickedRangeButton(int x, int y){
  y = glutGet(GLUT_WINDOW_HEIGHT)-y;
  
  if((x >= rangeButton_x) && ( x <= rangeButton_x + 80) &&
     (y >= rangeButton_y) && ( y <= rangeButton_y + 10)){
    localValues = !localValues;
  }
		     

}
void VisualPlane::moveSliderDown(){
  slider -= 0.01;
  if(slider < 0.0) slider = 0.0;
}
void VisualPlane::moveSliderUp(){
  slider += 0.01;
  if(slider > 1.0) slider = 1.0;
}
