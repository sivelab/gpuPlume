#include "VisualPlane.h"
#include "glErrorUtil.h"
#include <math.h>


static char number[128];
VisualPlane::VisualPlane(ParticleControl* pc, float* TMax,
			 float * TMin, float* max, float* min){

  texType = GL_TEXTURE_3D;
  format = GL_RGBA32F_ARB;

  tau = pc->tau;

  nx = pc->nx;
  ny = pc->ny;
  nz = pc->nz;
  nxdx = pc->nxdx;
  nydy = pc->nydy;
  nzdz = pc->nzdz;

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

  plane_normal = 0;
  if(plane_normal == 0)
    max_layer = nz;
  else if(plane_normal == 1)
    max_layer = nx;
  else
    max_layer = ny;

  
  getLocalTauValues();
	
  plane_layer = -1;
  plane_layer_z = -1.0;
  plane_layer_x = -1.0;
  plane_layer_y = -1.0;
  
  //GLfloat data2[nz][ny][nx][4];

  GLfloat* data = new GLfloat[nxdx*nydy*nzdz*4];
  
  for(int k=0; k < nzdz; k++){
    for(int i=0; i < nydy; i++){
      for(int j = 0; j < nxdx; j++){
	
	  int idx = k*nydy*nxdx + i*nxdx + j;
	
	  int texidx = k*nydy*nxdx*4 + i*nxdx*4 + j*4;
	  /*
	  data2[k][i][j][0] =  tau[idx].t11;
	  data2[k][i][j][1] =  tau[idx].t22;
	  data2[k][i][j][2] =  tau[idx].t33;
	  data2[k][i][j][3] =  tau[idx].t13;*/

	  
	  data[texidx] = tau[idx].t11;
	  data[texidx+1] = tau[idx].t22;
	  data[texidx+2] = tau[idx].t33;
	  data[texidx+3] = tau[idx].t13;

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
	
  glTexImage3D(texType, 0, format, nxdx, nydy, nzdz,0,GL_RGBA, GL_FLOAT,data); 
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
  

  //Initialize angle for rotating plane
  thetaX = M_PI_2;
  thetaX1 = 0.0;

  thetaY = M_PI_2;
  thetaY1 = 0.0;

  //Initialize Point in the domain
  r1.x = nx/2;
  r1.y = 0.0;
  r1.z = 0.0; 

  //Initialize normals and offsets of the domain planes
  n0.x = 1.0; n0.y = 0.0; n0.z = 0.0;
  n1.x = 0.0; n1.y = 1.0; n1.z = 0.0;
  n2.x = 0.0; n2.y = 0.0; n2.z = 1.0;
  n3.x = 0.0; n3.y = -1.0; n3.z = 0.0;
  n4.x = -1.0; n4.y = 0.0; n4.z = 0.0;
  n5.x = 0.0; n5.y = 0.0; n5.z = -1.0;
  
  d0 = 0.0;
  d1 = 0.0;
  d2 = 0.0;
  d3 = (ny);
  d4 = (nx);
  d5 = (nz);
  
  num_Points = 0;
  
  eps = 0.00001;


  yaw = 0.0;
  pitch = 0.0;
  roll = 0.0;

  n.x = 0.0;
  n.y = 0.0;
  n.z = 1.0;

  calculateNormal();
  
  /*
  //Determine normal of rotating plane
  n.x = cos(thetaX);
  n.y = cos(thetaY);
  n.z = sin(thetaX);
  
  if(fabs(n.x) < eps)
    n.x = 0.0;
  if(fabs(n.z) < eps)
    n.z = 0.0;
  if(fabs(n.y) < eps)
  n.y = 0.0;*/

  getIntersectionPoints();

}
void VisualPlane::drawPlane(){
  
  glDisable(GL_TEXTURE_RECTANGLE_ARB);
  glDisable(GL_TEXTURE_2D);
  
  if(plane_normal == 0){
    max_layer = nz;
    plane_layer = (int)plane_layer_z;
  }
  else if(plane_normal == 1){
    max_layer = nx;
    plane_layer = (int)plane_layer_x;
  }
  else{
    max_layer = ny;
    plane_layer = (int)plane_layer_y;
  }

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
	
	if(plane_normal == 0){
	  s = 0.0;
	  t = 0.0;
	  r = plane_layer_z/((float)nzdz - 1.0);
	}
	//Vertical Texture coordinates
	else if(plane_normal == 1){
	  s = plane_layer_x/((float)nxdx - 1.0);
	  t = 0.0;
	  r = 0.0;
	}
	else{
	  s = 0.0;
	  t = plane_layer_y/((float)nydy - 1.0);
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

	if(plane_normal == 0){
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
	else if(plane_normal == 1){
	  glBegin(GL_QUADS);
	  glNormal3f(1.0,0.0,0.0);
	  glTexCoord3f( s, t, r);         glVertex3f(plane_layer, 0, 0);
	  glTexCoord3f( s, t, r+1);       glVertex3f(plane_layer, 0, nz);
	  glTexCoord3f( s, t+1, r+1);     glVertex3f(plane_layer, ny, nz);
	  glTexCoord3f( s, t+1, r);       glVertex3f(plane_layer, ny, 0);

	  glEnd();
	}
	else{
	  glBegin(GL_QUADS);
	  glNormal3f(0.0,1.0,0.0);
	  glTexCoord3f( s, t, r);         glVertex3f(0, plane_layer, 0);
	  glTexCoord3f( s+1, t, r);       glVertex3f(nx, plane_layer, 0);
	  glTexCoord3f( s+1, t, r+1);     glVertex3f(nx, plane_layer, nz);
	  glTexCoord3f( s, t, r+1);       glVertex3f(0, plane_layer, nz);

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
    if(plane_normal == 0){
      idx = k*nydy*nxdx;
    }
    else if(plane_normal == 1){
      idx = k;
    }
    else{
      idx = k*nxdx;
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

  for(int k=0; k<nzdz; k++){
    for(int i=0; i<nydy; i++){
      for(int j=0; j<nxdx; j++){

	int idx = k*nydy*nxdx + i*nxdx + j;
	int tidx;
	if(plane_normal == 0){
	  //idx = k*ny*nx + i*nx + j;
	  tidx = k*4;
	}
	else if(plane_normal == 1){
	  tidx = j*4;	    
	}
	else{
	  tidx = i*4;
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
  plane_normal++;
  if(plane_normal > 2){
    plane_normal = 0;
  }

  if(plane_normal == 0){
    max_layer = nz;
  }
  else if(plane_normal == 1){
    max_layer = nx;
  }
  else
    max_layer = ny;

  delete [] tauLocalMax;
  delete [] tauLocalMin;

  getLocalTauValues();

  //if(plane_layer > max_layer)
  //plane_layer = max_layer;

}
void VisualPlane::increaseYaw(){
  yaw = M_PI_2/60.0;
 
  pitch = 0.0;
  roll = 0.0;

  calculateNormal();

  getIntersectionPoints();
}
void VisualPlane::increasePitch(){
  pitch = M_PI_2/60.0;

  yaw = 0.0;
  roll = 0.0;

  calculateNormal();

  getIntersectionPoints();
}
void VisualPlane::increaseRoll(){
  roll = M_PI_2/60.0;
 
  yaw = 0.0;
  pitch = 0.0;

  calculateNormal();
  getIntersectionPoints();
}
void VisualPlane::decreaseYaw(){
  yaw = -M_PI_2/60.0;

  pitch = 0.0;
  roll = 0.0;

  calculateNormal();

  getIntersectionPoints();
}
void VisualPlane::decreasePitch(){
  pitch = -M_PI_2/60.0;
 
  yaw = 0.0;
  roll = 0.0;

  calculateNormal();

  getIntersectionPoints();
}
void VisualPlane::decreaseRoll(){
  roll = -M_PI_2/60.0;
 
  yaw = 0.0;
  pitch = 0.0;

  calculateNormal();
  getIntersectionPoints();
  
}
void VisualPlane::calculateNormal(){

  float a11,a12,a13;
  float a21,a22,a23;
  float a31,a32,a33;
  
  //Transformation matrix
  a11 = cos(pitch)*cos(yaw);
  a12 = sin(roll)*sin(pitch)*cos(yaw) + cos(roll)*sin(yaw);
  a13 = (-cos(roll))*sin(pitch)*cos(yaw) + sin(roll)*sin(yaw);
  a21 = (-cos(pitch))*sin(yaw);
  a22 = (-sin(roll))*sin(pitch)*sin(yaw) + cos(roll)*cos(yaw);
  a23 = cos(roll)*sin(pitch)*sin(yaw) + sin(roll)*cos(yaw);
  a31 = sin(pitch);
  a32 = (-sin(roll))*cos(pitch);
  a33 = cos(roll)*cos(pitch);

  /*
  if((a11 < e) && (a11 > -e))
    a11 = 0.0;
  if((a12 < e) && (a12 > -e))
    a12 = 0.0;
  if((a13 < e) && (a13 > -e))
    a13 = 0.0;
  if((a21 < e) && (a21 > -e))
    a21 = 0.0;
  if((a22 < e) && (a22 > -e))
    a22 = 0.0;
  if((a23 < e) && (a23 > -e))
    a23 = 0.0;
  if((a31 < e) && (a31 > -e))
    a31 = 0.0;
  if((a32 < e) && (a32 > -e))
    a32 = 0.0;
  if((a33 < e) && (a33 > -e))
  a33 = 0.0;*/

  n.x = a11*n.x + a12*n.y + a13*n.z;
  n.y = a21*n.x + a22*n.y + a23*n.z;
  n.z = a31*n.x + a32*n.y + a33*n.z;
  
  /*
  float e = 0.1;
  if((n.x < e) && (n.x > -e))
    n.x = 0.0;
  if((n.z < e) && (n.z > -e))
    n.z = 0.0;
  if((n.y < e) && (n.y > -e))
    n.y = 0.0;
  */

  //normalize the normal
  float length;
  length = sqrt((n.x*n.x)+(n.y*n.y)+(n.z*n.z));
   
  n.x = n.x/length;
  n.y = n.y/length;
  n.z = n.z/length;
  

  std::cout << "normal is " << n.x << " " << n.y << " " << n.z << std::endl;

}

void VisualPlane::increaseAngle(){
  thetaX += M_PI_2/6.0;

  //if(thetaX > (M_PI+M_PI_2))
  //thetaX = M_PI+M_PI_2;

  thetaX1 += M_PI_2/6.0;
  
  //if(thetaX1 > M_PI)
  //thetaX1 = M_PI;
 
  //Determine normal of rotating plane
  n.x = cos(thetaX);
  n.y = cos(thetaY);
  n.z = sin(thetaX);
  
  if(fabs(n.x) < eps)
    n.x = 0.0;
  if(fabs(n.z) < eps)
    n.z = 0.0;
  if(fabs(n.y) < eps)
    n.y = 0.0;

  getIntersectionPoints();
}
void VisualPlane::decreaseAngle(){
  thetaX -= M_PI_2/6.0;

  //if(thetaX < M_PI_2)
  //thetaX = M_PI_2;

  thetaX1 -= M_PI_2/6.0;
  
  //if(thetaX1 < 0.0)
  //thetaX1 = 0.0;

  //Determine normal of rotating plane
  n.x = cos(thetaX);
  n.y = cos(thetaY);
  n.z = sin(thetaX);
  
  if(fabs(n.x) < eps)
    n.x = 0.0;
  if(fabs(n.z) < eps)
    n.z = 0.0;
  if(fabs(n.y) < eps)
    n.y = 0.0;

  getIntersectionPoints();
}
void VisualPlane::increaseAngle_Y(){
  thetaY += M_PI_2/6.0;

  //if(thetaY > (M_PI+M_PI_2))
  // thetaY = M_PI+M_PI_2;

  thetaY1 += M_PI_2/6.0;
  
  //if(thetaY1 > M_PI)
  //thetaY1 = M_PI;
 
  //Determine normal of rotating plane
  n.x = cos(thetaX);
  n.y = cos(thetaY);
  n.z = sin(thetaX);
  
  if(fabs(n.x) < eps)
    n.x = 0.0;
  if(fabs(n.z) < eps)
    n.z = 0.0;
  if(fabs(n.y) < eps)
    n.y = 0.0;

  getIntersectionPoints();
}
void VisualPlane::decreaseAngle_Y(){
  thetaY -= M_PI_2/6.0;

  //if(thetaY < M_PI_2)
  //thetaY = M_PI_2;

  thetaY1 -= M_PI_2/6.0;
  
  //if(thetaY1 < 0.0)
  //thetaY1 = 0.0;

  //Determine normal of rotating plane
  n.x = cos(thetaX);
  n.y = cos(thetaY);
  n.z = sin(thetaX);
  
  if(fabs(n.x) < eps)
    n.x = 0.0;
  if(fabs(n.z) < eps)
    n.z = 0.0;
  if(fabs(n.y) < eps)
    n.y = 0.0;

  getIntersectionPoints();
}
void VisualPlane::draw(){
  
  //glDisable(GL_TEXTURE_RECTANGLE_ARB);
  //glDisable(GL_TEXTURE_2D);
  /*
  if(plane_normal == 0){
    max_layer = nz;
    plane_layer = (int)plane_layer_z;
  }
  else if(plane_normal == 1){
    max_layer = nx;
    plane_layer = (int)plane_layer_x;
  }
  else{
    max_layer = ny;
    plane_layer = (int)plane_layer_y;
    }*/

  //getIntersectionPoints();

  if(visual_field > 0){

    //if (plane_layer >= 0 && plane_layer < max_layer)
    //{
	
	glPushMatrix();
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);     

	glEnable(texType);
	glBindTexture(texType, tex_id[0]);
	//glTexParameteri(texType, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	//glTexParameteri(texType, GL_TEXTURE_MIN_FILTER, GL_LINEAR);    

	
	plane_shader.activate();
      
	glUniform1iARB(u_controlTau, visual_field);
	glUniform1fARB(u_sliderTurb, slider);

	//int tidx = plane_layer*4;
	//if(localValues){
	  	  
	/*glUniform1fARB(u_max11, tauLocalMax[tidx]);
	  glUniform1fARB(u_max22, tauLocalMax[tidx+1]);
	  glUniform1fARB(u_max33, tauLocalMax[tidx+2]);
	  glUniform1fARB(u_max13, tauLocalMax[tidx+3]);
  
	  glUniform1fARB(u_min11, tauLocalMin[tidx]);
	  glUniform1fARB(u_min22, tauLocalMin[tidx+1]);
	  glUniform1fARB(u_min33, tauLocalMin[tidx+2]);
	  glUniform1fARB(u_min13, tauLocalMin[tidx+3]);*/
	  //}
	  //else{
	glUniform1fARB(u_max11, TauMax[0]);
	glUniform1fARB(u_max22, TauMax[1]);
	glUniform1fARB(u_max33, TauMax[2]);
	glUniform1fARB(u_max13, TauMax[3]);
  
	glUniform1fARB(u_min11, TauMin[0]);
	glUniform1fARB(u_min22, TauMin[1]);
	glUniform1fARB(u_min33, TauMin[2]);
	glUniform1fARB(u_min13, TauMin[3]);
	  //}

	glUniform1iARB(u_tauTex, 0); 
	glColor4f(1.0,1.0,1.0,0.8);

	glBegin(GL_QUADS);
	{
	  glNormal3f(n.x, n.y, n.z);

	  glTexCoord3f(p1.x/(float)(nxdx), p1.y/(float)(nydy), p1.z/(float)(nzdz-1));  
	  glVertex3f(p1.x, p1.y, p1.z);

	  glTexCoord3f(p2.x/(float)(nxdx), p2.y/(float)(nydy), p2.z/(float)(nzdz-1));  
	  glVertex3f(p2.x, p2.y, p2.z);

	  glTexCoord3f(p3.x/(float)(nxdx), p3.y/(float)(nydy), p3.z/(float)(nzdz-1));  
	  glVertex3f(p3.x, p3.y, p3.z);

	  glTexCoord3f(p0.x/(float)(nxdx), p0.y/(float)(nydy), p0.z/(float)(nzdz-1));  
	  glVertex3f(p0.x, p0.y, p0.z);
	}
	glEnd();
      	
	plane_shader.deactivate();

	glTexParameteri(texType, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(texType, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glDisable(texType);
	glDisable(GL_BLEND);
	glDisable(GL_COLOR_MATERIAL);

	glPopMatrix();

	//}
  }

  //Draw point of rotation on plane
  glPushMatrix();

  glColor3f(1.0,1.0,1.0);
 
  glBegin(GL_POINTS);
  {
    glVertex3f(p.x,p.y,p.z);
  }
  glEnd();


  glBegin(GL_LINES);
  {
    glVertex3f(p.x,p.y,p.z);
    glVertex3f(p.x+n.x,p.y+n.y,p.z+n.z);
    
  }
  glEnd();

  glPopMatrix();

  //glEnable(GL_TEXTURE_RECTANGLE_ARB);

}
void VisualPlane::findEdgePoints(){

  //vec3 p;
  
  float Z = (float)nz;
  float Y = (float)ny;
  float X = (float)nx;

  //p.x = 80; 
  //p.y = 25.0; 
  //p.z = 0.0;

  p.x = plane_layer_x; 
  p.y = plane_layer_y; 
  p.z = plane_layer_z;

  float z = 0.0;
  float z1 = 0.0;
  float y1 = 0.0;

  //First Point

  //Case: plane is parallel to z-y plane
  if(fabs(cos(thetaY1)) <= eps){
    e1.x = p.x;
    e1.y = p.y;
    e1.z = Z;
  }
  //Case: plane is parallel to x-y plane
  else if(fabs(sin(thetaY1)) <= eps){  
    e1.z = p.z;   
    e1.x = p.x;
    e1.y = 0.0;    
  }
  else{

    z = p.y*tan(thetaY1);
    
    if(z < 0.0)
      z1 = z + p.z;
    else if(z > Z)
      z1 = z - Z + p.z;
    else
      z1 = z;
 
    y1 = z1*(cos(thetaY1)/sin(thetaY1));  
    
    e1.x = p.x;
  
    if(y1 == p.y)
      e1.y = 0.0;
    else
      e1.y = y1;

    if(z >= Z)
      e1.z = Z;
    else if(z1 == z)
      e1.z = z + p.z;
    else
      e1.z = 0.0;

    //Corrections
    if(e1.z < 0.0)
      e1.z = 0.0;
    if(e1.z == 0.0 && e1.y == 0.0)
      e1.y = p.y;
  }
 
  std::cout << "e1.x = " << e1.x << " e1.y = " << e1.y << " e1.z = " << e1.z
      <<std::endl;

  //Second point
  
  //Case: plane is parallel to z-y plane
  if(fabs(cos(thetaY1)) <= eps){
    e2.x = p.x;
    e2.y = p.y;
    e2.z = 0.0;
  }
  //Case: plane is parallel to x-y plane
  else if(fabs(sin(thetaY1)) <= eps){  
    e2.z = p.z;   
    e2.x = p.x;
    e2.y = Y;    
  }
  else{
    z = (p.y-Y)*tan(thetaY1);
  
    if(z < 0.0)
      z1 = z + p.z;
    else if(z > Z)
      z1 = z - Z + p.z;
    else
      z1 = z;
    
    y1 = z1*(cos(thetaY1)/sin(thetaY1));
  
    e2.x = p.x;
      
    if(y1 == (p.y-Y))
      e2.y = Y;
    else
      e2.y = Y+y1;

    if(z >= Z)
      e2.z = Z;
    else if(z1 == z)
      e2.z = z + p.z;
    else
      e2.z = 0.0;
     
    //Corrections
    if(e2.z < 0.0)
      e2.z = 0.0;
    //This applies when starting point is on z == 0
    if(e2.z == 0.0 && e2.y == Y)
      e2.y = p.y;
  }
  std::cout << "e2.x = " << e2.x << " e2.y = " << e2.y << " e2.z = " << e2.z
    <<std::endl;
  
}
void VisualPlane::getIntersectionPoints(){
  /*
  //Determine normal of rotating plane
  n.x = cos(thetaX);
  n.y = cos(thetaY);
  n.z = sin(thetaX);
  
  if(fabs(n.x) < eps)
    n.x = 0.0;
  if(fabs(n.z) < eps)
    n.z = 0.0;
  if(fabs(n.y) < eps)
  n.y = 0.0;*/
    
  std::cout << "Plane's normal" << std::endl;
  std::cout << n.x << " " << n.y << " " << n.z << std::endl;

  //vec3 p; 
  float d;

  //Point on rotating plane
  p.x = plane_layer_x;
  p.y = plane_layer_y;
  p.z = plane_layer_z;

  //offset of rotating plane
  d = -dotProduct(p,n);
  //std::cout << "d = " << d << std::endl;

  //p, d, and n define the rotating plane

  vec3 a,b,c;
  float denom;
    
  //Make sure list is empty first
  pList.clear();
  num_Points = 0;

  //Now find out where all the plane intersection points are
  if( (dotProduct(n,(crossProduct(n0,n1)))) != 0.0){

    denom = dotProduct(n,(crossProduct(n0,n1)));

    a = crossProduct(n0,n1);
    b = crossProduct(n1,n);
    c = crossProduct(n,n0);
    
    p1.x = ((-d)*a.x)/denom;
    p1.y = ((-d)*a.y)/denom;
    p1.z = ((-d)*a.z)/denom;
    if(!(p1.x < 0.0 || p1.y < 0.0 || p1.z < 0.0 ||
	 p1.x > (nx) || p1.y > (ny) || p1.z > (nz))){
      pList.push_back(p1);
      num_Points++;

      std::cout << "0 and 1" << std::endl;
    }
  }
  if( ((dotProduct(n,(crossProduct(n0,n2)))) != 0.0) ) {
    denom = dotProduct(n,(crossProduct(n0,n2)));

    a = crossProduct(n0,n2);
    b = crossProduct(n2,n);
    c = crossProduct(n,n0);

    p1.x = ((-d)*a.x)/denom;
    p1.y = ((-d)*a.y)/denom;
    p1.z = ((-d)*a.z)/denom;
    if(!(p1.x < 0.0 || p1.y < 0.0 || p1.z < 0.0 ||
	 p1.x > (nx) || p1.y > (ny) || p1.z > (nz))){
      pList.push_back(p1);
      num_Points++;

      std::cout << "0 and 2" << std::endl;
    }
  }
  if( (dotProduct(n,(crossProduct(n1,n2))) != 0.0) ){
    
    denom = dotProduct(n,(crossProduct(n1,n2)));

    a = crossProduct(n1,n2);
    b = crossProduct(n2,n);
    c = crossProduct(n,n1);
       
    p1.x = ((-d)*a.x)/denom;
    p1.y = ((-d)*a.y)/denom;
    p1.z = ((-d)*a.z)/denom;
    if(!(p1.x < 0.0 || p1.y < 0.0 || p1.z < 0.0 ||
	 p1.x > (nx) || p1.y > (ny) || p1.z > (nz))){
      pList.push_back(p1);
      num_Points++;
   
      std::cout << "1 and 2" << std::endl;
    }
  }
  if( (dotProduct(n,(crossProduct(n1,n5))) != 0.0) ){
    
    denom = dotProduct(n,(crossProduct(n1,n5)));

    a = crossProduct(n1,n5);
    b = crossProduct(n5,n);
    c = crossProduct(n,n1);

    p1.x = (((-d)*a.x) - (d1*b.x) - (d5*c.x))/denom;
    p1.y = (((-d)*a.y) - (d1*b.y) - (d5*c.y))/denom;
    p1.z = (((-d)*a.z) - (d1*b.z) - (d5*c.z))/denom;
    if(!(p1.x < 0.0 || p1.y < 0.0 || p1.z < 0.0 ||
	 p1.x > (nx) || p1.y > (ny) || p1.z > (nz))){
      pList.push_back(p1);
      num_Points++;
      std::cout << "1 and 5" << std::endl;
    }
  }
  if( (dotProduct(n,(crossProduct(n0,n5))) != 0.0) ){
    
    denom = dotProduct(n,(crossProduct(n0,n5)));

    a = crossProduct(n0,n5);
    b = crossProduct(n5,n);
    c = crossProduct(n,n0);

    p1.x = (((-d)*a.x) - (d0*b.x) - (d5*c.x))/denom;
    p1.y = (((-d)*a.y) - (d0*b.y) - (d5*c.y))/denom;
    p1.z = (((-d)*a.z) - (d0*b.z) - (d5*c.z))/denom;
    if(!(p1.x < 0.0 || p1.y < 0.0 || p1.z < 0.0 ||
	 p1.x > (nx) || p1.y > (ny) || p1.z > (nz))){ 
      pList.push_back(p1);
      num_Points++;

      std::cout << "0 and 5" << std::endl;
    }
  }
  if( (dotProduct(n,(crossProduct(n0,n3))) != 0.0) ){
    
    denom = dotProduct(n,(crossProduct(n0,n3)));

    a = crossProduct(n0,n3);
    b = crossProduct(n3,n);
    c = crossProduct(n,n0);

    p1.x = (((-d)*a.x) - (d0*b.x) - (d3*c.x))/denom;
    p1.y = (((-d)*a.y) - (d0*b.y) - (d3*c.y))/denom;
    p1.z = (((-d)*a.z) - (d0*b.z) - (d3*c.z))/denom;
    if(!(p1.x < 0.0 || p1.y < 0.0 || p1.z < 0.0 ||
	 p1.x > (nx) || p1.y > (ny) || p1.z > (nz))){
      pList.push_back(p1);
      num_Points++;
   
      std::cout << "0 and 3" << std::endl;
    }
  }
  if( (dotProduct(n,(crossProduct(n3,n5))) != 0.0) ){
    
    denom = dotProduct(n,(crossProduct(n3,n5)));

    a = crossProduct(n3,n5);
    b = crossProduct(n5,n);
    c = crossProduct(n,n3);

    p1.x = (((-d)*a.x) - (d3*b.x) - (d5*c.x))/denom;
    p1.y = (((-d)*a.y) - (d3*b.y) - (d5*c.y))/denom;
    p1.z = (((-d)*a.z) - (d3*b.z) - (d5*c.z))/denom;
    if(!(p1.x < 0.0 || p1.y < 0.0 || p1.z < 0.0 ||
	 p1.x > (nx) || p1.y > (ny) || p1.z > (nz))){
      pList.push_back(p1);
      num_Points++;
    
      std::cout << "3 and 5" << std::endl;
    }
    
  }
  if( (dotProduct(n,(crossProduct(n3,n2))) != 0.0) ){
    
    denom = dotProduct(n,(crossProduct(n3,n2)));

    a = crossProduct(n3,n2);
    b = crossProduct(n2,n);
    c = crossProduct(n,n3);

    p1.x = (((-d)*a.x) - (d3*b.x) - (d2*c.x))/denom;
    p1.y = (((-d)*a.y) - (d3*b.y) - (d2*c.y))/denom;
    p1.z = (((-d)*a.z) - (d3*b.z) - (d2*c.z))/denom;
    if(!(p1.x < 0.0 || p1.y < 0.0 || p1.z < 0.0 ||
	 p1.x > (nx) || p1.y > (ny) || p1.z > (nz))){
      pList.push_back(p1);
      num_Points++;
      std::cout << "3 and 2" << std::endl;
    }
  }
  if( (dotProduct(n,(crossProduct(n3,n4))) != 0.0) ){
    
    denom = dotProduct(n,(crossProduct(n3,n4)));

    a = crossProduct(n3,n4);
    b = crossProduct(n4,n);
    c = crossProduct(n,n3);

    p1.x = (((-d)*a.x) - (d3*b.x) - (d4*c.x))/denom;
    p1.y = (((-d)*a.y) - (d3*b.y) - (d4*c.y))/denom;
    p1.z = (((-d)*a.z) - (d3*b.z) - (d4*c.z))/denom;
    if(!(p1.x < 0.0 || p1.y < 0.0 || p1.z < 0.0 ||
	 p1.x > (nx) || p1.y > (ny) || p1.z > (nz))){
      pList.push_back(p1);
      num_Points++;
      std::cout << "3 and 4" << std::endl;
    }
  }
  if( (dotProduct(n,(crossProduct(n5,n4))) != 0.0) ){
    
    denom = dotProduct(n,(crossProduct(n5,n4)));

    a = crossProduct(n5,n4);
    b = crossProduct(n4,n);
    c = crossProduct(n,n5);

    p1.x = (((-d)*a.x) - (d5*b.x) - (d4*c.x))/denom;
    p1.y = (((-d)*a.y) - (d5*b.y) - (d4*c.y))/denom;
    p1.z = (((-d)*a.z) - (d5*b.z) - (d4*c.z))/denom;
        
    if(!(p1.x < 0.0 || p1.y < 0.0 || p1.z < 0.0 ||
	 p1.x > (nx) || p1.y > (ny) || p1.z > (nz))){
      pList.push_back(p1);
      num_Points++;
      std::cout << "5 and 4" << std::endl;
    }
  }
  if( (dotProduct(n,(crossProduct(n2,n4))) != 0.0) ){
    
    denom = dotProduct(n,(crossProduct(n2,n4)));

    a = crossProduct(n2,n4);
    b = crossProduct(n4,n);
    c = crossProduct(n,n2);

    p1.x = (((-d)*a.x) - (d2*b.x) - (d4*c.x))/denom;
    p1.y = (((-d)*a.y) - (d2*b.y) - (d4*c.y))/denom;
    p1.z = (((-d)*a.z) - (d2*b.z) - (d4*c.z))/denom;
    if(!(p1.x < 0.0 || p1.y < 0.0 || p1.z < 0.0 ||
	 p1.x > (nx) || p1.y > (ny) || p1.z > (nz))){ 
      pList.push_back(p1);
      num_Points++;
      std::cout << "2 and 4" << std::endl;

    }
  }
  if( (dotProduct(n,(crossProduct(n1,n4))) != 0.0) ){
    
    denom = dotProduct(n,(crossProduct(n1,n4)));

    a = crossProduct(n1,n4);
    b = crossProduct(n4,n);
    c = crossProduct(n,n1);

    p1.x = (((-d)*a.x) - (d1*b.x) - (d4*c.x))/denom;
    p1.y = (((-d)*a.y) - (d1*b.y) - (d4*c.y))/denom;
    p1.z = (((-d)*a.z) - (d1*b.z) - (d4*c.z))/denom;
        
    if(!(p1.x < 0.0 || p1.y < 0.0 || p1.z < 0.0 ||
	 p1.x > (nx) || p1.y > (ny) || p1.z > (nz))){

      pList.push_back(p1);
      num_Points++;
      std::cout << "1 and 4" << std::endl;

    }
  }

  //Now we have all the intersection points, but some 
  //will be outside the domain.  
  //Need to get rid off the points outside the domain

  listIter = pList.begin();
  vec3 t;

  
  int n = 0;

  while(listIter != pList.end()){
    t = *listIter;
        
    if(t.x < 0.0 || t.y < 0.0 || t.z < 0.0 ||
       t.x > (nx) || t.y > (ny) || t.z > (nz)){
      
      //pList.pop_back();
      num_Points--;
    }
    else{
      if(t.x == -0.0)
	t.x = 0.0;
      if(t.y == -0.0)
	t.y = 0.0;
      if(t.z == -0.0)
	t.z = 0.0;
      
      switch(n){
      case 0:
	p1.x = t.x;
	p1.y = t.y;
	p1.z = t.z;
	break;
      case 1:
	p2.x = t.x;
	p2.y = t.y;
	p2.z = t.z;
	break;
      case 2:
	p3.x = t.x;
	p3.y = t.y;
	p3.z = t.z;
	break;
      case 3:
	p0.x = t.x;
	p0.y = t.y;
	p0.z = t.z;
	break;
	
      }
      n++;

      //piList.push_back(t);
      std::cout << "t.x = " << t.x << " t.y = " << 
	t.y << " t.z = " << t.z <<std::endl;
      
    }
    listIter++;
  }

  std::cout << "Num points = " << num_Points << std::endl;

  //Now we have the four intersection points of the domain
  //They need to be sorted so that they can be specified in
  //counter-clockwise order for drawing the plane
  
  //Might not have to do this???
  //sortIntersectionPoints();
  
}
void VisualPlane::sortIntersectionPoints(){
  


}

/*
void VisualPlane::getIntersectionPoints(){
  vec3 p;
  
  float eps = 0.00001;

  //Determine normal of plane
  n.x = cos(thetaX);
  n.y = cos(thetaY);
  n.z = sin(thetaX);

  if(fabs(n.x) < eps)
    n.x = 0.0;
  if(fabs(n.z) < eps)
    n.z = 0.0;
  if(fabs(n.y) < eps)
    n.y = 0.0;
  

  //Find the two points on the edge of planes y = 0 and y = Y
  findEdgePoints();
  //Using these two points, we can now get the intersection points of the domain
  float z = 0.0;
  float z1 = 0.0;
  float y1 = 0.0;
  float Z = (float)nz;
  float Y = (float)ny;
  float X = (float)nx;

  p.x = e1.x; 
  p.y = e1.y; 
  p.z = e1.z;
  
  //First Point

  //Case: plane is parallel to z-y plane
  if(fabs(cos(thetaX1)) <= eps){
    p2.x = p.x;
    p2.y = p.y;
    p2.z = Z;
  }
  //Case: plane is parallel to x-y plane
  else if(fabs(sin(thetaX1)) <= eps){  
    p2.z = p.z;   
    p2.x = X;
    p2.y = p.y;    
  }
  else{

    z = (X-p.x)*tan(thetaX1);
       
    if(z < 0.0)
      z1 = z + p.z;
    else if(z > Z)
      z1 = z - Z + p.z;
    else
      z1 = z;
 
    y1 = z1*(cos(thetaX1)/sin(thetaX1));  
    
    p2.y = p.y;
  
    if(y1 == (X-p.x))
      p2.x = X;
    else
      p2.x = X-y1;
    
    if(z >= Z)
      p2.z = Z;
    else if(z1 == z)
      p2.z = z + p.z;
    else
      p2.z = 0.0;

    //Corrections
    if(p2.z < 0.0)
      p2.z = 0.0;
    if(p2.z == 0.0 && p2.x == X)
      p2.x = p.x;
  }
 
  std::cout << "p2.x = " << p2.x << " p2.y = " << p2.y << " p2.z = " << p2.z
	    <<std::endl;
  
  //Case: plane is parallel to z-y plane
  if(fabs(cos(thetaX1)) <= eps){
    p1.x = p.x;
    p1.y = p.y;
    p1.z = 0.0;
  }
  //Case: plane is parallel to x-y plane
  else if(fabs(sin(thetaX1)) <= eps){  
    p1.z = p.z;   
    p1.x = 0.0;
    p1.y = p.y;    
  }
  else{

    z = -p.x*tan(thetaX1);
       
    if(z < 0.0)
      z1 = z + p.z;
    else if(z > Z)
      z1 = z - Z + p.z;
    else
      z1 = z;
 
    y1 = z1*(cos(thetaX1)/sin(thetaX1));  
    
    p1.y = p.y;
  
    if(y1 == -p.x)
      p1.x = 0.0;
    else
      p1.x = -y1;

    if(z >= (float)Z)
      p1.z = Z;
    else if(z1 == z)
      p1.z = z + p.z;
    else
      p1.z = 0.0;

    //Corrections
    if(p1.z < 0.0)
      p1.z = 0.0;
    if(p1.z == 0.0 && p1.x == 0.0)
      p1.x = p.x;
  }
 
  std::cout << "p1.x = " << p1.x << " p1.y = " << p1.y << " p1.z = " << p1.z
	    <<std::endl;

  ///////////////////////////////////////////////////////////////////////////
  // Now find p0 and p3
  ///////////////////////////////////////////////////////////////////////////
  //Use e2 as the starting point.
  p.x = e2.x; 
  p.y = e2.y; 
  p.z = e2.z;
  
  //First Point
  
  //Case: plane is parallel to z-y plane
  if(fabs(cos(thetaX1)) <= eps){
    p3.x = p.x;
    p3.y = p.y;
    p3.z = Z;
  }
  //Case: plane is parallel to x-y plane
  else if(fabs(sin(thetaX1)) <= eps){  
    p3.z = p.z;   
    p3.x = X;
    p3.y = p.y;    
  }
  else{

    z = (X-p.x)*tan(thetaX1);
       
    if(z < 0.0)
      z1 = z + p.z;
    else if(z > Z)
      z1 = z - Z + p.z;
    else
      z1 = z;
 
    y1 = z1*(cos(thetaX1)/sin(thetaX1));  
    
    p3.y = p.y;
  
    if(y1 == (X-p.x))
      p3.x = X;
    else
      p3.x = X-y1;
    
    if(z >= (float)Z)
      p3.z = Z;
    else if(z1 == z)
      p3.z = z + p.z;
    else
      p3.z = 0.0;

    //Corrections
    if(p3.z < 0.0)
      p3.z = 0.0;
    if(p3.z == 0.0 && p3.x == X)
      p3.x = p.x;
  }
 
  std::cout << "p3.x = " << p3.x << " p3.y = " << p3.y << " p3.z = " << p3.z
	    <<std::endl;
  
  //Case: plane is parallel to z-y plane
  if(fabs(cos(thetaX1)) <= eps){
    p0.x = p.x;
    p0.y = p.y;
    p0.z = 0.0;
  }
  //Case: plane is parallel to x-y plane
  else if(fabs(sin(thetaX1)) <= eps){  
    p0.z = p.z;   
    p0.x = 0.0;
    p0.y = p.y;    
  }
  else{

    z = -p.x*tan(thetaX1);
       
    if(z < 0.0)
      z1 = z + p.z;
    else if(z > Z)
      z1 = z - Z + p.z;
    else
      z1 = z;
 
    y1 = z1*(cos(thetaX1)/sin(thetaX1));  
    
    p0.y = p.y;
  
    if(y1 == -p.x)
      p0.x = 0.0;
    else
      p0.x = -y1;

    if(z >= Z)
      p0.z = nz;
    else if(z1 == z)
      p0.z = z + p.z;
    else
      p0.z = 0.0;

    //Corrections
    if(p0.z < 0.0)
      p0.z = 0.0;
    if(p0.z == 0.0 && p0.x == 0.0)
      p0.x = p.x;
  }
 
  std::cout << "p0.x = " << p0.x << " p0.y = " << p0.y << " p0.z = " << p0.z
	    <<std::endl;
  
  std::cout << std::endl;
  }*/

/*
void VisualPlane::getIntersectionPoints(){
  
  vec3 PI;
  PI.x = 0.0; PI.y = 0.0; PI.z = 0.0;
  float s1,s2,s3,s4;
  vec3 u;
  vec3 w;

  float eps = 0.0001;
  
  //Determine normal of plane
  n.x = cos(thetaX);
  n.y = cos(thetaY);
  n.z = sin(thetaX);

  if(fabs(n.x) < eps)
    n.x = 0.0;
  if(fabs(n.z) < eps)
    n.z = 0.0;
  if(fabs(n.y) < eps)
    n.y = 0.0;
  
  //Find the two points on the edge of planes y = 0 and y = ny
  findEdgePoints();
  //Using these two points, we can now get the intersection points of the domain

  //Point on the domain edge
  r1.x = e1.x;
  r1.y = e1.y;
  r1.z = e1.z; 
  //Find Second point to make line
  r2.x = r1.x + (2*nx*cos(thetaX1));
  r2.y = r1.y;
  r2.z = r1.z + (2*nx*sin(thetaX1));
 
  if(fabs(r2.x) < eps) 
    r2.x = 0.0;
  if(fabs(r2.z) < eps)
    r2.z = 0.0;
  if(fabs(r2.y) < eps)
    r2.y = 0.0;

  std::cout << "r2.x = " << r2.x << " r2.y = " << r2.y 
	    << " r2.z = " << r2.z << std::endl;
  //Now find intersection Point on plane

  //difference between two points
  u.x = r2.x-r1.x; u.y = r2.y-r1.y; u.z = r2.z-r1.z;
  //1st Plane

  //point on domain plane
  V.x = nx; V.y = 0.0; V.z = 0.0;
  w.x = r1.x-V.x; w.y = r1.y-V.y; w.z = r1.z-V.z;
  //normal of domain plane
  N.x = -1.0; N.y = 0.0; N.z = 0.0; 
  s1 = ((-N.x)*w.x + (-N.y)*w.y + (-N.z)*w.z)/(N.x*u.x + N.y*u.y + N.z*u.z);
 
  //2nd Plane
  //V stays the same for this plane
  N.x = 0.0; N.y = 0.0; N.z = 1.0;
  s2 = ((-N.x)*w.x + (-N.y)*w.y + (-N.z)*w.z)/(N.x*u.x + N.y*u.y + N.z*u.z);
 
  //3rd Plane
  V.x = 0.0; V.y = 0.0; V.z = nz;
  N.z = -1.0;
  w.x = r1.x-V.x; w.y = r1.y-V.y; w.z = r1.z-V.z;
  s3 = ((-N.x)*w.x + (-N.y)*w.y + (-N.z)*w.z)/(N.x*u.x + N.y*u.y + N.z*u.z);
 
  //4th Plane
  V.x = 0.0; V.y = 0.0; V.z = 0.0;
  w.x = r1.x-V.x; w.y = r1.y-V.y; w.z = r1.z-V.z;
  //normal of domain plane
  N.x = 1.0; N.y = 0.0; N.z = 0.0;
  s4 = ((-N.x)*w.x + (-N.y)*w.y + (-N.z)*w.z)/(N.x*u.x + N.y*u.y + N.z*u.z);
 

  //Find Smallest s - this means that intersection is the one to use
  float s = 100.0;
  if(r2.z <= 0.0){    
    if((s1 < s) && (s1 >= 0.0)){
      s = s1;
    }
    if((s2 < s) && (s2 >= 0.0)){
      s = s2;
    }
    if((s3 < s) && (s3 >= 0.0)){
      s = s3;
    }
    if((s4 < s) && (s4 >= 0.0)){
      s = s4;
    }
  }
  else{
    if((s1 < s) && (s1 > 0.0)){
      s = s1;
    }
    if((s2 < s) && (s2 > 0.0)){
      s = s2;
    }
    if((s3 < s) && (s3 > 0.0)){
      s = s3;
    }
    if((s4 < s) && (s4 > 0.0)){
      s = s4;
    }
  }
  
  //Find point of intersection
  PI.x = s*u.x + r1.x;
  PI.y = s*u.y + r1.y;
  PI.z = s*u.z + r1.z;
  
  std::cout << "PI.x = " << PI.x << " PI.y = " << PI.y << " PI.z = " << PI.z
	    <<std::endl;

  //We now have the point of Intersection, which is p2
  p2.x = PI.x; p2.y = PI.y; p2.z = PI.z;
  p3.x = PI.x; p3.y = PI.y+ny; p3.z = PI.z;
  

  //Now find p0 and p1
  //Need to change second point on plane
  r2.x = r1.x - (2*nx*cos(thetaX1));
  r2.y = r1.y;
  r2.z = r1.z - (2*nx*sin(thetaX1));
  if(fabs(r2.x) < eps) 
    r2.x = 0.0;
  if(fabs(r2.z) < eps)
    r2.z = 0.0;
  if(fabs(r2.y) < eps)
    r2.y = 0.0;
  
  std::cout << "r2.x = " << r2.x << " r2.y = " << r2.y 
	    << " r2.z = " << r2.z << std::endl;
  //point on domain plane
  V.x = 0.0; V.y = 0.0; V.z = 0.0;
  w.x = r1.x-V.x; w.y = r1.y-V.y; w.z = r1.z-V.z;
  //normal of domain plane
  N.x = 1.0; N.y = 0.0; N.z = 0.0;
  u.x = r2.x-r1.x; u.y = r2.y-r1.y; u.z = r2.z-r1.z;
  //s = dot(-N,w)/dot(N,u)
  s1 = ((-N.x)*w.x + (-N.y)*w.y + (-N.z)*w.z)/(N.x*u.x + N.y*u.y + N.z*u.z);
 
  //2nd Plane
  //V and w stay the same for this plane
  N.x = 0.0; N.y = 0.0; N.z = 1.0;
  s2 = ((-N.x)*w.x + (-N.y)*w.y + (-N.z)*w.z)/(N.x*u.x + N.y*u.y + N.z*u.z);
 
  //3rd Plane
  V.x = 0.0; V.y = 0.0; V.z = nz;
  w.x = r1.x-V.x; w.y = r1.y-V.y; w.z = r1.z-V.z;
  N.z = -1.0;
  s3 = ((-N.x)*w.x + (-N.y)*w.y + (-N.z)*w.z)/(N.x*u.x + N.y*u.y + N.z*u.z);
  
  //4th Plane
  V.x = nx; V.y = 0.0; V.z = 0.0;
  w.x = r1.x-V.x; w.y = r1.y-V.y; w.z = r1.z-V.z;
  //normal of domain plane
  N.x = -1.0; N.y = 0.0; N.z = 0.0;
  s4 = ((-N.x)*w.x + (-N.y)*w.y + (-N.z)*w.z)/(N.x*u.x + N.y*u.y + N.z*u.z);
 
  //Find Smallest s - this means that intersection is the one to use
  
  s = 100.0;
  if(r2.z <= 0.0){    
    if((s1 < s) && (s1 >= 0.0)){
      s = s1;
    }
    if((s2 < s) && (s2 >= 0.0)){
      s = s2;
    }
    if((s3 < s) && (s3 >= 0.0)){
      s = s3;
    }
    if((s4 < s) && (s4 >= 0.0)){
      s = s4;
    }
  }
  else{
    if((s1 < s) && (s1 > 0.0)){
      s = s1;
    }
    if((s2 < s) && (s2 > 0.0)){
      s = s2;
    }
    if((s3 < s) && (s3 > 0.0)){
      s = s3;
    }
    if((s4 < s) && (s4 > 0.0)){
      s = s4;
    }
  }
  //Find point of intersection
  PI.x = s*u.x + r1.x;
  PI.y = s*u.y + r1.y;
  PI.z = s*u.z + r1.z;

  //We now have the point of Intersection and can find p0 and p1
  p0.x = PI.x; p0.y = PI.y+ny; p0.z = PI.z;
  p1.x = PI.x; p1.y = PI.y; p1.z = PI.z;
  
  std::cout << p0.x << " " << p0.y << " " << p0.z << std::endl;
  std::cout << p1.x << " " << p1.y << " " << p1.z << std::endl;
  std::cout << p2.x << " " << p2.y << " " << p2.z << std::endl;
  std::cout << p3.x << " " << p3.y << " " << p3.z << std::endl;
  }*/
void VisualPlane::increasePlaneLayer(){
  if(plane_normal == 0){
    plane_layer_z += 1.0;
    if(plane_layer_z > (float)max_layer)
      plane_layer_z = (float)max_layer;
  }
  else if(plane_normal == 1){
    plane_layer_x += 1.0;
    if(plane_layer_x > (float)max_layer)
      plane_layer_x = (float)max_layer;
  }
  else{
    plane_layer_y += 1.0;
    if(plane_layer_y > (float)max_layer)
      plane_layer_y = (float)max_layer;
  }

  getIntersectionPoints();
  //plane_layer++;
  //if(plane_layer > max_layer) plane_layer = max_layer;

}
void VisualPlane::decreasePlaneLayer(){
  if(plane_normal == 0){
    plane_layer_z -= 1.0;
    if(plane_layer_z < -1.0)
      plane_layer_z = -1.0;
  }
  else if(plane_normal == 1){
    plane_layer_x -= 1.0;
    if(plane_layer_x < -1.0)
      plane_layer_x = -1.0;
  }
  else{
    plane_layer_y -= 1.0;
    if(plane_layer_y < -1.0)
      plane_layer_y = -1.0;
  }
  
  getIntersectionPoints();
  //plane_layer--;
  //if(plane_layer < -1) plane_layer = -1;
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
/*void VisualPlane::clickedDisplayTauButton(int x, int y){
  y = glutGet(GLUT_WINDOW_HEIGHT) - y;

  

  }*/
void VisualPlane::moveSliderDown(){
  slider -= 0.01;
  if(slider < 0.0) slider = 0.0;
}
void VisualPlane::moveSliderUp(){
  slider += 0.01;
  if(slider > 1.0) slider = 1.0;
}
vec3 VisualPlane::crossProduct(vec3 a, vec3 b){
  vec3 c;
  
  c.x = (a.y*b.z)-(a.z*b.y);
  c.y = (a.z*b.x)-(a.x*b.z);
  c.z = (a.x*b.y)-(a.y*b.x);

  return c;  
}
float VisualPlane::dotProduct(vec3 a, vec3 b){
  float c;

  c = (a.x*b.x) + (a.y*b.y) + (a.z*b.z);
  
  return c;
  
}

/*
void VisualPlane::findEdgePoints(){
  vec3 u;
  vec3 w;
  
  float s1,s2,s3,s4;
  float eps = 0.0001;

  vec3 p;
  p.x = plane_layer_x; 
  p.y = 0.0; 
  p.z = plane_layer_z;
  vec3 p1;
  p1.x = p.x; 
  p1.y = p.y - 2*ny*cos(thetaY1);
  p1.z = p.z + 2*ny*sin(thetaY1); 
  
  if(fabs(p1.z) < eps)
    p1.z = 0.0;
  if(fabs(p1.y) < eps)
    p1.y = 0.0;

  //difference between two points
  u.x = p1.x-p.x; u.y = p1.y-p.y; u.z = p1.z-p.z;
  //1st Plane

  //point on domain plane
  V.x = 0.0; V.y = 0.0; V.z = 0.0;
  w.x = p.x-V.x; w.y = p.y-V.y; w.z = p.z-V.z;
  //normal of domain plane
  N.x = 0.0; N.y = 1.0; N.z = 0.0; 
  s1 = ((-N.x)*w.x + (-N.y)*w.y + (-N.z)*w.z)/(N.x*u.x + N.y*u.y + N.z*u.z);
 
  //2nd Plane
  //V stays the same for this plane
  N.x = 0.0; N.y = 0.0; N.z = 1.0;
  s2 = ((-N.x)*w.x + (-N.y)*w.y + (-N.z)*w.z)/(N.x*u.x + N.y*u.y + N.z*u.z);
 
  //3rd Plane
  V.x = 0.0; V.y = 0.0; V.z = nz;
  N.z = -1.0;
  w.x = p.x-V.x; w.y = p.y-V.y; w.z = p.z-V.z;
  s3 = ((-N.x)*w.x + (-N.y)*w.y + (-N.z)*w.z)/(N.x*u.x + N.y*u.y + N.z*u.z);
 
  //4th Plane
  V.x = 0.0; V.y = ny; V.z = 0.0;
  w.x = p.x-V.x; w.y = p.y-V.y; w.z = p.z-V.z;
  //normal of domain plane
  N.x = 0.0; N.y = -1.0; N.z = 0.0;
  s4 = ((-N.x)*w.x + (-N.y)*w.y + (-N.z)*w.z)/(N.x*u.x + N.y*u.y + N.z*u.z);
 

  //Find Smallest s - this means that intersection is the one to use
  float s = 100.0;
    
  //+y plane
  if((s1 < s) && (s1 >= 0.0)){
    s = s1;
  }
  //+z plane
  if((s2 < s) && (s2 >= 0.0) && (p.z < p1.z)){
    s = s2;
    //std::cout << "s2" << std::endl;
  }
  //-z plane
  if((s3 < s) && (s3 >= 0.0)){
    s = s3;
  }
  //-y plane
  if((s4 < s) && (s4 >= 0.0) && (p1.y <= 0.0)){
    s = s4;
  }
 
  s=s1;
    
  //Find point of intersection, which is e1
  e1.x = s*u.x + p.x;
  e1.y = s*u.y + p.y;
  e1.z = s*u.z + p.z;
  
  //std::cout << "s = " << s << " and u.z =  " << u.z <<std::endl;

  std::cout << "e1.x = " << e1.x << " e1.y = " << e1.y << " e1.z = " << e1.z
	    <<std::endl;

  //Now find e2, which is the point on the opposite edge of plane
  p1.x = p.x; 
  p1.y = p.y + 2*ny*cos(thetaY1);
  p1.z = p.z - 2*ny*sin(thetaY1); 
  
  if(fabs(p1.z) < eps)
    p1.z = 0.0;
  if(fabs(p1.y) < eps)
    p1.y = 0.0;

  //difference between two points
  u.x = p1.x-p.x; u.y = p1.y-p.y; u.z = p1.z-p.z;
  //1st Plane

  //point on domain plane
  V.x = 0.0; V.y = 0.0; V.z = 0.0;
  w.x = p.x-V.x; w.y = p.y-V.y; w.z = p.z-V.z;
  //normal of domain plane
  N.x = 0.0; N.y = 1.0; N.z = 0.0; 
  s1 = ((-N.x)*w.x + (-N.y)*w.y + (-N.z)*w.z)/(N.x*u.x + N.y*u.y + N.z*u.z);
 
  //2nd Plane
  //V stays the same for this plane
  N.x = 0.0; N.y = 0.0; N.z = 1.0;
  s2 = ((-N.x)*w.x + (-N.y)*w.y + (-N.z)*w.z)/(N.x*u.x + N.y*u.y + N.z*u.z);
 
  //3rd Plane
  V.x = 0.0; V.y = 0.0; V.z = nz;
  N.z = -1.0;
  w.x = p.x-V.x; w.y = p.y-V.y; w.z = p.z-V.z;
  s3 = ((-N.x)*w.x + (-N.y)*w.y + (-N.z)*w.z)/(N.x*u.x + N.y*u.y + N.z*u.z);
 
  //4th Plane
  V.x = 0.0; V.y = ny; V.z = 0.0;
  w.x = p.x-V.x; w.y = p.y-V.y; w.z = p.z-V.z;
  //normal of domain plane
  N.x = 0.0; N.y = -1.0; N.z = 0.0;
  s4 = ((-N.x)*w.x + (-N.y)*w.y + (-N.z)*w.z)/(N.x*u.x + N.y*u.y + N.z*u.z);
 

  //Find Smallest s - this means that intersection is the one to use
  s = 100.0;
 
  //+y plane
  if((s1 < s) && (s1 >= 0.0) && (p1.y <= 0.0)){
    s = s1;
  }
  //+z plane
  if((s2 < s) && (s2 >= 0.0)){
    s = s2;
  }
  //-z plane
  if((s3 < s) && (s3 >= 0.0)){
    s = s3;
  }
  //-y plane
  if((s4 < s) && (s4 >= 0.0)){
    s = s4;
  }
    
    
  //Find point of intersection, which is e1
  e2.x = s*u.x + p.x;
  e2.y = s*u.y + p.y;
  e2.z = s*u.z + p.z;
  
  std::cout << "e2.x = " << e2.x << " e2.y = " << e2.y << " e2.z = " << e2.z
	    <<std::endl;

	    }*/
