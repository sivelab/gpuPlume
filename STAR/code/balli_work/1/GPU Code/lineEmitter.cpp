#include <iostream>
#include "lineEmitter.h"
#include "Random.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

LineEmitter::LineEmitter(float x, float y, float z,
			 float xe, float ye, float ze,
			 float rate,
			 int w,int h,std::list<int>* ind, GLSLObject* emit_shader,
			 std::vector<float>* randoms, wind* sig,
			 int dx, int dy, int dz)
{
  xpos = x;
  ypos = y;
  zpos = z;

  xpos_end = xe;
  ypos_end = ye;
  zpos_end = ze;

  nxdx = dx;
  nydy = dy;
  nzdz = dz;

  reuse = false;
  lifeTime = -1.0;

  releaseRate = rate;

  numToEmit = 0;

  emitTime = 0;
  remTime = 0;
  remAmount = 0;

  twidth = w;
  theight = h;

  indices = ind;
 
  shader = emit_shader;

  random_values = randoms;
  sigma = sig;
  
  //counter used to step through the random values
  curr = 0;

}
LineEmitter::~LineEmitter(){}

void LineEmitter::getReleasedPosition(float*x,float*y,float*z){
  *x = xpos + offsetx;
  *y = ypos + offsety;
  *z = zpos + offsetz;
}

void LineEmitter::Draw(){
  glDisable(GL_TEXTURE_RECTANGLE_ARB);

  glEnable(GL_COLOR_MATERIAL);

  GLint lwidth;
  glGetIntegerv(GL_LINE_WIDTH, &lwidth);
  glLineWidth(4.0);
  
  glPushMatrix();
  glColor4f(0.0, 0.0, 1.0, 1.0);
  glBegin(GL_LINES);
  glVertex3f(xpos, ypos, zpos);
  glVertex3f(xpos_end, ypos_end, zpos_end);
  glEnd();
  glPopMatrix();

  glLineWidth(lwidth);

  glDisable(GL_BLEND);
  glDisable(GL_COLOR_MATERIAL);

  glEnable(GL_TEXTURE_RECTANGLE_ARB);
}

void LineEmitter::setVertices(){
  posCoord.clear();

  // to generate a random point along a line, we need only query a
  // single random variable between 0 and 1 that can be used as a
  // parameter for a parameterized line function.

  float t;
  for(int i=0; i < numToEmit; i++){

    t = Random::uniform();

    offsetx = xpos + t * (xpos_end - xpos);
    offsety = ypos + t * (ypos_end - ypos);
    offsetz = zpos + t * (zpos_end - zpos);
    
    posCoord.push_back(offsetx);
    posCoord.push_back(offsety);
    posCoord.push_back(offsetz);
    posCoord.push_back(lifeTime);

  }
}
int LineEmitter::EmitParticle(bool odd,GLuint pos0, GLuint pos1,
			      float time_step,GLuint prime0,GLuint prime1){
  int p_index;
  float px,py,pz;

  setEmitAmount(time_step);

  int count = 0;
  if(!indices->empty()){
    //THIS Method *seems* to work now!
    //Punch Hole method. Need to set drawbuffer and activate shader.
    if(Punch_Hole){

      if(odd){
	GLenum buffers[] = {GL_COLOR_ATTACHMENT0_EXT,GL_COLOR_ATTACHMENT2_EXT};
	glDrawBuffers(2,buffers);
      }
      else{ 
	GLenum buffers[] = {GL_COLOR_ATTACHMENT1_EXT,GL_COLOR_ATTACHMENT3_EXT};
	glDrawBuffers(2,buffers);
      }

      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      gluOrtho2D(0, twidth, 0, theight);
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();     

    }
    	
    //Do this for each particle that is being emitted.
    for(int i = 0; i < numToEmit; i++){
      if(!indices->empty()){
	count++;
	//First get available index
	p_index = indices->back();
	indices->pop_back();
	if(reuse){
	  pIndex newIndex;
	  newIndex.id = p_index;
	  newIndex.time = 0;
	  indicesInUse->push_back(newIndex);	   
	}	 

	float t = Random::uniform();
	offsetx = xpos + t * (xpos_end - xpos);
	offsety = ypos + t * (ypos_end - ypos);
	offsetz = zpos + t * (zpos_end - zpos);
    
	//Determine the coordinates into the position texture
	s = (p_index%twidth);
	t = (p_index/twidth);	  

	//Determine initial prime value
	int p2idx = ((int)offsetz)*nydy*nxdx + ((int)offsety)*nxdx + (int)offsetx;
	px = random_values->at(curr)*(sigma[p2idx].u);
	py = random_values->at(curr+1)*(sigma[p2idx].v);
	pz = random_values->at(curr+2)*(sigma[p2idx].w);
	//std::cout << random_values->at(curr) << std::endl;
	curr += 3;
	if(curr >= random_values->size())
	  curr = 0;

	if(Punch_Hole){
	  glPointSize(1.0);

	  shader->activate();

	  glViewport(s,t,1,1);
	  //std::cout << "particle num= " << p_index << "  s = " << s << "  t = " << t << std::endl;
	  glBegin(GL_POINTS);
	  {
	    //passes initial prime into shader
	    glNormal3f(px,py,pz);
	    //passes initial position into shader
	    glColor4f(offsetx, offsety, offsetz, lifeTime);
	    glVertex2f(0.5, 0.5);
	    //glVertex2f(s,t);
	  }
	  glEnd();
	  shader->deactivate();

	}
	else{
	  // Second mechanism to release particles.  Uses a texture
	  // copy which may likely be more expensive even though we're
	  // doing 1x1 pixels (but many times).  This operation
	  // appears to work consistently.
	  GLfloat value[4];
	  value[0] = offsetx;
	  value[1] = offsety;
	  value[2] = offsetz;
	  value[3] = lifeTime;

	  // write there via a glTexSubImage2D
	  if(odd)
	    // will read from this texture next
	    glBindTexture(GL_TEXTURE_RECTANGLE_ARB, pos0);
	  else 
	    // will read from this texture next
	    glBindTexture(GL_TEXTURE_RECTANGLE_ARB, pos1);

	  glTexSubImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, s, t, 1, 1, GL_RGBA, GL_FLOAT, value);
	}
  
      }
    }

    if(Punch_Hole){
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      gluOrtho2D(-1, 1, -1, 1);
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
    }
    
  }
  return count;
   
}
