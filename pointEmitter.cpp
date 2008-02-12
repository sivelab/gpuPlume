#include "pointEmitter.h"

PointEmitter::PointEmitter(float x,float y,float z,float rate, int w, 
			   int h,std::list<int>* ind,GLSLObject* emit_shader,
			   std::vector<float>* randoms,wind* sig,
			   int dx,int dy,int dz){

  xpos = x;
  ypos = y;
  zpos = z;

  nxdx = dx;
  nydy = dy;
  nzdz = dz;

  reuse = false;
  lifeTime = -1.0;

  releaseRate = rate;
  numToEmit = 0;

  emitTime = 0;
  remTime = 0;

  twidth = w;
  theight = h;

  indices = ind; 

  shader = emit_shader;

  random_values = randoms;
  sigma = sig;
  
  //counter used to step through the random values
  curr = 0;

}
PointEmitter::~PointEmitter(){}

void PointEmitter::setVertices(){

  posCoord.clear();

  for(int i=0; i<numToEmit; i++){
    posCoord.push_back(xpos);
    posCoord.push_back(ypos);
    posCoord.push_back(zpos);
    posCoord.push_back(lifeTime);
    
  }
}
int PointEmitter::EmitParticle(bool odd,GLuint pos0, GLuint pos1,
			       float time_step,GLuint prime0,GLuint prime1){
  int p_index;

  float px,py,pz;

  switch(releaseType){
  case perSecond:
    if(!timeToEmit(time_step))
	numToEmit = 0;
    break;
  case onePerKeyPress:
    setNumToEmit(1);
    emit = false;
    break;
  case instantaneous:
    setNumToEmit(twidth*theight);
    emit = false;
    break;
  default:
    //This will release per time step with numToEmit already defined!!!
    break;
  }

  if(!indices->empty()){
    //THIS Method *seems* to work now!
    //Punch Hole method. Need to set drawbuffer and activate shader.
    if(Punch_Hole){

      /*if(odd)
	glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
      else 
      glDrawBuffer(GL_COLOR_ATTACHMENT1_EXT);*/
      
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
	
	//First get available index
	p_index = indices->back();
	indices->pop_back();
	if(reuse){
	  pIndex newIndex;
	  newIndex.id = p_index;
	  newIndex.time = 0;
	  indicesInUse->push_back(newIndex);	   
	}	 

	//Determine the coordinates into the position and prime texture
	s = (p_index%twidth);
	t = (p_index/twidth);	  

	//Determine initial prime value
	int p2idx = ((int)zpos)*nydy*nxdx + ((int)ypos)*nxdx + (int)xpos;
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
	    //passes particle position into shader
	    glColor4f(xpos, ypos, zpos, lifeTime);
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
	  value[0] = xpos;
	  value[1] = ypos;
	  value[2] = zpos;
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
  return numToEmit;
   
}
