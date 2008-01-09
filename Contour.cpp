#include "Contour.h"
#include <math.h>

Contour::Contour(ParticleControl* pc){
  nx = pc->nx;
  ny = pc->ny;
  nz = pc->nz;

  //n = number of regions, so there will be n-1 contour values
  n = 5;
  num_cValue = n-1;
  cValue = new float[num_cValue];
  
  contourValues = new float[num_cValue*nz*4];
  
  layer=-1;
  tauValue = 0;

  numPoints = new int[4*nz];
  for(int i=0; i < 4*nz; i++)
    numPoints[i] = 0;

  contourLayer = new GLuint[4*nz];
  glGenBuffers(4*nz, contourLayer);

  
  //Create tau from pc->tau structure
  tau = new float[nx*nz*ny*4];

  for(int k=0;k<nz;k++){
    for(int i=0; i < ny; i++){
      for(int j=0; j < nx; j++){
	int idx = k*ny*nx + i*nx + j;

	int tidx = k*ny*nx*4 + i*nx*4 + j*4;
	
	tau[tidx] = pc->tau[idx].t11;
	tau[tidx+1] = pc->tau[idx].t22;
	tau[tidx+2] = pc->tau[idx].t33;
	tau[tidx+3] = pc->tau[idx].t13;
      }
    }
  }

  //Get Contour lines for Tau11,Tau22,Tau33,and Tau13
  for(int i=0; i < 4; i++){
    tauValue = i;
    //findContours_Averaging(pc);
    find_Multiple_Contours(pc);
  }

  //Delete tau
  delete [] tau;

  

  //Prints out values of array contourValues
  /*
  for(int i=0;i<4; i++){
    tauValue = i;

    for(int k=0; k<nz; k++){
      int cidx = (k*4*num_cValue) + tauValue*4;

      //Find cValues based on layer value
      for(int i=0; i < num_cValue; i++){
	cValue[i] = contourValues[i+cidx];
	//std::cout << cValue[i] << std::endl;
	//std::cout << "C value " << i << " equals " << cValue[i] << std::endl;
      }

    }
  }
  */

  tauValue = -1;

  //Set up shader to color contour areas
  contour_shader.addShader("Shaders/contours_vp.glsl",GLSLObject::VERTEX_SHADER);
  contour_shader.addShader("Shaders/contours_fp.glsl",GLSLObject::FRAGMENT_SHADER);
  contour_shader.createProgram();
  uniform_numContours = contour_shader.createUniform("numContours");
  uniform_tauTex = contour_shader.createUniform("tau");
  uniform_tauValue = contour_shader.createUniform("tauValue");
  uniform_contourTex = contour_shader.createUniform("contourTex");
  uniform_height = contour_shader.createUniform("height");
  
  
  glDisable(GL_TEXTURE_2D);
  glEnable(pc->texType);
  glGenTextures(1,tex_id);

  glBindTexture(pc->texType, tex_id[0]);
  
  glTexParameteri(pc->texType, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(pc->texType, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(pc->texType, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(pc->texType, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  
  /*
  int cidx = 0;
  GLfloat *data = new GLfloat[num_cValue*nz*4];
  for(int k=0; k < nz; k++){
    for(int i=0; i < num_cValue; i++){
      for(int j=0; j < 4; j++){
	int idx = k*4*num_cValue + i*4 + j;
	data[idx] = contourValues[cidx];
	cidx++;
	//std::cout << data[idx] << std::endl;
      }
    }
    }*/

  
  glTexImage2D(pc->texType, 0, GL_RGBA32F_ARB, num_cValue*nz*4,1,0, GL_RGBA, GL_FLOAT, contourValues);

  glBindTexture(pc->texType, 0);

  //delete [] data;
  glDisable(pc->texType);
  
}


void Contour::setContourValuesLocally(ParticleControl* pc,int k){
  int idx = (k*4)+tauValue;
  int cidx;
  
  //std::cout << "Local Max = " << pc->tauLocalMax[idx] << std::endl;
  //std::cout << "Local Min = " << pc->tauLocalMin[idx] << std::endl;
  
  //Values are stored in array from min to max

  for(int i=0; i < num_cValue; i++){
    cValue[i] = (i+1)*(((pc->tauLocalMax[idx]-pc->tauLocalMin[idx])/(float)n)) + pc->tauLocalMin[idx];
    //std::cout << "C value " << i << " equals " << cValue[i] << std::endl;
    
    cidx = k*4*num_cValue + 4*i + tauValue;
    contourValues[cidx] = cValue[i];

  }

}
void Contour::setContourValuesGlobally(ParticleControl* pc,int k){
  
  int cidx;

  for(int i=0; i < n-1; i++){
    cValue[i] = (i+1)*(((pc->tauMax[tauValue]-pc->tauMin[tauValue])/(float)n)) + pc->tauMin[tauValue];
    //std::cout << "C value " << i << " equals " << cValue[i] << std::endl;
    
    cidx = k*4*num_cValue + i*4 + tauValue;
    contourValues[cidx] = cValue[i];

  }

}

void Contour::find_Multiple_Contours(ParticleControl* pc){
  int tauIdx = 0;

  for(int k=0;k<nz;k++){
    //Set Contour values
    setContourValuesLocally(pc,k);
    //setContourValuesGlobally(pc,k);
    
    tauIdx = (k*4) + tauValue;

    for(int i=0; i < ny-1; i++){
      for(int j=0; j < nx-1; j++){
	int idx = k*ny*nx + i*nx + j;
	int idxAbove = k*ny*nx + (i+1)*nx + j;
		
	
	
	//These are the verices of the cell.
	//Lower left corner
	p0.x = (float)j+0.5f;
	p0.y = (float)i+0.5f;
	//Lower right corner
	p1.x = (float)j+1.5f;
	p1.y = (float)i+0.5f;
	//Upper left corner
	p2.x = (float)j+0.5f;
	p2.y = (float)i+1.5f;
	//Upper right corner
	p3.x = (float)j+1.5f;
	p3.y = (float)i+1.5f;


	switch(tauValue){
	  
	case 0:
	  v0 = pc->tau[idx].t11;
	  v1 = pc->tau[idx+1].t11;
	  v2 = pc->tau[idxAbove].t11;
	  v3 = pc->tau[idxAbove+1].t11;
	  break;
	case 1:
	  v0 = pc->tau[idx].t22;
	  v1 = pc->tau[idx+1].t22;
	  v2 = pc->tau[idxAbove].t22;
	  v3 = pc->tau[idxAbove+1].t22;
	  break;
	case 2:
	  v0 = pc->tau[idx].t33;
	  v1 = pc->tau[idx+1].t33;
	  v2 = pc->tau[idxAbove].t33;
	  v3 = pc->tau[idxAbove+1].t33;
	  break;
	case 3:
	  v0 = pc->tau[idx].t13;
	  v1 = pc->tau[idx+1].t13;
	  v2 = pc->tau[idxAbove].t13;
	  v3 = pc->tau[idxAbove+1].t13;
	  break;
	default:
	  v0 = pc->tau[idx].t11;
	  v1 = pc->tau[idx+1].t11;
	  v2 = pc->tau[idxAbove].t11;
	  v3 = pc->tau[idxAbove+1].t11;

	}

	

	for(int c=0; c < num_cValue; c++){

	  localPoints=0;
	  //If isocontour value is between the tau values
	  if(((v0 <= cValue[c]) && (cValue[c] <= v1)) || ((v1 <= cValue[c]) && (cValue[c] <= v0))){
	    //Interpolate to find where value is cValue[c] at point p
	    vec2 p;
	    p.x = p0.x + ((cValue[c]-v0)/(v1-v0))*(p1.x-p0.x);
	    p.y = p0.y + ((cValue[c]-v0)/(v1-v0))*(p1.y-p0.y);
	    numPoints[tauIdx]++;
	    localPoints++;
	    c1List.push_back(p);
	  }
	  if(((v0 <= cValue[c]) && (cValue[c] <= v2)) || ((v2 <= cValue[c]) && (cValue[c] <= v0))){
	    vec2 p;
	    p.x = p0.x + ((cValue[c]-v0)/(v2-v0))*(p2.x-p0.x);
	    p.y = p0.y + ((cValue[c]-v0)/(v2-v0))*(p2.y-p0.y);
	    numPoints[tauIdx]++;
	    localPoints++;
	    c1List.push_back(p);
	  }
	  if(((v1 <= cValue[c]) && (cValue[c] <= v3)) || ((v3 <= cValue[c]) && (cValue[c] <= v1))){
	    vec2 p;
	    p.x = p1.x + ((cValue[c]-v1)/(v3-v1))*(p3.x-p1.x);
	    p.y = p1.y + ((cValue[c]-v1)/(v3-v1))*(p3.y-p1.y);
	    numPoints[tauIdx]++;
	    localPoints++;
	    c1List.push_back(p);
	  }
	  if(((v2 <= cValue[c]) && (cValue[c] <= v3)) || ((v3 <= cValue[c]) && (cValue[c] <= v2))){
	    vec2 p;
	    p.x = p2.x + ((cValue[c]-v2)/(v3-v2))*(p3.x-p2.x);
	    p.y = p2.y + ((cValue[c]-v2)/(v3-v2))*(p3.y-p2.y);
	    numPoints[tauIdx]++;
	    localPoints++;
	    c1List.push_back(p);
	  }
	  if(localPoints%2 != 0){
	    c1List.pop_back();
	    numPoints[tauIdx] -= localPoints;
	    std::cout << "There are " << localPoints << " points in a cell" << std::endl;
	  }
	}

      }
    }
    //Store List as vbo; Clear List;
    putListinVBO(k); 
    c1List.clear();

    //std::cout << "Num of Points : " << numPoints[tauIdx] << " for layer " << k << std::endl;
  }
  
}


void Contour::findContours_Averaging(ParticleControl* pc){
  float div;

  //will need to change it to k < nz
  for(int k=0;k<nz;k++){
    //Set Contour values
    setContourValuesLocally(pc,k);


    int tauIdx = (k*4) + tauValue;

    for(int i=0; i < ny; i++){
      for(int j=0; j < nx; j++){
	
	int tidx = k*ny*nx*4 + i*nx*4 + j*4 + tauValue;
	int tidxAbove = k*ny*nx*4 + (i+1)*nx*4 + j*4 + tauValue;
	int tidxBelow = k*ny*nx*4 + (i-1)*nx*4 + j*4 + tauValue;
	int tidxLeft = tidx - 4;
	int tidxRight = tidx + 4;

	int tidxBelowLeft = tidxBelow - 4;
	int tidxBelowRight = tidxBelow + 4;
	int tidxAboveLeft = tidxAbove - 4;
	int tidxAboveRight = tidxAbove + 4;

	if(tau[tidx] != 0.0){


	  //These are the verices of the cell.
	  //Lower left corner
	  p0.x = j;
	  p0.y = i;
	  //Lower right corner
	  p1.x = j+1;
	  p1.y = i;
	  //Upper left corner
	  p2.x = j;
	  p2.y = i+1;
	  //Upper right corner
	  p3.x = j+1;
	  p3.y = i+1;

	  //Determine value at point p0 (lower left corner of cell)
	  if(p0.x==0 && p0.y==0)
	    v0 = tau[tidx];
	  else if(p0.x==0 && p0.y != 0){
	    if(tau[tidxBelow] == 0.0)
	      v0 = tau[tidx];
	    else
	      v0 = (tau[tidx] + tau[tidxBelow])/2.0;
	  }
	  else if(p0.x!=0 && p0.y==0){
	    if(tau[tidxLeft] == 0.0)
	      v0 = tau[tidx];
	    else
	      v0 = (tau[tidx] + tau[tidxLeft])/2.0;
	  }
	  else{
	    div = 4.0;
	    if(tau[tidxLeft]== 0.0)
	      div -= 1.0;
	    if(tau[tidxBelow] == 0)
	      div -= 1.0;
	    if(tau[tidxBelowLeft] == 0)
	      div -= 1.0;

	    v0 = (tau[tidx] + tau[tidxLeft]
		  + tau[tidxBelow] + tau[tidxBelowLeft])/div;
	  }
	  //Determine value at point p1 (lower right corner of cell)
	  if(p1.x == nx && p1.y==0)
	    v1 = tau[tidx];
	  else if(p1.x != nx && p1.y==0){
	    if(tau[tidxRight] == 0)
	      v1 = tau[tidx];
	    else
	      v1 = (tau[tidx] + tau[tidxRight])/2.0;
	  }
	  else if(p1.x == nx && p1.y!=0){
	    if(tau[tidxBelow] == 0)
	      v1 = tau[tidx];
	    else
	      v1 = (tau[tidx] + tau[tidxBelow])/2.0;
	  }
	  else{
	    div = 4.0;
	    if(tau[tidxRight] == 0)
	      div -= 1.0;
	    if(tau[tidxBelow] == 0)
	      div -= 1.0;
	    if(tau[tidxBelowRight] == 0)
	      div -= 1.0;
	    
	    v1 = (tau[tidx] + tau[tidxRight]
		  + tau[tidxBelow] + tau[tidxBelowRight])/div;
	  }
	  //Determine value at point p2 (upper left corner of cell)
	  if(p2.x == 0 && p2.y == ny)
	    v2 = tau[tidx];
	  else if(p2.x == 0 && p2.y != ny){
	    if(tau[tidxAbove] == 0)
	      v2 = tau[tidx];
	    else
	      v2 = (tau[tidx] + tau[tidxAbove])/2.0;
	  }
	  else if(p2.x != 0 && p2.y == ny){
	    if(tau[tidxLeft]== 0)
	      v2 = tau[tidx];
	    else
	      v2 = (tau[tidx] + tau[tidxLeft])/2.0;
	  }
	  else{
	    div = 4.0;
	    if(tau[tidxLeft]== 0)
	      div -= 1.0;
	    if(tau[tidxAbove] == 0)
	      div -= 1.0;
	    if(tau[tidxAboveLeft] == 0)
	      div -= 1.0;

	    v2 = (tau[tidx] + tau[tidxLeft]
		  + tau[tidxAbove] + tau[tidxAboveLeft])/div;
	  }
	  //Determine value at point p3 (upper right corner of cell)
	  if(p3.x == nx && p3.y == ny)
	    v3 = tau[tidx];
	  else if(p3.x ==nx && p3.y != ny){
	    if(tau[tidxAbove] == 0)
	      v3 = tau[tidx];
	    else
	      v3 = (tau[tidx] + tau[tidxAbove])/2.0;
	  }
	  else if(p3.x != nx && p3.y == ny){
	    if(tau[tidxRight] == 0)
	      v3 = tau[tidx];
	    else
	      v3 = (tau[tidx] + tau[tidxRight])/2.0;
	  }
	  else{
	    div = 4.0;
	    if(tau[tidxRight] == 0)
	      div -= 1.0;
	    if(tau[tidxAbove] == 0)
	      div -= 1.0;
	    if(tau[tidxAboveRight] == 0)
	      div -= 1.0;

	    v3 = (tau[tidx] + tau[tidxRight] 
		  + tau[tidxAbove] + tau[tidxAboveRight])/div;
	  }

	  //////////////////////////////////////////////
	  //Now check to see where the contour lines are
	  //////////////////////////////////////////////

	  //vec2 p;
	
	  for(int c=0; c < num_cValue; c++){

	    localPoints=0;
	    //If isocontour value is between the tau values
	    if(((v0 <= cValue[c]) && (cValue[c] <= v1)) || ((v1 <= cValue[c]) && (cValue[c] <= v0))){
	      //Interpolate to find where value is cValue[c] at point p
	      vec2 p;
	      p.x = p0.x + ((cValue[c]-v0)/(v1-v0))*(p1.x-p0.x);
	      p.y = p0.y + ((cValue[c]-v0)/(v1-v0))*(p1.y-p0.y);
	      numPoints[tauIdx]++;
	      localPoints++;
	      c1List.push_back(p);
	    }
	    if(((v0 <= cValue[c]) && (cValue[c] <= v2)) || ((v2 <= cValue[c]) && (cValue[c] <= v0))){
	      vec2 p;
	      p.x = p0.x + ((cValue[c]-v0)/(v2-v0))*(p2.x-p0.x);
	      p.y = p0.y + ((cValue[c]-v0)/(v2-v0))*(p2.y-p0.y);
	      numPoints[tauIdx]++;
	      localPoints++;
	      c1List.push_back(p);
	    }
	    if(((v1 <= cValue[c]) && (cValue[c] <= v3)) || ((v3 <= cValue[c]) && (cValue[c] <= v1))){
	      vec2 p;
	      p.x = p1.x + ((cValue[c]-v1)/(v3-v1))*(p3.x-p1.x);
	      p.y = p1.y + ((cValue[c]-v1)/(v3-v1))*(p3.y-p1.y);
	      numPoints[tauIdx]++;
	      localPoints++;
	      c1List.push_back(p);
	    }
	    if(((v2 <= cValue[c]) && (cValue[c] <= v3)) || ((v3 <= cValue[c]) && (cValue[c] <= v2))){
	      vec2 p;
	      p.x = p2.x + ((cValue[c]-v2)/(v3-v2))*(p3.x-p2.x);
	      p.y = p2.y + ((cValue[c]-v2)/(v3-v2))*(p3.y-p2.y);
	      numPoints[tauIdx]++;
	      localPoints++;
	      c1List.push_back(p);
	    }
	    if(localPoints%2 != 0){
	      c1List.pop_back();
	      numPoints[tauIdx] -= localPoints;
	      std::cout << "There are " << localPoints << " points in a cell" << std::endl;
	    }
	  }

	}

      }
    }
    //Store List as vbo; Clear List;
    putListinVBO(k); 
    c1List.clear();

    //std::cout << "Num of Points : " << numPoints[tauIdx] << " for layer " << k << std::endl; 
  }
  

}
void Contour::putListinVBO(int k){
  GLfloat* data = new GLfloat[(numPoints[(k*4)+tauValue])*4];
 

  listIter = c1List.begin();
  vec2 p;

  int i=0;
  int idx;
  while(listIter != c1List.end()){
    idx = i*4;
    p = *listIter;
    listIter++;
    i++;

    data[idx] = p.x;
    data[idx+1] = p.y;
    data[idx+2] = k+0.1;
    data[idx+3] = 1.0;

  }
  
  glBindBuffer(GL_ARRAY_BUFFER, contourLayer[(k*4)+tauValue]);
  glBufferData(GL_ARRAY_BUFFER, numPoints[(k*4)+tauValue]*4*sizeof(GLfloat), data, GL_STATIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  delete [] data;
}
void Contour::draw(){
  if(layer >= 0 && layer < nz && tauValue >= 0){

    int idx = 4*layer + tauValue;
    
    glDisable(GL_COLOR_ARRAY);
    glEnableClientState(GL_VERTEX_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, contourLayer[idx]);
    glVertexPointer(4,GL_FLOAT,0,0);

    glColor4f(1.0,1.0,1.0,1.0);
    glDrawArrays(GL_LINES,0,numPoints[idx]);

  }


}
void Contour::displayContourLayer(ParticleControl* pc,GLuint texId, int numInRow){

  if(layer >= 0 && layer < nz && tauValue >= 0){

    glPushMatrix();  

    //glEnable(GL_COLOR_MATERIAL);
    //glEnable(GL_BLEND);
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);     
    
    glEnable(pc->texType);
    
    glTexParameteri(pc->texType, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(pc->texType, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    int s = 0;
    int t = 0;

    s = (int)(layer % numInRow) * nx;
    t = (int)(floor(layer/(float)numInRow) * ny);   

    contour_shader.activate();
  
    glUniform1iARB(uniform_numContours, num_cValue);
    glUniform1iARB(uniform_tauValue, tauValue);
    glUniform1iARB(uniform_height, layer);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(pc->texType,tex_id[0]);
    glUniform1iARB(uniform_contourTex,1);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(pc->texType, texId);
    glUniform1iARB(uniform_tauTex, 0); 

    glBegin(GL_QUADS);
    {
      glNormal3f(0.0, 0.0, 1.0);
      glTexCoord2f(s,t);              glVertex3f(0, 0, layer);
      glTexCoord2f(s+nx, t);          glVertex3f(nx, 0, layer);
      glTexCoord2f(s+nx, t+ny);       glVertex3f(nx, ny, layer);
      glTexCoord2f(s, t+ny);          glVertex3f(0, ny, layer);
    }
    glEnd();
    contour_shader.deactivate();



    glTexParameteri(pc->texType, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(pc->texType, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    //glBindTexture(pc->texType,0);

    glDisable(pc->texType);
    //glDisable(GL_BLEND);
    //glDisable(GL_COLOR_MATERIAL);

    glPopMatrix();
  
  }
}

void Contour::decreaseContourLayer(){
  layer--;
  if(layer < -1) layer = -1;

}
void Contour::increaseContourLayer(){
  layer++;
  if(layer > nz) layer = nz;
}
