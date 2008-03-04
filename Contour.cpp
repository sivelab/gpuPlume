#include "Contour.h"
#include <math.h>

Contour::Contour(ParticleControl* partcont, int num){
  pc = partcont;
  
  nx = pc->nx;
  ny = pc->ny;
  nz = pc->nz;
  
  nzdz = pc->nzdz;
  nydy = pc->nydy;
  nxdx = pc->nxdx;
  
  
  //n = number of regions, so there will be n-1 contour values
  n = num;
  num_cValue = n-1;
  cValue = new float[num_cValue];
     
  layer=-1;
  plane_layer_z = -1;
  plane_layer_x = -1;
  plane_layer_y = -1;
    
  contourLayer = new GLuint[4*nzdz];
  glGenBuffers(4*nzdz, contourLayer);
  contourLayer_nx = new GLuint[4*nxdx];
  glGenBuffers(4*nxdx, contourLayer_nx);
  contourLayer_ny = new GLuint[4*nydy];
  glGenBuffers(4*nydy, contourLayer_ny);
  
  //Create tau from pc->tau structure
  tau = new float[nxdx*nzdz*nydy*4];

  for(int k=0;k<nzdz;k++){
    for(int i=0; i < nydy; i++){
      for(int j=0; j < nxdx; j++){
	int idx = k*nydy*nxdx + i*nxdx + j;

	int tidx = k*nydy*nxdx*4 + i*nxdx*4 + j*4;
	
	tau[tidx] = pc->tau[idx].t11;
	tau[tidx+1] = pc->tau[idx].t22;
	tau[tidx+2] = pc->tau[idx].t33;
	tau[tidx+3] = pc->tau[idx].t13;
      }
    }
  }

  //Need to generate contour lines for each axis

  tauValue = 0;

  int width;
  
  glDisable(GL_TEXTURE_2D);
  glEnable(pc->texType);
  glGenTextures(3,tex_id);

  for(plane_normal = 0; plane_normal <=2; plane_normal++){

    if(plane_normal == 0){
      
      max_layer = nz;
      contourValues = new float[num_cValue*nzdz*4];
      numPoints = new int[4*nzdz];
      for(int i=0; i < 4*nzdz; i++)
	numPoints[i] = 0;
 
      width = num_cValue*nzdz;

    }
    else if(plane_normal == 1){
      
      max_layer = nx;
      contourValues = new float[num_cValue*nxdx*4];
      numPoints_nx = new int[4*nxdx];
      for(int i=0; i < 4*nxdx; i++)
	numPoints_nx[i] = 0;
      
      width = num_cValue*nxdx;
      
    }
    else{
      
      max_layer = ny;
      contourValues = new float[num_cValue*nydy*4];
      numPoints_ny = new int[4*nydy];
      for(int i=0; i < 4*nydy; i++)
	numPoints_ny[i] = 0;

      width = num_cValue*nydy;
    }

    setLocalTauValues();
  
    //Get Contour lines for Tau11,Tau22,Tau33,and Tau13
    for(int i=0; i < 4; i++){
      tauValue = i;
      //findContours_Averaging(pc);

      //This Method works a lot better!
      //Need to generalize this method still for each axis aligned plane
      if(plane_normal == 0)
	find_Multiple_Contours();
      else if(plane_normal == 1)
	find_Multiple_Contours_nx();
      else
	find_Multiple_Contours_ny();
    }

    
    glBindTexture(pc->texType, tex_id[plane_normal]);
  
    glTexParameteri(pc->texType, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(pc->texType, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(pc->texType, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(pc->texType, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
          
    glTexImage2D(pc->texType, 0, GL_RGBA32F_ARB, width,1,0, GL_RGBA, GL_FLOAT, contourValues);

    glBindTexture(pc->texType, 0);
    
    delete [] contourValues;
  }
  
  plane_normal = 0;
  tauValue = -1;

  //Set up shader to color contour areas
  contour_shader.addShader("Shaders/contours_vp.glsl",GLSLObject::VERTEX_SHADER);
  contour_shader.addShader("Shaders/contours_fp.glsl",GLSLObject::FRAGMENT_SHADER);
  contour_shader.createProgram();
  uniform_numContours = contour_shader.createUniform("numContours");
  //uniform_tauTex = contour_shader.createUniform("tau");
  uniform_tauValue = contour_shader.createUniform("tauValue");
  uniform_contourTex = contour_shader.createUniform("contourTex");
  uniform_height = contour_shader.createUniform("height");
  uniform_3Dtau = contour_shader.createUniform("tau3d");
    

  glDisable(pc->texType);
  glDisable(GL_TEXTURE_2D);
  glEnable(GL_TEXTURE_3D);

  glGenTextures(1,tex_3D);

  glBindTexture(GL_TEXTURE_3D, tex_3D[0]);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  //glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  //glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
	
  glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA32F_ARB, nxdx, nydy, nzdz,0,GL_RGBA, GL_FLOAT,tau); 
  //CheckErrorsGL("3D texture\n");

  glBindTexture(GL_TEXTURE_3D, 0);

  glDisable(GL_TEXTURE_3D);


  delete [] tau;
  
  /*glDisable(GL_TEXTURE_2D);
  glEnable(pc->texType);
  glGenTextures(3,tex_id);

  glBindTexture(pc->texType, tex_id[0]);
  
  glTexParameteri(pc->texType, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(pc->texType, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(pc->texType, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(pc->texType, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  
   
  //Need to create three texture like this, with nxdx,nydy and nzdz
  glTexImage2D(pc->texType, 0, GL_RGBA32F_ARB, num_cValue*nxdx,1,0, GL_RGBA, GL_FLOAT, contourValues);

  glBindTexture(pc->texType, 0);*/

  glDisable(pc->texType);
   
}

void Contour::setLocalTauValues(){

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

    tauLocalMax[tidx] = pc->tau[idx].t11;
    tauLocalMax[tidx+1] = pc->tau[idx].t22;
    tauLocalMax[tidx+2] = pc->tau[idx].t33;
    tauLocalMax[tidx+3] = pc->tau[idx].t13;

    
    tauLocalMin[tidx] = pc->tau[idx].t11;
    tauLocalMin[tidx+1] = pc->tau[idx].t22;
    tauLocalMin[tidx+2] = pc->tau[idx].t33;
    tauLocalMin[tidx+3] = pc->tau[idx].t13;
    
  }

  //Find max and min tau values for each height value

  for(int k=0; k<nzdz; k++){
    for(int i=0; i<nydy; i++){
      for(int j=0; j<nxdx; j++){

	int idx = k*nydy*nxdx + i*nxdx + j;
	int tidx;
	if(plane_normal == 0){
	  tidx = k*4;
	}
	else if(plane_normal == 1){
	  tidx = j*4;	    
	}
	else{
	  tidx = i*4;
	}
	

	if(pc->tau[idx].t11 > tauLocalMax[tidx])
	  tauLocalMax[tidx] = pc->tau[idx].t11;
	if(pc->tau[idx].t22 > tauLocalMax[tidx+1])
	  tauLocalMax[tidx+1] = pc->tau[idx].t22;
	if(pc->tau[idx].t33 > tauLocalMax[tidx+2])
	  tauLocalMax[tidx+2] = pc->tau[idx].t33;
	if(pc->tau[idx].t13 > tauLocalMax[tidx+3])
	  tauLocalMax[tidx+3] = pc->tau[idx].t13;

	if(pc->tau[idx].t11 < tauLocalMin[tidx])
	  tauLocalMin[tidx] = pc->tau[idx].t11;
	if(pc->tau[idx].t22 < tauLocalMin[tidx+1])
	  tauLocalMin[tidx+1] = pc->tau[idx].t22;
	if(pc->tau[idx].t33 < tauLocalMin[tidx+2])
	  tauLocalMin[tidx+2] = pc->tau[idx].t33;
	if(pc->tau[idx].t13 < tauLocalMin[tidx+3])
	  tauLocalMin[tidx+3] = pc->tau[idx].t13;

      }
    }
  }


}

void Contour::setContourValuesLocally(int k){
  int idx = (k*4)+tauValue;
  int cidx;
  
  //std::cout << "Local Max = " << pc->tauLocalMax[idx] << std::endl;
  //std::cout << "Local Min = " << pc->tauLocalMin[idx] << std::endl;
  
  //Values are stored in array from min to max

  for(int i=0; i < num_cValue; i++){
    cValue[i] = (i+1)*(((tauLocalMax[idx]-tauLocalMin[idx])/(float)n)) + tauLocalMin[idx];
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

void Contour::find_Multiple_Contours(){
  int tauIdx = 0;

  for(int k=0;k<nzdz;k++){
    //Set Contour values
    setContourValuesLocally(k);
    //setContourValuesGlobally(pc,k);
        
    tauIdx = (k*4) + tauValue;
    
    for(int i=0; i < nydy-1; i++){
      for(int j=0; j < nxdx-1; j++){

	int idx = k*nydy*nxdx + i*nxdx + j;
	int idxAbove = k*nydy*nxdx + (i+1)*nxdx + j;
		
	
	
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
	    numPoints[tauIdx]--;
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
void Contour::find_Multiple_Contours_nx(){
  int tauIdx = 0;

  for(int k=0;k<nxdx;k++){
    //Set Contour values
    setContourValuesLocally(k);
    //setContourValuesGlobally(pc,k);
        
    tauIdx = (k*4) + tauValue;
    
    for(int i=0; i < nzdz-1; i++){
      for(int j=0; j < nydy-1; j++){

	int idx = i*nydy*nxdx + j*nxdx + k;

	int idxLeft = i*nydy*nxdx + (j+1)*nxdx + k;

	int idxAbove = (i+1)*nydy*nxdx + j*nxdx + k;

	int idxAboveLeft = (i+1)*nydy*nxdx + (j+1)*nxdx + k;
	
	
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
	  v1 = pc->tau[idxLeft].t11;
	  v2 = pc->tau[idxAbove].t11;
	  v3 = pc->tau[idxAboveLeft].t11;
	  break;
	case 1:
	  v0 = pc->tau[idx].t22;
	  v1 = pc->tau[idxLeft].t22;
	  v2 = pc->tau[idxAbove].t22;
	  v3 = pc->tau[idxAboveLeft].t22;
	  break;
	case 2:
	  v0 = pc->tau[idx].t33;
	  v1 = pc->tau[idxLeft].t33;
	  v2 = pc->tau[idxAbove].t33;
	  v3 = pc->tau[idxAboveLeft].t33;
	  break;
	case 3:
	  v0 = pc->tau[idx].t13;
	  v1 = pc->tau[idxLeft].t13;
	  v2 = pc->tau[idxAbove].t13;
	  v3 = pc->tau[idxAboveLeft].t13;
	  break;
	default:
	  v0 = pc->tau[idx].t11;
	  v1 = pc->tau[idxLeft].t11;
	  v2 = pc->tau[idxAbove].t11;
	  v3 = pc->tau[idxAboveLeft].t11;

	}

	

	for(int c=0; c < num_cValue; c++){

	  localPoints=0;
	  //If isocontour value is between the tau values
	  if(((v0 <= cValue[c]) && (cValue[c] <= v1)) || ((v1 <= cValue[c]) && (cValue[c] <= v0))){
	    //Interpolate to find where value is cValue[c] at point p
	    vec2 p;
	    p.x = p0.x + ((cValue[c]-v0)/(v1-v0))*(p1.x-p0.x);
	    p.y = p0.y + ((cValue[c]-v0)/(v1-v0))*(p1.y-p0.y);
	    numPoints_nx[tauIdx]++;
	    localPoints++;
	    c1List.push_back(p);
	  }
	  if(((v0 <= cValue[c]) && (cValue[c] <= v2)) || ((v2 <= cValue[c]) && (cValue[c] <= v0))){
	    vec2 p;
	    p.x = p0.x + ((cValue[c]-v0)/(v2-v0))*(p2.x-p0.x);
	    p.y = p0.y + ((cValue[c]-v0)/(v2-v0))*(p2.y-p0.y);
	    numPoints_nx[tauIdx]++;
	    localPoints++;
	    c1List.push_back(p);
	  }
	  if(((v1 <= cValue[c]) && (cValue[c] <= v3)) || ((v3 <= cValue[c]) && (cValue[c] <= v1))){
	    vec2 p;
	    p.x = p1.x + ((cValue[c]-v1)/(v3-v1))*(p3.x-p1.x);
	    p.y = p1.y + ((cValue[c]-v1)/(v3-v1))*(p3.y-p1.y);
	    numPoints_nx[tauIdx]++;
	    localPoints++;
	    c1List.push_back(p);
	  }
	  if(((v2 <= cValue[c]) && (cValue[c] <= v3)) || ((v3 <= cValue[c]) && (cValue[c] <= v2))){
	    vec2 p;
	    p.x = p2.x + ((cValue[c]-v2)/(v3-v2))*(p3.x-p2.x);
	    p.y = p2.y + ((cValue[c]-v2)/(v3-v2))*(p3.y-p2.y);
	    numPoints_nx[tauIdx]++;
	    localPoints++;
	    c1List.push_back(p);
	  }
	  if(localPoints%2 != 0){
	    c1List.pop_back();
	    numPoints_nx[tauIdx]--;
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
void Contour::find_Multiple_Contours_ny(){
  int tauIdx = 0;

  for(int k=0;k<nydy;k++){
    //Set Contour values
    setContourValuesLocally(k);
    //setContourValuesGlobally(pc,k);
        
    tauIdx = (k*4) + tauValue;
    
    for(int i=0; i < nzdz-1; i++){
      for(int j=0; j < nxdx-1; j++){

	int idx = i*nydy*nxdx + k*nxdx + j;

	int idxLeft = i*nydy*nxdx + k*nxdx + j+1;

	int idxAbove = (i+1)*nydy*nxdx + k*nxdx + j;

	int idxAboveLeft = (i+1)*nydy*nxdx + k*nxdx + j+1;
	
	
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
	  v1 = pc->tau[idxLeft].t11;
	  v2 = pc->tau[idxAbove].t11;
	  v3 = pc->tau[idxAboveLeft].t11;
	  break;
	case 1:
	  v0 = pc->tau[idx].t22;
	  v1 = pc->tau[idxLeft].t22;
	  v2 = pc->tau[idxAbove].t22;
	  v3 = pc->tau[idxAboveLeft].t22;
	  break;
	case 2:
	  v0 = pc->tau[idx].t33;
	  v1 = pc->tau[idxLeft].t33;
	  v2 = pc->tau[idxAbove].t33;
	  v3 = pc->tau[idxAboveLeft].t33;
	  break;
	case 3:
	  v0 = pc->tau[idx].t13;
	  v1 = pc->tau[idxLeft].t13;
	  v2 = pc->tau[idxAbove].t13;
	  v3 = pc->tau[idxAboveLeft].t13;
	  break;
	default:
	  v0 = pc->tau[idx].t11;
	  v1 = pc->tau[idxLeft].t11;
	  v2 = pc->tau[idxAbove].t11;
	  v3 = pc->tau[idxAboveLeft].t11;

	}

	

	for(int c=0; c < num_cValue; c++){

	  localPoints=0;
	  //If isocontour value is between the tau values
	  if(((v0 <= cValue[c]) && (cValue[c] <= v1)) || ((v1 <= cValue[c]) && (cValue[c] <= v0))){
	    //Interpolate to find where value is cValue[c] at point p
	    vec2 p;
	    p.x = p0.x + ((cValue[c]-v0)/(v1-v0))*(p1.x-p0.x);
	    p.y = p0.y + ((cValue[c]-v0)/(v1-v0))*(p1.y-p0.y);
	    numPoints_ny[tauIdx]++;
	    localPoints++;
	    c1List.push_back(p);
	  }
	  if(((v0 <= cValue[c]) && (cValue[c] <= v2)) || ((v2 <= cValue[c]) && (cValue[c] <= v0))){
	    vec2 p;
	    p.x = p0.x + ((cValue[c]-v0)/(v2-v0))*(p2.x-p0.x);
	    p.y = p0.y + ((cValue[c]-v0)/(v2-v0))*(p2.y-p0.y);
	    numPoints_ny[tauIdx]++;
	    localPoints++;
	    c1List.push_back(p);
	  }
	  if(((v1 <= cValue[c]) && (cValue[c] <= v3)) || ((v3 <= cValue[c]) && (cValue[c] <= v1))){
	    vec2 p;
	    p.x = p1.x + ((cValue[c]-v1)/(v3-v1))*(p3.x-p1.x);
	    p.y = p1.y + ((cValue[c]-v1)/(v3-v1))*(p3.y-p1.y);
	    numPoints_ny[tauIdx]++;
	    localPoints++;
	    c1List.push_back(p);
	  }
	  if(((v2 <= cValue[c]) && (cValue[c] <= v3)) || ((v3 <= cValue[c]) && (cValue[c] <= v2))){
	    vec2 p;
	    p.x = p2.x + ((cValue[c]-v2)/(v3-v2))*(p3.x-p2.x);
	    p.y = p2.y + ((cValue[c]-v2)/(v3-v2))*(p3.y-p2.y);
	    numPoints_ny[tauIdx]++;
	    localPoints++;
	    c1List.push_back(p);
	  }
	  if(localPoints%2 != 0){
	    c1List.pop_back();
	    numPoints_ny[tauIdx]--;
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
  for(int k=0;k<nzdz;k++){
    //Set Contour values
    setContourValuesLocally(k);


    int tauIdx = (k*4) + tauValue;

    for(int i=0; i < nydy; i++){
      for(int j=0; j < nxdx; j++){
	
	int tidx = k*nydy*nxdx*4 + i*nxdx*4 + j*4 + tauValue;
	int tidxAbove = k*nydy*nxdx*4 + (i+1)*nxdx*4 + j*4 + tauValue;
	int tidxBelow = k*nydy*nxdx*4 + (i-1)*nxdx*4 + j*4 + tauValue;
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
	      numPoints[tauIdx]--;
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
  GLfloat* data;
  
  if(plane_normal == 0)
    data = new GLfloat[(numPoints[(k*4)+tauValue])*4];
  else if(plane_normal == 1)
    data = new GLfloat[(numPoints_nx[(k*4)+tauValue])*4];
  else
    data = new GLfloat[(numPoints_ny[(k*4)+tauValue])*4];

  listIter = c1List.begin();
  vec2 p;

  int i=0;
  int idx;
  while(listIter != c1List.end()){
    idx = i*4;
    p = *listIter;
    listIter++;
    i++;

    if(plane_normal == 0){
      data[idx] = p.x;
      data[idx+1] = p.y;
      data[idx+2] = k+0.5;
      data[idx+3] = 1.0;
    }
    else if(plane_normal == 1){
      data[idx] = k-0.5;
      data[idx+1] = p.x;
      data[idx+2] = p.y;
      data[idx+3] = 1.0;

    }
    else{
      data[idx] = p.x;
      data[idx+1] = k-0.5;
      data[idx+2] = p.y;
      data[idx+3] = 1.0;
    }

  }
  if(plane_normal == 0){
    glBindBuffer(GL_ARRAY_BUFFER, contourLayer[(k*4)+tauValue]);
    glBufferData(GL_ARRAY_BUFFER, numPoints[(k*4)+tauValue]*4*sizeof(GLfloat), data, GL_STATIC_DRAW);
  }
  else if(plane_normal == 1){
    glBindBuffer(GL_ARRAY_BUFFER, contourLayer_nx[(k*4)+tauValue]);
    glBufferData(GL_ARRAY_BUFFER, numPoints_nx[(k*4)+tauValue]*4*sizeof(GLfloat), data, GL_STATIC_DRAW);
  }
  else{
    glBindBuffer(GL_ARRAY_BUFFER, contourLayer_ny[(k*4)+tauValue]);
    glBufferData(GL_ARRAY_BUFFER, numPoints_ny[(k*4)+tauValue]*4*sizeof(GLfloat), data, GL_STATIC_DRAW);
  }

  //glBufferData(GL_ARRAY_BUFFER, numPoints[(k*4)+tauValue]*4*sizeof(GLfloat), data, GL_STATIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  delete [] data;
}
void Contour::draw(){

  if(plane_normal == 0){
    max_layer = nz;
    layer = plane_layer_z;
  }
  else if(plane_normal == 1){
    max_layer = nx;
    layer = plane_layer_x;
    //plane_layer = (int)plane_layer_x;
  }
  else{
    max_layer = ny;
    layer = plane_layer_y;
    //plane_layer = (int)plane_layer_y;
  }

  if(layer >= 0 && layer < max_layer && tauValue >= 0){

    int idx = 4*layer + tauValue;
    
    glDisable(GL_COLOR_ARRAY);
    glEnableClientState(GL_VERTEX_ARRAY);
    glLineWidth(1.0);

    glColor4f(1.0,1.0,1.0,1.0);

    if(plane_normal == 0){
      glBindBuffer(GL_ARRAY_BUFFER, contourLayer[idx]);
      glVertexPointer(4,GL_FLOAT,0,0);
      glDrawArrays(GL_LINES,0,numPoints[idx]);
    }
    else if(plane_normal == 1){
      glBindBuffer(GL_ARRAY_BUFFER, contourLayer_nx[idx]);
      glVertexPointer(4,GL_FLOAT,0,0);
      glDrawArrays(GL_LINES,0,numPoints_nx[idx]);
    }
    else{
      glBindBuffer(GL_ARRAY_BUFFER, contourLayer_ny[idx]);
      glVertexPointer(4,GL_FLOAT,0,0);
      glDrawArrays(GL_LINES,0,numPoints_ny[idx]);
    }

    
    //glDrawArrays(GL_LINES,0,numPoints[idx]);

    glLineWidth(1.0);
  }


  

}
void Contour::displayContourLayer(ParticleControl* pc,GLuint texId, int numInRow){

  if(plane_normal == 0){
    max_layer = nz;
    layer = plane_layer_z;
    //plane_layer = (int)plane_layer_z;
  }
  else if(plane_normal == 1){
    max_layer = nx;
    //plane_layer = (int)plane_layer_x;
    layer = plane_layer_x;
  }
  else{
    max_layer = ny;
    layer = plane_layer_y;
    //plane_layer = (int)plane_layer_y;
  }


  if(layer >= 0 && layer < max_layer && tauValue >= 0){

    //glEnable(GL_COLOR_MATERIAL);
    //glEnable(GL_BLEND);
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);     
    
    //glEnable(pc->texType);

    //int s = 0;
    //int t = 0;

    //s = (int)(layer % numInRow) * nxdx;
    //t = (int)(floor(layer/(float)numInRow) * nydy);   
    float s,t,r;
	
    if(plane_normal == 0){
      s = 0.0;
      t = 0.0;
      r = layer/((float)nzdz - 1.0);
    }
    //Vertical Texture coordinates
    else if(plane_normal == 1){
      s = layer/((float)nxdx - 1.0);
      t = 0.0;
      r = 0.0;
    }
    else{
      s = 0.0;
      t = layer/((float)nydy - 1.0);
      r = 0.0;
    }

    contour_shader.activate();

    //the number of contour values
    glUniform1iARB(uniform_numContours, num_cValue);

    //Selects which tau value to display in contours
    //glUniform1iARB(uniform_tauValue, tauValue);

    if(tauValue == 0){
      glUniform4fARB(uniform_tauValue,1.0f,0.0f,0.0f,0.0f);
    }
    else if(tauValue == 1){
      glUniform4fARB(uniform_tauValue,0.0f,1.0f,0.0f,0.0f);
    }
    else if(tauValue == 2){
      glUniform4fARB(uniform_tauValue,0.0f,0.0f,1.0f,0.0f);
    }
    else{
      glUniform4fARB(uniform_tauValue,0.0f,0.0f,0.0f,1.0f);
    }
    //Selects the layer to display
    glUniform1iARB(uniform_height, layer);


    glEnable(pc->texType);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(pc->texType,tex_id[plane_normal]);
    glUniform1iARB(uniform_contourTex,1);
    glDisable(pc->texType);

    glEnable(GL_TEXTURE_3D);
    glActiveTexture(GL_TEXTURE0);
    //glBindTexture(pc->texType, texId);
    glBindTexture(GL_TEXTURE_3D,tex_3D[0]);
    //glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    //glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glUniform1iARB(uniform_3Dtau,0);
    //glUniform1iARB(uniform_tauTex, 0); 
 
    /*glBegin(GL_QUADS);
    {
      glNormal3f(0.0, 0.0, 1.0);
      glTexCoord2f(s,t);              glVertex3f(0, 0, layer);
      glTexCoord2f(s+nxdx, t);          glVertex3f(nx, 0, layer);
      glTexCoord2f(s+nxdx, t+nydy);       glVertex3f(nx, ny, layer);
      glTexCoord2f(s, t+nydy);          glVertex3f(0, ny, layer);
    }
    glEnd();*/

    if(plane_normal == 0){
      glBegin(GL_QUADS);
      {
	glNormal3f(0.0, 0.0, 1.0);
	glTexCoord3f( s, t, r);         glVertex3f(0, 0, layer);
	glTexCoord3f( s+1, t, r);       glVertex3f(nx, 0,layer);
	glTexCoord3f( s+1, t+1, r);     glVertex3f(nx, ny,layer);
	glTexCoord3f( s,  t+1, r);      glVertex3f(0, ny,layer);
      }
      glEnd();
    }
    else if(plane_normal == 1){
      glBegin(GL_QUADS);
      glNormal3f(1.0,0.0,0.0);
      glTexCoord3f( s, t, r);         glVertex3f(layer, 0, 0);
      glTexCoord3f( s, t, r+1);       glVertex3f(layer, 0, nz);
      glTexCoord3f( s, t+1, r+1);     glVertex3f(layer, ny, nz);
      glTexCoord3f( s, t+1, r);       glVertex3f(layer, ny, 0);

      glEnd();
    }
    else{
      glBegin(GL_QUADS);
      glNormal3f(0.0,1.0,0.0);
      glTexCoord3f( s, t, r);         glVertex3f(0, layer, 0);
      glTexCoord3f( s+1, t, r);       glVertex3f(nx, layer, 0);
      glTexCoord3f( s+1, t, r+1);     glVertex3f(nx, layer, nz);
      glTexCoord3f( s, t, r+1);       glVertex3f(0, layer, nz);

      glEnd();
    }


    contour_shader.deactivate();

    
    //??Do I need to do this??
    //glTexParameteri(pc->texType, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    //glTexParameteri(pc->texType, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    glDisable(GL_TEXTURE_3D);
    //glDisable(pc->texType);
    //glDisable(GL_BLEND);
    //glDisable(GL_COLOR_MATERIAL);
  
    

  }
}

void Contour::decreaseContourLayer(){
  //layer--;
  //if(layer < -1) layer = -1;
  if(plane_normal == 0){
    plane_layer_z--;
    if(plane_layer_z < -1)
      plane_layer_z = -1;
  }
  else if(plane_normal == 1){
    plane_layer_x--;
    if(plane_layer_x < -1)
      plane_layer_x = -1;
  }
  else{
    plane_layer_y--;
    if(plane_layer_y < -1)
      plane_layer_y = -1;
  }
  


}
void Contour::increaseContourLayer(){
  //layer++;
  //if(layer > max_layer) layer = max_layer;

  if(plane_normal == 0){
    plane_layer_z++;
    if(plane_layer_z > max_layer)
      plane_layer_z = max_layer;
  }
  else if(plane_normal == 1){
    plane_layer_x++;
    if(plane_layer_x > max_layer)
      plane_layer_x = max_layer;
  }
  else{
    plane_layer_y++;
    if(plane_layer_y > max_layer)
      plane_layer_y = max_layer;
  }

}
void Contour::switchPlane(){
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
   

}
