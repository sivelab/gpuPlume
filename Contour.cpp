#include "Contour.h"

Contour::Contour(ParticleControl* pc){
  nx = pc->nx;
  ny = pc->ny;
  nz = pc->nz;

  //n = number of regions, so there will be n-1 contour values
  int n = 10;
  num_cValue = n-1;
  cValue = new float[n-1];
  for(int i=0; i < n-1; i++){
    cValue[i] = (i+1)*((pc->tauMax[0]-pc->tauMin[0])/(float)n);
    std::cout << "C value " << i << " equals " << cValue[i] << std::endl;
  }
  
  c1 = 0.1;
  numPoints = 0;

  //findContours_Averaging(pc);
  //findContours(pc);
  find_Multiple_Contours(pc);

}
//void Contour::setContourValues(){

//}
void Contour::find_Multiple_Contours(ParticleControl* pc){
  for(int k=0;k<1;k++){
    for(int i=0; i < ny-1; i++){
      for(int j=0; j < nx-1; j++){
	int idx = k*ny*nx + i*nx + j;
	int idxAbove = k*ny*nx + (i+1)*nx + j;
	

	//Set Contour values
	//setContourValues();

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

	v0 = pc->tau[idx].t11;
	v1 = pc->tau[idx+1].t11;
	v2 = pc->tau[idxAbove].t11;
	v3 = pc->tau[idxAbove+1].t11;

	for(int c=0; c < num_cValue; c++){

	  localPoints=0;
	  //If isocontour value is between the tau values
	  if(((v0 <= cValue[c]) && (cValue[c] <= v1)) || ((v1 <= cValue[c]) && (cValue[c] <= v0))){
	    //Interpolate to find where value is cValue[c] at point p
	    vec2 p;
	    p.x = p0.x + ((cValue[c]-v0)/(v1-v0))*(p1.x-p0.x);
	    p.y = p0.y + ((cValue[c]-v0)/(v1-v0))*(p1.y-p0.y);
	    numPoints++;
	    localPoints++;
	    c1List.push_back(p);
	  }
	  if(((v0 <= cValue[c]) && (cValue[c] <= v2)) || ((v2 <= cValue[c]) && (cValue[c] <= v0))){
	    vec2 p;
	    p.x = p0.x + ((cValue[c]-v0)/(v2-v0))*(p2.x-p0.x);
	    p.y = p0.y + ((cValue[c]-v0)/(v2-v0))*(p2.y-p0.y);
	    numPoints++;
	    localPoints++;
	    c1List.push_back(p);
	  }
	  if(((v1 <= cValue[c]) && (cValue[c] <= v3)) || ((v3 <= cValue[c]) && (cValue[c] <= v1))){
	    vec2 p;
	    p.x = p1.x + ((cValue[c]-v1)/(v3-v1))*(p3.x-p1.x);
	    p.y = p1.y + ((cValue[c]-v1)/(v3-v1))*(p3.y-p1.y);
	    numPoints++;
	    localPoints++;
	    c1List.push_back(p);
	  }
	  if(((v2 <= cValue[c]) && (cValue[c] <= v3)) || ((v3 <= cValue[c]) && (cValue[c] <= v2))){
	    vec2 p;
	    p.x = p2.x + ((cValue[c]-v2)/(v3-v2))*(p3.x-p2.x);
	    p.y = p2.y + ((cValue[c]-v2)/(v3-v2))*(p3.y-p2.y);
	    numPoints++;
	    localPoints++;
	    c1List.push_back(p);
	  }
	  if(localPoints%2 != 0){
	    c1List.pop_back();
	    numPoints--;
	    std::cout << "There are " << localPoints << " points in a cell" << std::endl;
	  }
	}

      }
    }
  }
  std::cout << "Num of Points: " << numPoints << std::endl;
}

void Contour::findContours(ParticleControl* pc){
  for(int k=0;k<1;k++){
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

	v0 = pc->tau[idx].t11;
	v1 = pc->tau[idx+1].t11;
	v2 = pc->tau[idxAbove].t11;
	v3 = pc->tau[idxAbove+1].t11;

	localPoints=0;
	//If isocontour value is between the tau values
	if(((v0 <= c1) && (c1 <= v1)) || ((v1 <= c1) && (c1 <= v0))){
	  //Interpolate to find where value is c1 at point p
	  vec2 p;
	  p.x = p0.x + ((c1-v0)/(v1-v0))*(p1.x-p0.x);
	  p.y = p0.y + ((c1-v0)/(v1-v0))*(p1.y-p0.y);
	  numPoints++;
	  localPoints++;
	  c1List.push_back(p);
	}
	if(((v0 <= c1) && (c1 <= v2)) || ((v2 <= c1) && (c1 <= v0))){
	  vec2 p;
	  p.x = p0.x + ((c1-v0)/(v2-v0))*(p2.x-p0.x);
	  p.y = p0.y + ((c1-v0)/(v2-v0))*(p2.y-p0.y);
	  numPoints++;
	  localPoints++;
	  c1List.push_back(p);
	}
	if(((v1 <= c1) && (c1 <= v3)) || ((v3 <= c1) && (c1 <= v1))){
	  vec2 p;
	  p.x = p1.x + ((c1-v1)/(v3-v1))*(p3.x-p1.x);
	  p.y = p1.y + ((c1-v1)/(v3-v1))*(p3.y-p1.y);
	  numPoints++;
	  localPoints++;
	  c1List.push_back(p);
	}
	if(((v2 <= c1) && (c1 <= v3)) || ((v3 <= c1) && (c1 <= v2))){
	  vec2 p;
	  p.x = p2.x + ((c1-v2)/(v3-v2))*(p3.x-p2.x);
	  p.y = p2.y + ((c1-v2)/(v3-v2))*(p3.y-p2.y);
	  numPoints++;
	  localPoints++;
	  c1List.push_back(p);
	}
	if(localPoints%2 != 0){
	  c1List.pop_back();
	  numPoints--;
	  std::cout << "There are " << localPoints << " points in a cell" << std::endl;
	}

      }
    }
  }
  std::cout << "Num of Points: " << numPoints << std::endl;
}
void Contour::findContours_Averaging(ParticleControl* pc){
  float div;

  //will need to change it to k < nz
  for(int k=0;k<1;k++){
    for(int i=0; i < ny; i++){
      for(int j=0; j < nx; j++){
	int idx = k*ny*nx + i*nx + j;
	int idxAbove = k*ny*nx + (i+1)*nx + j;
	int idxBelow = k*ny*nx + (i-1)*nx + j;

	//v0 = pc->tau[idx].t11;
	//v1 = pc->tau[idx+1].t11;
	//v2 = pc->tau[idxAbove].t11;
	//v3 = pc->tau[idxAbove+1].t11;
	if(pc->tau[idx].t11 != 0.0){

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
	    v0 = pc->tau[idx].t11;
	  else if(p0.x==0 && p0.y != 0){
	    if(pc->tau[idxBelow].t11 == 0.0)
	      v0 = pc->tau[idx].t11;
	    else
	      v0 = (pc->tau[idx].t11 + pc->tau[idxBelow].t11)/2.0;
	  }
	  else if(p0.x!=0 && p0.y==0){
	    if(pc->tau[idx-1].t11 == 0.0)
	      v0 = pc->tau[idx].t11;
	    else
	      v0 = (pc->tau[idx].t11 + pc->tau[idx-1].t11)/2.0;
	  }
	  else{
	    div = 4.0;
	    if(pc->tau[idx-1].t11 == 0.0)
	      div -= 1.0;
	    if(pc->tau[idxBelow].t11 == 0)
	      div -= 1.0;
	    if(pc->tau[idxBelow-1].t11 == 0)
	      div -= 1.0;

	    v0 = (pc->tau[idx].t11 + pc->tau[idx-1].t11 
		  + pc->tau[idxBelow].t11 + pc->tau[idxBelow-1].t11)/div;
	  }
	  //Determine value at point p1 (lower right corner of cell)
	  if(p1.x == nx && p1.y==0)
	    v1 = pc->tau[idx].t11;
	  else if(p1.x != nx && p1.y==0){
	    if(pc->tau[idx+1].t11 == 0)
	      v1 = pc->tau[idx].t11;
	    else
	      v1 = (pc->tau[idx].t11 + pc->tau[idx+1].t11)/2.0;
	  }
	  else if(p1.x == nx && p1.y!=0){
	    if(pc->tau[idxBelow].t11 == 0)
	      v1 = pc->tau[idx].t11;
	    else
	      v1 = (pc->tau[idx].t11 + pc->tau[idxBelow].t11)/2.0;
	  }
	  else{
	    div = 4.0;
	    if(pc->tau[idx+1].t11 == 0)
	      div -= 1.0;
	    if(pc->tau[idxBelow].t11 == 0)
	      div -= 1.0;
	    if(pc->tau[idxBelow+1].t11 == 0)
	      div -= 1.0;
	    
	    v1 = (pc->tau[idx].t11 + pc->tau[idx+1].t11
		  + pc->tau[idxBelow].t11 + pc->tau[idxBelow+1].t11)/div;
	  }
	  //Determine value at point p2 (upper left corner of cell)
	  if(p2.x == 0 && p2.y == ny)
	    v2 = pc->tau[idx].t11;
	  else if(p2.x == 0 && p2.y != ny){
	    if(pc->tau[idxAbove].t11 == 0)
	      v2 = pc->tau[idx].t11;
	    else
	      v2 = (pc->tau[idx].t11 + pc->tau[idxAbove].t11)/2.0;
	  }
	  else if(p2.x != 0 && p2.y == ny){
	    if(pc->tau[idx-1].t11 == 0)
	      v2 = pc->tau[idx].t11;
	    else
	      v2 = (pc->tau[idx].t11 + pc->tau[idx-1].t11)/2.0;
	  }
	  else{
	    div = 4.0;
	    if(pc->tau[idx-1].t11 == 0)
	      div -= 1.0;
	    if(pc->tau[idxAbove].t11 == 0)
	      div -= 1.0;
	    if(pc->tau[idxAbove-1].t11 == 0)
	      div -= 1.0;

	    v2 = (pc->tau[idx].t11 + pc->tau[idx-1].t11 
		  + pc->tau[idxAbove].t11 + pc->tau[idxAbove-1].t11)/div;
	  }
	  //Determine value at point p3 (upper right corner of cell)
	  if(p3.x == nx && p3.y == ny)
	    v3 = pc->tau[idx].t11;
	  else if(p3.x ==nx && p3.y != ny){
	    if(pc->tau[idxAbove].t11 == 0)
	      v3 = pc->tau[idx].t11;
	    else
	      v3 = (pc->tau[idx].t11 + pc->tau[idxAbove].t11)/2.0;
	  }
	  else if(p3.x != nx && p3.y == ny){
	    if(pc->tau[idx+1].t11 == 0)
	      v3 = pc->tau[idx].t11;
	    else
	      v3 = (pc->tau[idx].t11 + pc->tau[idx+1].t11)/2.0;
	  }
	  else{
	    div = 4.0;
	    if(pc->tau[idx+1].t11 == 0)
	      div -= 1.0;
	    if(pc->tau[idxAbove].t11 == 0)
	      div -= 1.0;
	    if(pc->tau[idxAbove+1].t11 == 0)
	      div -= 1.0;

	    v3 = (pc->tau[idx].t11 + pc->tau[idx+1].t11 
		  + pc->tau[idxAbove].t11 + pc->tau[idxAbove+1].t11)/div;
	  }

	  //////////////////////////////////////////////
	  //Now check to see where the contour lines are
	  //////////////////////////////////////////////

	  //vec2 p;
	
	  localPoints=0;
	  if(((v0 <= c1) && (c1 <= v1)) || ((v1 <= c1) && (c1 <= v0))){
	    //Interpolate to find where value is c1 at point p
	    vec2 p;
	    p.x = p0.x + ((c1-v0)/(v1-v0))*(p1.x-p0.x);
	    p.y = p0.y + ((c1-v0)/(v1-v0))*(p1.y-p0.y);
	    numPoints++;
	    localPoints++;
	    c1List.push_back(p);
	  }
	  if(((v0 <= c1) && (c1 <= v2)) || ((v2 <= c1) && (c1 <= v0))){
	    vec2 p;
	    p.x = p0.x + ((c1-v0)/(v2-v0))*(p2.x-p0.x);
	    p.y = p0.y + ((c1-v0)/(v2-v0))*(p2.y-p0.y);
	    numPoints++;
	    localPoints++;
	    c1List.push_back(p);
	  }
	  if(((v1 <= c1) && (c1 <= v3)) || ((v3 <= c1) && (c1 <= v1))){
	    vec2 p;
	    p.x = p1.x + ((c1-v1)/(v3-v1))*(p3.x-p1.x);
	    p.y = p1.y + ((c1-v1)/(v3-v1))*(p3.y-p1.y);
	    numPoints++;
	    localPoints++;
	    c1List.push_back(p);
	  }
	  if(((v2 <= c1) && (c1 <= v3)) || ((v3 <= c1) && (c1 <= v2))){
	    vec2 p;
	    p.x = p2.x + ((c1-v2)/(v3-v2))*(p3.x-p2.x);
	    p.y = p2.y + ((c1-v2)/(v3-v2))*(p3.y-p2.y);
	    numPoints++;
	    localPoints++;
	    c1List.push_back(p);
	  }
	  if(localPoints%2 != 0){
	    c1List.pop_back();
	    numPoints--;
	    std::cout << "There are " << localPoints << " points in a cell" << std::endl;
	  }
	 

	}

      }
    }
  }
  std::cout << "Num of Points: " << numPoints << std::endl;

}
void Contour::draw(){
  
  listIter = c1List.begin();
  vec2 p1;
  vec2 p2;

  while(listIter != c1List.end()){
   
    p1 = *listIter;
    listIter++;
    p2 = *listIter;
    listIter++;
   

    glColor3f(1.0,1.0,1.0);
    glBegin(GL_LINES);
    {
      glVertex3f(p1.x,p1.y,0.2);
      glVertex3f(p2.x,p2.y,0.2);
    }
    glEnd();
    

  }

}
