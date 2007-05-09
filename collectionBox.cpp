#include "collectionBox.h"
#include <iostream>

CollectionBox::CollectionBox(int x,int y,int z,float* bounds,double time){
  numBox_x = x;
  numBox_y = y;
  numBox_z = z;
  
  //Converts from quic-plume coordinate system to openGL: Don't do anymore!
  /*lx = bounds[1];
  ly = bounds[2];
  lz = bounds[0];
  ux = bounds[4];
  uy = bounds[5];
  uz = bounds[3];*/
  lx = bounds[0];
  ly = bounds[1];
  lz = bounds[2];
  ux = bounds[3];
  uy = bounds[4];
  uz = bounds[5];

  volume = (double)(ux-lx)*(uy-ly)*(uz-lz);
  TotRel = (double)1.0;
  concAvgTime = time;

  cBox = new double[x*y*z];
  cPos = new BoxPos[x*y*z];
  cDis = new BoxDis[x*y*z];
  calcBoxPositions();

  //This is used for drawing the collection boxes.
  //It is the distance from center to edge, minus
  //a small spacing.
  dx = ((ux-lx)/(float)numBox_x)/2.0 - 0.1;
  dy = ((uy-ly)/(float)numBox_y)/2.0 - 0.1;
  dz = ((uz-lz)/(float)numBox_z)/2.0 - 0.1;


  clear();

  alreadyOpen = false;
  
}
void CollectionBox::calcBoxPositions(){
  int idx,xBox,yBox,zBox;
  float x,y,z,offsetx,offsety,offsetz;

  for(int k=0; k < numBox_z; k++)
    for(int i=0; i < numBox_y; i++)
      for(int j=0; j < numBox_x; j++)
	{
	  idx = k*numBox_x*numBox_y + i*numBox_x + j;
	  xBox = idx%numBox_x;
	  yBox = (idx/numBox_x)%numBox_y;
	  zBox = idx/(numBox_y*numBox_x);
	  
	  
	  x = ((ux-lx)/(float)numBox_x)*xBox + lx;
	  y = ((uy-ly)/(float)numBox_y)*yBox + ly;
	  z = ((uz-lz)/(float)numBox_z)*zBox + lz;
	  offsetx = ((ux-lx)/(float)numBox_x)/2.0;
	  offsety = ((uy-ly)/(float)numBox_y)/2.0;
	  offsetz = ((uz-lz)/(float)numBox_z)/2.0;
	  
	  cPos[idx].x = x+offsetx;
	  cPos[idx].y = y+offsety;
	  cPos[idx].z = z+offsetz;
	  cDis[idx].idx = idx;
	 
	}

}
void CollectionBox::sort(float eyeX, float eyeY, float eyeZ){
  float xdiff, ydiff, zdiff;
  int size = numBox_z*numBox_y*numBox_x;
  //int max;
  double temp;
  int t;

  //Calculate and store distance from eye to each box
  for(int i=0; i < size; i++){
    xdiff = fabs(cPos[i].x - eyeX);
    ydiff = fabs(cPos[i].y - eyeY);
    zdiff = fabs(cPos[i].z - eyeZ);

    cDis[i].d = sqrt((xdiff*xdiff)+(ydiff*ydiff)+(zdiff*zdiff));
  }
  //Sort Distances from farthest to closest.
  for(int i=0; i < size-1; i++){
    for(int j=0; j < size-1-i; j++){
      if(cDis[j+1].d > cDis[j].d){
	 temp = cDis[j].d;
	 cDis[j].d = cDis[j+1].d;
	 cDis[j+1].d = temp;
	 t = cDis[j].idx;
	 cDis[j].idx = cDis[j+1].idx;
	 cDis[j].idx = t;
      }
    }
   
   
  }
 

}

void CollectionBox::draw(double timeStepNum){
  float x,y,z;
  float alpha;
  float colorx,colory,colorz;

  glDisable(GL_TEXTURE_RECTANGLE_ARB);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


  for(int i=0; i < numBox_z*numBox_y*numBox_x; i++){
    int idx = cDis[i].idx;
    
     x = cPos[idx].x;
     y = cPos[idx].y;
     z = cPos[idx].z;
     
     colorx=1.0;
     colory=0.0;
     colorz=0.0;

     if(timeStepNum == 0)
       alpha = 0.0;
     else
       alpha = (cBox[idx]/(timeStepNum))/100.0;

     glBegin(GL_QUADS);
     {
	glColor4f(colorx,colory,colorz,alpha);
	glVertex3f(x+dx,y-dy,z-dz);
	glVertex3f(x+dx,y+dy,z-dz);
	glVertex3f(x+dx,y+dy,z+dz);
	glVertex3f(x+dx,y-dy,z+dz);

	glVertex3f(x-dx,y+dy,z-dz);
	glVertex3f(x-dx,y+dy,z+dz);
	glVertex3f(x-dx,y-dy,z+dz);
	glVertex3f(x-dx,y-dy,z-dz);

	glVertex3f(x+dx,y+dy,z-dz);
	glVertex3f(x-dx,y+dy,z-dz);
	glVertex3f(x-dx,y+dy,z+dz);
	glVertex3f(x+dx,y+dy,z+dz);

	glVertex3f(x+dx,y-dy,z-dz);
	glVertex3f(x+dx,y-dy,z+dz);
	glVertex3f(x-dx,y-dy,z+dz);
	glVertex3f(x-dx,y-dy,z-dz);
	      
	glVertex3f(x+dx,y-dy,z+dz);
	glVertex3f(x+dx,y+dy,z+dz);
	glVertex3f(x-dx,y+dy,z+dz);
	glVertex3f(x-dx,y-dy,z+dz);

	glVertex3f(x+dx,y-dy,z-dz);
	glVertex3f(x+dx,y+dy,z-dz);
	glVertex3f(x-dx,y+dy,z-dz);
	glVertex3f(x-dx,y-dy,z-dz);	    

      }
      glEnd();
  }

  glDisable(GL_BLEND);
  glDisable(GL_COLOR_MATERIAL);
  glEnable(GL_TEXTURE_RECTANGLE_ARB);

}

void CollectionBox::calcSimpleConc(float x, float y, float z){
  if( (lx <= x)&&(x <= ux)&&(ly <= y)&&(y <= uy)&&(lz <= z)&&(z <= uz) ){
        
    int xBox = (int)floor((x-lx)/((ux-lx)/(float)numBox_x));
    int yBox = (int)floor((y-ly)/((uy-ly)/(float)numBox_y));
    int zBox = (int)floor((z-lz)/((uz-lz)/(float)numBox_z));
    int idx = zBox*numBox_y*numBox_x + yBox*numBox_x + xBox;
       
    cBox[idx] = cBox[idx] + 1.0;   
  }

}

void CollectionBox::calculateConc(float x,float y,float z,float timeStep,double totalNumPar){

   if( (lx <= x)&&(x <= ux)&&(ly <= y)&&(y <= uy)&&(lz <= z)&&(z <= uz) ){
        
    int xBox = (int)floor((x-lx)/((ux-lx)/(float)numBox_x));
    int yBox = (int)floor((y-ly)/((uy-ly)/(float)numBox_y));
    int zBox = (int)floor((z-lz)/((uz-lz)/(float)numBox_z));
    int idx = zBox*numBox_x*numBox_y + yBox*numBox_x + xBox;
       
    constant = (timeStep*TotRel)/(volume*concAvgTime*totalNumPar);
    cBox[idx] = cBox[idx] + constant;   
  }

}
void CollectionBox::clear(){
  for(int i=0; i < numBox_x*numBox_y*numBox_z; i++){
    cBox[i] = 0.0; 
  }

}
void CollectionBox::outputConc(std::string file,double totalTime,double totalTimeSteps){
  std::ofstream output;

  int idx,xBox,yBox,zBox;
  float x,y,z,offsetx,offsety,offsetz;

  if(!alreadyOpen){
    output.open(file.c_str());
    alreadyOpen = true;
  }
  else{
    output.open(file.c_str(),std::ios::app);
  }

  output << "Variables: X Y Z C" << "\n";
  output << "Average Time: " << totalTime << "\n";
  //std::cout << "Variables: X Y Z C" << std::endl;
  //std::cout << "Average Time: " << totalTime << std::endl;
 
  for(int k=0; k < numBox_z; k++)
    for(int i=0; i < numBox_y; i++)
      for(int j=0; j < numBox_x; j++)
	{
	  idx = k*numBox_x*numBox_y + i*numBox_x + j;
	  xBox = idx%numBox_x;
	  yBox = (idx/numBox_x)%numBox_y;
	  zBox = idx/(numBox_y*numBox_x);
	  
	  
	  x = ((ux-lx)/(float)numBox_x)*xBox + lx;
	  y = ((uy-ly)/(float)numBox_y)*yBox + ly;
	  z = ((uz-lz)/(float)numBox_z)*zBox + lz;
	  offsetx = ((ux-lx)/(float)numBox_x)/2.0;
	  offsety = ((uy-ly)/(float)numBox_y)/2.0;
	  offsetz = ((uz-lz)/(float)numBox_z)/2.0;
	  
	  output << x+offsetx << "  " << y+offsety << "  " << z+offsetz <<
	    "  " << cBox[idx]/totalTimeSteps << "\n";

	  /*std::cout << z+offsetz << "  " << x+offsetx << "  " << y+offsety <<
	    "  " << cBox[idx] << std::endl;*/

	}

  //clear();

}
