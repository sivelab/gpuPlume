#include "collectionBox.h"
#include <iostream>

CollectionBox::CollectionBox(Util* util){
  numBox_x = util->numBox_x;
  numBox_y = util->numBox_y;
  numBox_z = util->numBox_z;

  twidth = util->twidth;
  theight = util->theight;
  
  lx = util->bounds[0];
  ly = util->bounds[1];
  lz = util->bounds[2];
  ux = util->bounds[3];
  uy = util->bounds[4];
  uz = util->bounds[5];

  volume = (double)(ux-lx)*(uy-ly)*(uz-lz)/(double)(numBox_x*numBox_y*numBox_z);
  TotRel = (double)1.0;
  concAvgTime = util->averagingTime;


  avgTime = util->averagingTime + util->startCBoxTime;
  //avgTime = util->endCBoxTime;
  endCBoxTime = util->endCBoxTime;
  startCBoxTime = util->startCBoxTime;

  cBox = new double[numBox_x*numBox_y*numBox_z];
  cPos = new BoxPos[numBox_x*numBox_y*numBox_z];
  cDis = new BoxDis[numBox_x*numBox_y*numBox_z];
  calcBoxPositions();

  //This is used for drawing the collection boxes.
  //It is the distance from center to edge, minus
  //a small spacing.
  dx = float(((ux-lx)/(float)numBox_x)/2.0 - 0.1);
  dy = float(((uy-ly)/(float)numBox_y)/2.0 - 0.1);
  dz = float(((uz-lz)/(float)numBox_z)/2.0 - 0.1);

  pos_buffer = new GLfloat[ twidth * theight * 4 ];

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
	  offsetx = float(((ux-lx)/(float)numBox_x)/2.0);
	  offsety = float(((uy-ly)/(float)numBox_y)/2.0);
	  offsetz = float(((uz-lz)/(float)numBox_z)/2.0);
	  
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
  float temp;
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
       alpha = float((cBox[idx]/timeStepNum)/10.0);

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
bool CollectionBox::findConc(Simulation* sim, bool* endCBox,bool odd){
  //if(sim->totalTime >= endCBoxTime){

  if(fabs(sim->totalTime - endCBoxTime) <= 0.000001){
     //sim->curr_timeStep += 1.0;
     *endCBox = true;
     if(endCBoxTime != 0)
       return true;
   }  
   if((sim->totalTime >= startCBoxTime) && !*endCBox){

     //sim->curr_timeStep +=1.0;

     if(odd)
       glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
     else
       glReadBuffer(GL_COLOR_ATTACHMENT1_EXT);


     glReadPixels(0, 0, twidth, theight, GL_RGBA, GL_FLOAT, pos_buffer); 

     for(int i = 3; i <= (theight*twidth*4); i+=4){
       //If particle has been emitted
       if(pos_buffer[i] == -1){
	
	 //Get the x,y,z position of the particle
	 float x = pos_buffer[i-3];
	 float y = pos_buffer[i-2];
	 float z = pos_buffer[i-1];

	 //Check to see if particle is inside a collection box
	 //if a particle is in a box the concentration value is updated
	 //cBoxes[j]->calculateConc(x,y,z,time_step,totalNumPar);
	 calcSimpleConc(x,y,z);	
       }
     }    

     if(sim->totalTime >= avgTime){
       avgTime += concAvgTime;
       return true;
     }

   }
   return false;

}

void CollectionBox::calcSimpleConc(float x, float y, float z){
  if( (lx <= x)&&(x < ux)&&(ly <= y)&&(y < uy)&&(lz <= z)&&(z < uz) ){   

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
    std::cout << "File Name: " << file.c_str() << std::endl;
    alreadyOpen = true;
  }
  else{
    output.open(file.c_str(),std::ios::app);
  }

  // The format of the output file has been changed to work directly
  // with matlab. -Pete

  output << "average_time = " << totalTime << ";\n";
  output << "total_time_steps = " << totalTimeSteps << ";\n\n";

  output << "% The following array contains the locations of the\n";
  output << "% collection box cells in X, Y, and Z followed by the \n";
  output << "% concentration in the cell." << std::endl;
  output << "concentration = [\n";

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
	  offsetx = float(((ux-lx)/(float)numBox_x)/2.0);
	  offsety = float(((uy-ly)/(float)numBox_y)/2.0);
	  offsetz = float(((uz-lz)/(float)numBox_z)/2.0);
	  
	  output << '\t' << x+offsetx << ' ' << y+offsety << ' ' << z+offsetz 
		 << ' ' << cBox[idx]/totalTimeSteps << ";\n";

	  /*std::cout << z+offsetz << "  " << x+offsetx << "  " << y+offsety <<
	    "  " << cBox[idx] << std::endl;*/
	}

  output << "];\n";
}

void CollectionBox::outputConcStd(std::string file,double averagingTime,double volume, float time_step, int numPar){ //standard Concentration Calc. - Balli
  std::ofstream output;

  int idx,xBox,yBox,zBox;
  float x,y,z,offsetx,offsety,offsetz;

  if(!alreadyOpen){
    output.open(file.c_str());
    std::cout << "File Name: " << file.c_str() << std::endl;
    alreadyOpen = true;
  }
  else{
    output.open(file.c_str(),std::ios::app);
  }

  // The format of the output file has been changed to work directly
  // with matlab. -Pete

  output << "average_time = " << averagingTime << ";\n";
  //output << "total_time_steps = " << totalTimeSteps << ";\n\n";

  output << "% The following array contains the locations of the\n";
  output << "% collection box cells in X, Y, and Z followed by the \n";
  output << "% concentration in the cell." << std::endl;
  output << "concentration = [\n";

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
	  offsetx = float(((ux-lx)/(float)numBox_x)/2.0);
	  offsety = float(((uy-ly)/(float)numBox_y)/2.0);
	  offsetz = float(((uz-lz)/(float)numBox_z)/2.0);
	  float concFactor= (time_step)/(averagingTime*volume*numPar);
	  
	  output << '\t' << x+offsetx << ' ' << y+offsety << ' ' << z+offsetz 
		 << ' ' << cBox[idx] *(concFactor) << ";\n";

	  /*std::cout << z+offsetz << "  " << x+offsetx << "  " << y+offsety <<
	    "  " << cBox[idx] << std::endl;*/
	}

  output << "];\n";
}