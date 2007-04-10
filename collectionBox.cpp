#include "collectionBox.h"
#include <iostream>

CollectionBox::CollectionBox(int x,int y,int z,float* bounds,double time){
  numBox_x = x;
  numBox_y = y;
  numBox_z = z;
  
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

  clear();

  alreadyOpen = false;
}

void CollectionBox::calculateConc(float x,float y,float z,float timeStep,double totalNumPar){

   if( (lx <= x)&&(x <= ux)&&(ly <= y)&&(y <= uy)&&(lz <= z)&&(z <= uz) ){
        
    int xBox = (int)floor((x-lx)/((ux-lx)/(float)numBox_x));
    int yBox = (int)floor((y-ly)/((uy-ly)/(float)numBox_y));
    int zBox = (int)floor((z-lz)/((uz-lz)/(float)numBox_z));
    int idx = yBox*numBox_x*numBox_z + xBox*numBox_z + zBox;
       
    constant = (timeStep*TotRel)/(volume*concAvgTime*totalNumPar);
    cBox[idx] = cBox[idx] + constant;   
  }

}
void CollectionBox::clear(){
  for(int i=0; i < numBox_x*numBox_y*numBox_z; i++){
    cBox[i] = 0.0; 
  }

}
void CollectionBox::outputConc(std::string file,double totalTime){
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
 
  for(int k=0; k < numBox_y; k++)
    for(int i=0; i < numBox_x; i++)
      for(int j=0; j < numBox_z; j++)
	{
	  idx = k*numBox_z*numBox_x + i*numBox_z + j;
	  xBox = idx%numBox_x;
	  yBox = idx/(numBox_x*numBox_z);
	  zBox = (idx/numBox_x)%numBox_z;
	  
	  x = ((ux-lx)/(float)numBox_x)*xBox + lx;
	  y = ((uy-ly)/(float)numBox_y)*yBox + ly;
	  z = ((uz-lz)/(float)numBox_z)*zBox + lz;
	  offsetx = ((ux-lx)/(float)numBox_x)/2.0;
	  offsety = ((uy-ly)/(float)numBox_y)/2.0;
	  offsetz = ((uz-lz)/(float)numBox_z)/2.0;
	  
	  output << z+offsetz << "  " << x+offsetx << "  " << y+offsety <<
	    "  " << cBox[idx] << "\n";

	  /*std::cout << z+offsetz << "  " << x+offsetx << "  " << y+offsety <<
	    "  " << cBox[idx] << std::endl;*/

	}

  clear();

}
