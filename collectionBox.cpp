#include "collectionBox.h"
#include <iostream>

CollectionBox::CollectionBox(int x,int y,int z,float* bounds,float time){
  numBox_x = x;
  numBox_y = y;
  numBox_z = z;
  
  lx = bounds[0];
  ly = bounds[1];
  lz = bounds[2];
  ux = bounds[3];
  uy = bounds[4];
  uz = bounds[5];

  volume = (ux-lx)*(uy-ly)*(uz-lz);
  TotRel = 1.0;
  concAvgTime = time;
  
  cBox = new cell[x*y*z];
  for(int i = 0; i < numBox_x*numBox_y*numBox_z; i++){
    cBox[i].concentration = 0.0;    
    cBox[i].count = 0;
  }

  alreadyOpen = false;
  
}
/*void CollectionBox::seeIfInBox(float x, float y, float z){

  if( (lx <= x)&&(x <= ux)&&(ly <= y)&&(y <= uy)&&(lz <= z)&&(z <= uz) ){
   
    int xBox = (int)floor((x-lx)/((ux-lx)/numBox_x));
    int yBox = (int)floor((y-ly)/((uy-ly)/numBox_y));
    int zBox = (int)floor((z-lz)/((uz-lz)/numBox_z));
    int idx = yBox*numBox_x*numBox_z + xBox*numBox_z + zBox;
    //cBox[idx].count++; 
    
    
  }
  }*/
void CollectionBox::calculateConc(float x,float y,float z,float timeStep,double totalNumPar){

   if( (lx <= x)&&(x <= ux)&&(ly <= y)&&(y <= uy)&&(lz <= z)&&(z <= uz) ){
   
    int xBox = (int)floor((x-lx)/((ux-lx)/numBox_x));
    int yBox = (int)floor((y-ly)/((uy-ly)/numBox_y));
    int zBox = (int)floor((z-lz)/((uz-lz)/numBox_z));
    int idx = yBox*numBox_x*numBox_z + xBox*numBox_z + zBox;
    
    constant = (timeStep*TotRel)/(volume*concAvgTime*totalNumPar);
    cBox[idx].concentration = cBox[idx].concentration + constant;
    
  }

}
void CollectionBox::clear(){
  for(int i=0; i < numBox_x*numBox_y*numBox_z; i++){
    cBox[i].concentration = 0;
  }

}
void CollectionBox::outputConc(std::string file,double totalTime){
  std::ofstream output;


  if(!alreadyOpen){
    output.open(file.c_str());
    alreadyOpen = true;
  }
  else{
    output.open(file.c_str(),std::ios::app);
  }

  output << "Variables: X Y Z C" << "\n";
  output << "Average Time: " << totalTime << "\n";
  //std::cout << totalTime << std::endl;
 
  for(int k=0; k < numBox_y; k++)
    for(int i=0; i < numBox_x; i++)
      for(int j=0; j < numBox_z; j++)
	{
	  int idx = k*numBox_z*numBox_x + i*numBox_z + j;
	  int xBox = idx%numBox_x;
	  int yBox = idx/(numBox_x*numBox_z);
	  int zBox = (idx/numBox_x)%numBox_z;
	  
	  float x = ((ux-lx)/(float)numBox_x)*xBox + lx;
	  float y = ((uy-ly)/(float)numBox_y)*yBox + ly;
	  float z = ((uz-lz)/(float)numBox_z)*zBox + lz;
	  float offsetx = ((ux-lx)/(float)numBox_x)/2.0;
	  float offsety = ((uy-ly)/(float)numBox_y)/2.0;
	  float offsetz = ((uz-lz)/(float)numBox_z)/2.0;
	  
	  output << z+offsetz << "  " << x+offsetx << "  " << y+offsety <<
	  "  " << cBox[idx].concentration << "\n";

	  /*std::cout << z+offsetz << "  " << x+offsetx << "  " << y+offsety <<
	    "  " << cBox[idx].concentration << std::endl;*/

	}

}
