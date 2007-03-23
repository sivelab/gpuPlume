#include "collectionBox.h"

CollectionBox::CollectionBox(int x,int y,int z,float* bounds){
  numBox_x = x;
  numBox_y = y;
  numBox_z = z;
  
  lx = bounds[0];
  ly = bounds[1];
  lz = bounds[2];
  ux = bounds[3];
  uy = bounds[4];
  uz = bounds[5];
  
  cBox = new cell[x*y*z];
  for(int i = 0; i < numBox_x*numBox_y*numBox_z; i++){
    cBox[i].moving_avg = 0;    
    cBox[i].count = 0;
  }
  n = 0;

  
}
void CollectionBox::seeIfInBox(float x, float y, float z){

  if( (lx <= x)&&(x <= ux)&&(ly <= y)&&(y <= uy)&&(lz <= z)&&(z <= uz) ){
   
    int xBox = (int)floor((x-lx)/((ux-lx)/numBox_x));
    int yBox = (int)floor((y-ly)/((uy-ly)/numBox_y));
    int zBox = (int)floor((z-lz)/((uz-lz)/numBox_z));
    int idx = yBox*numBox_x*numBox_z + xBox*numBox_z + zBox;
    cBox[idx].count++;     
  }

}
void CollectionBox::calculateAvg(){
  alpha = n/(n+1.0);
  //std::cout << alpha << std::endl;

  for(int i = 0; i < numBox_x*numBox_y*numBox_z; i++){
    cBox[i].moving_avg = alpha*cBox[i].moving_avg + (1.0-alpha)*cBox[i].count;
    
    cBox[i].count = 0;

  }
  n++;

}
void CollectionBox::outputAvg(){
   for(int i = 0; i < numBox_x*numBox_y*numBox_z; i++){
     std::cout << cBox[i].moving_avg << std::endl;

  }

}
