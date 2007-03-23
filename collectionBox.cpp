#include "collectionBox.h"

CollectionBox::CollectionBox(int x,int y,int z,float* bounds){
  numBox_x = x;
  numBox_y = y;
  numBox_z = z;
  
  lx = bounds[0];
  ux = bounds[1];
  ly = bounds[2];
  uy = bounds[3];
  lz = bounds[4];
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
  alpha = n*(n/1.0);

  for(int i = 0; i < numBox_x*numBox_y*numBox_z; i++){
    cBox[i].moving_avg = alpha*cBox[i].moving_avg + (1-alpha)*cBox[i].count;
    
    cBox[i].count = 0;

  }
  n++;

}
