#include "streamLine.h"
#include <math.h>

StreamLine::StreamLine(int w, int h,int x,int y,int z){
  twidth = w;
  theight = h;
  nx = x;
  ny = y;
  nz = z;

  streamNum = 0;
  startStream = false;
  pos_buffer = new GLfloat[ twidth * theight * 4 ];
  update = true;
}

void StreamLine::addNewStream(ParticleEmitter* pe){
   std::vector<partPos> posList;
   partPos p;
   int xindex;
   int yindex;
   streamIndex sIndex;

   pe->getReleasedPosition(&p.x,&p.y,&p.z);
   posList.push_back(p);
   pe->getIndex(&xindex,&yindex);
   sIndex.s = xindex;
   sIndex.t = yindex;
   sIndex.i = streamNum;
   sIndex.done = false;
   indexList.push_front(sIndex);
     
   streamList.push_back(posList);
   streamNum++;
   
   startStream = true;
   update = true;

}

void StreamLine::updateStreamPos(){
   glReadPixels(0, 0, twidth, theight, GL_RGBA, GL_FLOAT, pos_buffer); 
   partPos pos;
   streamIndex sIndex;

   indexIter = indexList.begin();
   while(indexIter != indexList.end()){
     sIndex = *indexIter;
     if(!sIndex.done){

       streamIndex &sIdx = *indexIter;

       int idx = sIndex.t*twidth*4 + sIndex.s*4;
       pos.x = pos_buffer[idx];
       pos.y = pos_buffer[idx+1];
       pos.z = pos_buffer[idx+2];
       
       partPos prev = streamList[sIndex.i].back();
       //if((prev.x != pos.x) || (prev.y != pos.y) || (prev.z != pos.z))
       if((pos.x >= nx) || (pos.y >= ny) || (pos.z >= nz)
	  || (pos.x < 0) || (pos.y < 0) || (pos.z < 0)){

	 sIdx.done = true;
       }
       else
	 streamList[sIndex.i].push_back(pos);

     }
     indexIter++;
   }
     
   update = true;
}

void StreamLine::draw(){

  partPos prevPos;
  partPos pos;
  bool first = true;

  //This is just so I don't get a warning after compiling
  pos.x = pos.y = pos.z = 0.0;
  
  glLineWidth(3.0);
  //glEnable(GL_LINE_SMOOTH);
  //glEnable(GL_BLEND);
  //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  //glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  

  for(int i =0;i < (int)streamList.size(); i++){
    for(int j = 0; j < (int)streamList[i].size(); j++){
      prevPos = pos;
      pos = streamList[i][j];


      colorx = fabs(prevPos.x - pos.x)*100;
      colory = fabs(prevPos.y-pos.y)*100;
      colorz = fabs(prevPos.z-pos.z)*100;
     //glColor4f(1.0,1.0,0.0,1.0);

      if(!first){
	glColor4f(colorx,colory,colorz,1.0);

	glBegin(GL_LINES);
	{		
	  glVertex3f(prevPos.x,prevPos.y,prevPos.z);
	  glVertex3f(pos.x,pos.y,pos.z);
	}
	glEnd();    
      }
      if(first)
	first = false;
      
    }
    first = true;
  }
  //glDisable(GL_BLEND);
  //glDisable(GL_LINE_SMOOTH);
  //glLineWidth(1.0);
  

}
bool StreamLine::doUpdate(){
  if(indexList.empty())
    return false;

  streamIndex sIndex;
  if(update){
    indexIter = indexList.begin();
    while(indexIter != indexList.end()){
      sIndex = *indexIter;
      if(!sIndex.done){
	return update;
      }
    
      indexIter++;
    }

  }
  update = false;
  return update;

}
void StreamLine::clear(){
  streamList.clear();
  indexList.clear();
  startStream = true;
  update = true;
  streamNum = 0;

}
