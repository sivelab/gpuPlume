#include "streamLine.h"

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
}

void StreamLine::updateStreamPos(){
   glReadPixels(0, 0, twidth, theight, GL_RGBA, GL_FLOAT, pos_buffer); 
   partPos pos;
   streamIndex sIndex;

   indexIter = indexList.begin();
   while(indexIter != indexList.end()){
     sIndex = *indexIter;

     int idx = sIndex.t*twidth*4 + sIndex.s*4;
     pos.x = pos_buffer[idx];
     pos.y = pos_buffer[idx+1];
     pos.z = pos_buffer[idx+2];
       
     partPos prev = streamList[sIndex.i].back();
     //if((prev.x != pos.x) || (prev.y != pos.y) || (prev.z != pos.z))
     if((pos.x >= nx) || (pos.y >= ny) || (pos.z >= nz))
       sIndex.done = true;
     else
       streamList[sIndex.i].push_back(pos);

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

  for(int i =0;i < (int)streamList.size(); i++){
    for(int j = 0; j < (int)streamList[i].size(); j++){
      prevPos = pos;
      pos = streamList[i][j];

      if(!first){
	glBegin(GL_LINES);
	glColor3f(1.0,1.0,0.0);
	glVertex3f(prevPos.x,prevPos.y,prevPos.z);
	glVertex3f(pos.x,pos.y,pos.z);
	glEnd();    
      }
      if(first)
	first = false;
      
    }
    first = true;
  }
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
