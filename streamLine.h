#ifndef __STREAMLINE_H__
#define __STREAMLINE_H__


#include <iostream>
#include <list>
#include <vector>
#include "particleEmitter.h"

typedef struct{
    int s;
    int t;
    int i;
    bool done;
}streamIndex;

typedef struct{
  float x;
  float y;
  float z;
}partPos;

class StreamLine{
 public:

  StreamLine(int,int,int,int,int);

  void addNewStream(ParticleEmitter*);

  void updateStreamPos();

  void draw();

  bool doUpdate();

  int streamNum;
  bool startStream;

 private:
  GLfloat* pos_buffer;
  //A list of the particle position lists
  std::vector<std::vector<partPos> > streamList;
  //A list of indices for streams
  std::list<streamIndex> indexList;

  std::list<streamIndex>::iterator indexIter;

  int twidth;
  int theight;
  int nx;
  int ny;
  int nz;

  bool update;

};

#endif //__STREAMLINE_H__
