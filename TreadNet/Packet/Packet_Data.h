#ifndef PACKET_DATA_H
#define PACKET_DATA_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "TRVector.h"
#include "TRQuaternion.h"

class Packet_Data {
 public:
  // classes which derive from packet data should
  //  have the following

  // destructor
  virtual ~Packet_Data() {}
  
  // constructor
  Packet_Data() {m_size = 0;}

  // copy constructor (create for each child class)
  Packet_Data(const Packet_Data& d) {m_size = d.GetSize();}

  // Assignment operator (create for each child class)
  void operator = (const Packet_Data& d) {m_size = d.GetSize();}

  // Return size of data packet
  inline ushort GetSize() const {return m_size;}

  // Access functions for data

 protected:
  // Data specific to this class
  
  // Size of packet data in # of bytes
  ushort m_size;
};

#endif
