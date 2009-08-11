#ifndef BYTEARRAY_DATA_H
#define BYTEARRAY_DATA_H

#include "Packet_Data.h"

#define ARRAY_SIZE 500

// This class is for sending byte arrays, which in
//  turn can be anything (you could store doubles, floats, etc. as
//  well as chars, bits).  One thing is make sure you put
//  any doubles/floats into network byte order.

class ByteArray_Data : public Packet_Data {
 public:
  friend class ByteArray_Packet;

  // Destructor
  ~ByteArray_Data();
  
  // default constructor
  ByteArray_Data();

  // copy constructor
  ByteArray_Data(const ByteArray_Data&);

  // Assignment operator
  void operator = (const ByteArray_Data&);

  // Accessor functions

  // Set ByteArray element
  int SetByte(int index, char value);

  // Get functions (this will actually be the most used for
  //  building the array - use the ptr to set variables larger
  //  than single byte - see Packet.cpp for examples).
  char * GetArrayPtr(int index = 0);

  // 2nd version is for read only 
  const char * GetArrayPtr(int index = 0) const;

  // Sets the last element index since no way to know for sure
  //  how many bytes to send.  w/o this the default size will be
  //  0 (or ind -1) - i.e. nothing.  Returns 0 if index out of bounds.
  int SetLastElementInd(int val);

  inline int GetLastElementInd() const { return m_LastElementInd; }
  
  // Already provides uint GetSize();

 protected:
  // Inherits variable ushort m_size

 private:
  // Array of bytes.  Setting a specific size so can match up
  //  between server/client side.  Also allows bounds checking to ensure
  //  using array correctly.
  char m_array[ARRAY_SIZE];

  // last element index, could use size, but it would be complicated
  //  because m_size now also has the 2 bytes for how large the array
  //  is since it's supposed to represent how many bytes this packet
  //  takes up in the message (and it sends the array size w/ it).
  int m_LastElementInd;
};

#endif
