#ifndef GETUSERSTATE_DATA_H
#define GETUSERSTATE_DATA_H

#include "Packet_Data.h"

// This class is for GetUserState packet data
//  There is nothing in the data, only the type matters
//  so all the functions (almost) can be null

class GetUserState_Data : public Packet_Data {
 public:
  // Destructor
  ~GetUserState_Data() {}
  
  // default constructor
  GetUserState_Data() { m_size = 0; }

  // copy constructor
  GetUserState_Data(const GetUserState_Data& r) : Packet_Data(r) {}

  // Assignment operator
  void operator = (const GetUserState_Data&) {}

  // Already provides uint GetSize();

 protected:
  // Inherits variable ushort m_size

 private:
  // None
};

#endif
