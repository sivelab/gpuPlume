#ifndef NEWPOS_DATA_H
#define NEWPOS_DATA_H

#include "VirtualState_Data.h"
#include "UserState_Data.h"

// This class is for New Position packet data
// Essentially it holds both a virtual state and user state
//  data class

class NewPos_Data : public Packet_Data {
 public:
  // Destructor
  ~NewPos_Data();
  
  // default constructor
  NewPos_Data();

  // copy constructor
  NewPos_Data(const NewPos_Data&);

  // Constructor with UserState and VirtualState passed in
  NewPos_Data(const UserState_Data&, const VirtualState_Data&);

  // Assignment operator
  void operator = (const NewPos_Data&);

  // Accessor functions
  inline UserState_Data GetUserState() const {return m_user_state;}
  inline VirtualState_Data GetVirtualState() const {return m_virtual_state;}

  // Set functions
  inline void SetUserState(const UserState_Data& usd) {m_user_state = usd;}
  inline void SetVirtualState(const VirtualState_Data& vsd) {m_virtual_state = vsd;}

  // Already provides uint GetSize();

 protected:
  // Inherits variable ushort m_size

 private:
  // User state information
  UserState_Data m_user_state;

  // Virtual state information
  VirtualState_Data m_virtual_state;
};

#endif
