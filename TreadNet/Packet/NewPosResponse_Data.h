#ifndef NEWPOSRESPONSE_DATA_H
#define NEWPOSRESPONSE_DATA_H

#include "Packet_Data.h"

enum TRNewPosFlag {
  NEW_POS_FAIL,
  NEW_POS_SUCCESS
};

// This class is for NewPosResponse packet's data

class NewPosResponse_Data : public Packet_Data {
 public:
  friend class NewPosResponse_Packet;

  // Destructor
  ~NewPosResponse_Data() {}
  
  // default constructor
  NewPosResponse_Data() {
    m_flag = NEW_POS_FAIL;
    m_size = sizeof(char); // using char
  }

  // copy constructor
  inline NewPosResponse_Data(const NewPosResponse_Data& right) : Packet_Data(right) {
    m_flag = right.GetFlag();
    m_size = sizeof(char); // using char
  }

  // Assignment operator
  inline void operator = (const NewPosResponse_Data& right) {
    m_flag = right.GetFlag();
  }

  // Accessor functions
  inline TRNewPosFlag GetFlag() const { return m_flag; }

  // Set functions
  inline void SetSuccessFlag() {m_flag = NEW_POS_SUCCESS;}
  inline void SetFailFlag() {m_flag = NEW_POS_FAIL;}

  // Already provides uint GetSize();

 protected:
  // Inherits variable ushort m_size
  inline void SetFlag(TRNewPosFlag flag) {m_flag = flag;}

 private:
  // NewPosResponse command
  TRNewPosFlag m_flag;
};

#endif
