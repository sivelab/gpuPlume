#ifndef REGISTRATION_DATA_H
#define REGISTRATION_DATA_H

#include "Packet_Data.h"

enum TRRegTag {
  UNKNOWN_CLASS,
  INVALID_REG_SIZE,
  VALID_REG
};

// This class is for registering packet types
//  with the controller (whether server or client)

class Registration_Data : public Packet_Data {
 public:
  friend class Registration_Packet;

  // Destructor
  ~Registration_Data();
  
  // default constructor
  Registration_Data();

  // copy constructor
  Registration_Data(const Registration_Data&);

  // Assignment operator
  void operator = (const Registration_Data&);

  // Accessor functions

  // Set Registration Class Name
  int SetRegClassName(const char *);

  // Set Registration Class's Data Size (in bytes)
  inline void SetRegClassSize(ushort val) {m_regclass_data_size = val;}

  //Set Registration Class's Bit type
  inline void SetRegClassType(uint val) {m_regclass_type = val;}

  // Set Tag values VALID_REG, UNKNOWN_CLASS, INVALID_REG_SIZE
  inline void SetValidRegTag() {m_tag = VALID_REG;}
  inline void SetUnknownClassRegTag() {m_tag = UNKNOWN_CLASS;}
  inline void SetIncorrectSizeRegTag() {m_tag = INVALID_REG_SIZE;}

  // Get functions
  inline char * GetRegClassName() const {return m_regclass_name;}
  inline ushort GetRegClassSize() const {return m_regclass_data_size;}
  inline uint GetRegClassType() const {return m_regclass_type;}
  inline TRRegTag GetRegTag() const {return m_tag;}

  // Already provides uint GetSize();

 protected:
  // Do not want user to set any tag
  inline void SetRegTag(TRRegTag tag) {m_tag = tag;}

  // Inherits variable ushort m_size

 private:
  // Function for calculating how many bytes this data packet
  //  will take up and sets m_size accordingly.
  void CalcDataSize();

  // Main part of data is the class_name to register
  char *m_regclass_name;

  // Registration Class's Data size for packet m_name.
  //  This is used to verify it is the correct registration
  //  (i.e. client and server do not have different Packet class
  //  descriptions).
  ushort m_regclass_data_size;

  // Type of packet m_name (sent from server to notify client) 
  uint m_regclass_type;
  
  // Tag for notifying client if was valid registration
  // possible values - VALID_REG, UNKNOWN_CLASS, INVALID_REG_SIZE
  TRRegTag m_tag;
};

#endif
