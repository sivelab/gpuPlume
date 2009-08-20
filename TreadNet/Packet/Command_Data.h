#ifndef COMMAND_DATA_H
#define COMMAND_DATA_H

#include "Packet_Data.h"

// This class is for stairs packet data

class Command_Data : public Packet_Data {
 public:
  friend class Command_Packet;

  // Destructor
  ~Command_Data();
  
  // default constructor
  Command_Data();

  // copy constructor
  Command_Data(const Command_Data&);

  // Assignment operator
  void operator = (const Command_Data&);

  // Accessor functions
  inline TRCommand GetCmd() const {return m_command;}

  // Set functions
  
  // Set Commands
  inline void SetEnableCmd() {m_command = ENABLE_CMD;}
  inline void SetDisableCmd() {m_command = DISABLE_CMD;}

  // Already provides uint GetSize();

 protected:
  // Do not want users to set any command value
  inline void SetCmd(TRCommand cmd) {m_command = cmd;}

  // Inherits variable ushort m_size

 private:
  // command
  TRCommand m_command;
};

#endif
