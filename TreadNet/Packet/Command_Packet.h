#ifndef COMMAND_PACKET_H
#define COMMAND_PACKET_H

#include "Packet.h"
#include "Command_Data.h"

// A quick definition to simplify the code below
#define CMD_DATA dynamic_cast<Command_Data*>(m_data)

class Command_Packet : public Packet {
 public:
  // Destructor
  ~Command_Packet();

  // Default Constructor
  Command_Packet();

  // Constructor w/ data passed in
  // Note: Data will be copied to a new stairs_data class
  Command_Packet(const Command_Data&);

  /*************************************************************/
  /* Functions already contained (see Packet.h for details)    */
  /*                                                           */
  /* int GetSize();                                            */
  /* uint GetType();                                           */
  /* BOOL IsRegistered();                                      */
  /*************************************************************/

  // Returns Class Name, for identifying packet
  const char * ClassName() {return Command_Packet::m_class_name;}

  // Packet_Data functions - these affect the values of m_data, but
  //  does not delete/create m_data.
  int GetCopyOfData(Packet_Data * CopyDat);
  int SetData(Packet_Data * NewDat);

  // Access Functions - just pass along to the data class
  inline TRCommand GetCmd() const { return CMD_DATA->GetCmd(); }

  // Set functions

  // Set Commands
  inline void SetEnableCmd() { CMD_DATA->SetEnableCmd(); }
  inline void SetDisableCmd() { CMD_DATA->SetDisableCmd(); }

 protected:
  // Do not want users to set any command value
  inline void SetCmd(TRCommand cmd) { CMD_DATA->SetCmd(cmd); }

  // Fill the message buffer up with the command data
  virtual int Fill(char *);

  // Read a command packet from message buffer
  virtual int ReadPacket(char *);

  /*************************************************************/
  /* Procedures from Packet                                    */
  /*                                                           */
  /* BOOL WaitingForResponse();                                */
  /* int FillShort(char *, const short&);                      */
  /* ... FillInt, FillDouble, FillFloat, and FillString        */
  /* short GrabShort(char *);                                  */
  /* ..... GrabInt, GrabDouble, GrabFloat                      */
  /*************************************************************/

  // Variable Section

  // Class Name
  static char* m_class_name;

  /*************************************************************/
  /*  Variables already contained in class from Packet         */
  /*                                                           */
  /* Packet_Data* m_data;                                      */
  /* unsigned int m_type;                                      */
  /* BOOL m_wait_response;                                     */
  /* BOOL m_registered;                                        */
  /*************************************************************/

};

#endif
