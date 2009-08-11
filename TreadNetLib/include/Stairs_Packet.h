#ifndef STAIRS_PACKET_H
#define STAIRS_PACKET_H

#include "Packet.h"
#include "Stairs_Data.h"

// A quick definition to simplify the code below
#define STAIR_DATA dynamic_cast<Stairs_Data*>(m_data)

class Stairs_Packet : public Packet {
 public:
  // Destructor
  ~Stairs_Packet();

  // Default Constructor
  Stairs_Packet();

  // Constructor w/ data passed in
  // Note: Data will be copied to a new stairs_data class
  Stairs_Packet(const Stairs_Data&);

  /*************************************************************/
  /* Functions already contained (see Packet.h for details)    */
  /*                                                           */
  /* int GetSize();                                            */
  /* uint GetType();                                           */
  /* BOOL IsRegistered();                                      */
  /*************************************************************/

  // Returns Class Name, for identifying packet
  const char * ClassName() {return Stairs_Packet::m_class_name;}

  // Packet_Data functions - these affect the values of m_data, but
  //  does not delete/create m_data.
  int GetCopyOfData(Packet_Data * CopyDat);
  int SetData(Packet_Data * NewDat);

  // Access Functions - just pass along to the data class
  inline TRStairs_Cmd GetCmd() const { return STAIR_DATA->GetCmd(); }
  inline TRStairs_Dir GetStairsDir() const { return STAIR_DATA->GetStairsDir(); }

  // Set functions

  // Set Commands
  inline void SetSwitchStairsCmd() { STAIR_DATA->SetSwitchStairsCmd(); }
  inline void SetIncConstFactCmd() { STAIR_DATA->SetIncConstFactCmd(); }
  inline void SetDecConstFactCmd() { STAIR_DATA->SetDecConstFactCmd(); }
  inline void SetNoStairsCmd()     { STAIR_DATA->SetNoStairsCmd();     }

  // Set Stairs Direction
  inline void SetUpStairsDir()         { STAIR_DATA->SetUpStairsDir();     }
  inline void SetDownStairsDir()       { STAIR_DATA->SetDownStairsDir();   }
  inline void SetIgnoreStairsDir()     { STAIR_DATA->SetIgnoreStairsDir(); }

 protected:
  // Do not want user to be able to set any value for these
  inline void SetCmd(TRStairs_Cmd cmd)       { STAIR_DATA->SetCmd(cmd);       }
  inline void SetStairsDir(TRStairs_Dir dir) { STAIR_DATA->SetStairsDir(dir); }

  // Fill the message buffer up with the stair data
  virtual int Fill(char *);

  // Read a stair packet from message buffer
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
