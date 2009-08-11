#ifndef NEWPOS_PACKET_H
#define NEWPOS_PACKET_H

#include "Packet.h"
#include "NewPos_Data.h"
#include "VirtualState_Packet.h"
#include "UserState_Packet.h"

// A quick definition to simplify the code below
#define NEWPOS_DATA dynamic_cast<NewPos_Data*>(m_data)

class NewPos_Packet : public Packet {
 public:
  // Destructor
  ~NewPos_Packet();

  // Default Constructor
  NewPos_Packet();

  // Constructor w/ data passed in
  // Note: Data will be copied to a new NewPos_data class
  NewPos_Packet(const NewPos_Data&);

  // Constructor with individual state data passed in
  NewPos_Packet(const UserState_Data&, const VirtualState_Data&);

  /*************************************************************/
  /* Functions already contained (see Packet.h for details)    */
  /*                                                           */
  /* int GetSize();                                            */
  /* uint GetType();                                           */
  /* BOOL IsRegistered();                                      */
  /*************************************************************/

  // Returns Class Name, for identifying packet
  const char * ClassName() {return NewPos_Packet::m_class_name;}

  // Packet_Data functions - these affect the values of m_data, but
  //  does not delete/create m_data.
  int GetCopyOfData(Packet_Data * CopyDat);
  int SetData(Packet_Data * NewDat);

  // Access Functions - just pass along to the data class
  inline UserState_Data GetUserState() const {return NEWPOS_DATA->GetUserState();}
  inline VirtualState_Data GetVirtualState() const {return NEWPOS_DATA->GetVirtualState();}

  // Set functions - again pass along
  inline void SetUserState(const UserState_Data& usd) {NEWPOS_DATA->SetUserState(usd);}
  inline void SetVirtualState(const VirtualState_Data& vsd) {NEWPOS_DATA->SetVirtualState(vsd);}

 protected:
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
