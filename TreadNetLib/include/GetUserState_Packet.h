#ifndef GETUSERSTATE_PACKET_H
#define GETUSERSTATE_PACKET_H

#include "Packet.h"
#include "GetUserState_Data.h"

// A quick definition to simplify the code below
#define GETUSTATE_DATA dynamic_cast<GetUserState_Data*>(m_data)

class GetUserState_Packet : public Packet {
 public:
  // Destructor
  ~GetUserState_Packet();

  // Default Constructor
  GetUserState_Packet();

  // Constructor w/ data passed in
  // Note: Data will be copied to a new GetUserState_data class
  GetUserState_Packet(const GetUserState_Data&);

  /*************************************************************/
  /* Functions already contained (see Packet.h for details)    */
  /*                                                           */
  /* int GetSize();                                            */
  /* uint GetType();                                           */
  /* BOOL IsRegistered();                                      */
  /*************************************************************/

  // Returns Class Name, for identifying packet
  const char * ClassName() {return GetUserState_Packet::m_class_name;}

  // Packet_Data functions - there is no data to copy, so just return 1
  int GetCopyOfData(Packet_Data * CopyDat) { return 1; }
  int SetData(Packet_Data * NewDat) { return 1; }

 protected:
  // Fill the message buffer up with the getuserstate data (none)
  virtual int Fill(char *) { return 0; }

  // Read a getuserstate packet from message buffer (return success)
  virtual int ReadPacket(char *) { return 1; }

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
