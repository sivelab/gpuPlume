#ifndef SERVERREQUEST_PACKET_H
#define SERVERREQUEST_PACKET_H

#include "Packet.h"
#include "ServerRequest_Data.h"

// A quick definition to simplify the code below
#define SERVERREQUEST_DATA dynamic_cast<ServerRequest_Data*>(m_data)

class ServerRequest_Packet : public Packet {
 public:
  // Destructor
  ~ServerRequest_Packet();

  // Default Constructor
  ServerRequest_Packet();

  // Constructor w/ data passed in
  // Note: Data will be copied to a new ServerRequest_data class
  ServerRequest_Packet(const ServerRequest_Data&);

  /*************************************************************/
  /* Functions already contained (see Packet.h for details)    */
  /*                                                           */
  /* int GetSize();                                            */
  /* uint GetType();                                           */
  /* BOOL IsRegistered();                                      */
  /*************************************************************/

  // Returns Class Name, for identifying packet
  const char * ClassName() {return ServerRequest_Packet::m_class_name;}

  // Packet_Data functions
  int GetCopyOfData(Packet_Data * CopyDat);
  int SetData(Packet_Data * NewDat);

  // Accessor functions to the data - just pass on
  // Allow user to remove and add packets to be requested
  inline int AddPacketRequest(Packet * req_pkt)    { return SERVERREQUEST_DATA->AddPacketRequest(req_pkt); }
  inline int RemovePacketRequest(Packet * req_pkt) { return SERVERREQUEST_DATA->RemovePacketRequest(req_pkt); }
  
  // if want to reset field
  inline void ZeroRequestField() { SERVERREQUEST_DATA->ZeroRequestField(); }
  
  // ask if packet is in request field (1 - yes, 0 - no)
  inline int IsPacketRequested(Packet * req_pkt)  { return SERVERREQUEST_DATA->IsPacketRequested(req_pkt); }


 protected:
  // Fill the message buffer up with the ServerRequest data
  virtual int Fill(char *);

  // Read a ServerRequest packet from message buffer
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
