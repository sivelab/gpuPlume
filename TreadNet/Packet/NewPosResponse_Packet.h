#ifndef NEWPOSRESPONSE_PACKET_H
#define NEWPOSRESPONSE_PACKET_H

#include "Packet.h"
#include "NewPosResponse_Data.h"

// A quick definition to simplify the code below
#define NEWPOSRESP_DATA dynamic_cast<NewPosResponse_Data*>(m_data)

class NewPosResponse_Packet : public Packet {
 public:
  // Destructor
  ~NewPosResponse_Packet();

  // Default Constructor
  NewPosResponse_Packet();

  // Constructor w/ data passed in
  // Note: Data will be copied to a new NewPosResponse_data class
  NewPosResponse_Packet(const NewPosResponse_Data&);

  /*************************************************************/
  /* Functions already contained (see Packet.h for details)    */
  /*                                                           */
  /* int GetSize();                                            */
  /* uint GetType();                                           */
  /* BOOL IsRegistered();                                      */
  /*************************************************************/

  // Returns Class Name, for identifying packet
  const char * ClassName() {return NewPosResponse_Packet::m_class_name;}

  // Packet_Data functions - these affect the values of m_data, but
  //  does not delete/create m_data.
  int GetCopyOfData(Packet_Data * CopyDat);
  int SetData(Packet_Data * NewDat);

  // Access Functions - just pass along to the data class
  inline TRNewPosFlag GetFlag() const { return NEWPOSRESP_DATA->GetFlag(); }

  // Set functions
  inline void SetSuccessFlag() { NEWPOSRESP_DATA->SetSuccessFlag(); }
  inline void SetFailFlag() { NEWPOSRESP_DATA->SetFailFlag(); }

 protected:
  // Do not want users setting any flag value
  inline void SetFlag(TRNewPosFlag flag) { NEWPOSRESP_DATA->SetFlag(flag); }

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
