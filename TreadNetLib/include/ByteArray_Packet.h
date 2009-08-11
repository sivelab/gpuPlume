#ifndef BYTEARRAY_PACKET_H
#define BYTEARRAY_PACKET_H

#include "Packet.h"
#include "ByteArray_Data.h"

/***************************************************************/
/* One note of warning, this packet's size can change each call*/
/*  This has the potential to confuse NetController when it    */
/*  checks for buffer overflow in GetPacketFromMsg.  Namely if */
/*  you send a huge array, then send a smaller one along w/    */
/*  alot of other data the function might think overflow even  */
/*  though it's not since the array shrunk from last instance. */
/*  This is only when reading an array, sending has no problems*/
/* I couldn't find a easy, non-hacked way to remove this, and  */
/*  it's a rare condition (large buffer, plus doubt change     */
/*  array size often), so left it alone.                       */
/***************************************************************/ 

// A quick definition to simplify the code below
#define BA_DATA dynamic_cast<ByteArray_Data*>(m_data)

class ByteArray_Packet : public Packet {
 public:
  // Destructor
  ~ByteArray_Packet();

  // Default Constructor
  ByteArray_Packet();

  // Constructor w/ data passed in
  // Note: Data will be copied to a new reg_data class
  ByteArray_Packet(const ByteArray_Data&);

  /*************************************************************/
  /* Functions already contained (see Packet.h for details)    */
  /*                                                           */
  /* int GetSize();                                            */
  /* uint GetType();                                           */
  /* BOOL IsRegistered();                                      */
  /*************************************************************/

  // Returns Class Name, for identifying packet
  const char * ClassName() {return ByteArray_Packet::m_class_name;}

  // Packet_Data functions - these affect the values of m_data, but
  //  does not delete/create m_data.
  int GetCopyOfData(Packet_Data * CopyDat);
  int SetData(Packet_Data * NewDat);

  // Access Functions - just pass along to the data class
  // Set ByteArray element
  inline int SetByte(int index, char value) { return BA_DATA->SetByte(index, value);  }

  // Set/get last element index
  inline int SetLastElementInd(int val) { return BA_DATA->SetLastElementInd(val); }
  inline int GetLastElementInd() { return BA_DATA->GetLastElementInd(); }

  // Get Array ptr (both versions)
  inline char* GetArrayPtr(int index = 0) { return BA_DATA->GetArrayPtr(index); }
  inline const char* GetArrayPtr(int index = 0) const { return BA_DATA->GetArrayPtr(index); }

 protected:
  // Fill the message buffer up with the bytearray data
  virtual int Fill(char *);

  // Read a bytearray packet from message buffer
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
