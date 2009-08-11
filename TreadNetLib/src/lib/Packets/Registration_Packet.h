#ifndef REGISTRATION_PACKET_H
#define REGISTRATION_PACKET_H

#include "Packet.h"
#include "Registration_Data.h"

// A quick definition to simplify the code below
#define REG_DATA dynamic_cast<Registration_Data*>(m_data)

class Registration_Packet : public Packet {
 public:
  // Destructor
  ~Registration_Packet();

  // Default Constructor
  Registration_Packet();

  // Constructor w/ data passed in
  // Note: Data will be copied to a new reg_data class
  Registration_Packet(const Registration_Data&);

  /*************************************************************/
  /* Functions already contained (see Packet.h for details)    */
  /*                                                           */
  /* int GetSize();                                            */
  /* uint GetType();                                           */
  /* BOOL IsRegistered();                                      */
  /*************************************************************/

  // Returns Class Name, for identifying packet
  const char * ClassName() {return Registration_Packet::m_class_name;}

  // Packet_Data functions - these affect the values of m_data, but
  //  does not delete/create m_data.
  int GetCopyOfData(Packet_Data * CopyDat);
  int SetData(Packet_Data * NewDat);

  // Access Functions - just pass along to the data class
  // Set Registration Class Name
  inline int SetRegClassName(const char* name) { return REG_DATA->SetRegClassName(name);  }

  // Set Registration Class Bit Type
  inline void SetRegClassType(uint type) { REG_DATA->SetRegClassType(type); }

  // Set Registration Class Size
  inline void SetRegClassSize(ushort size) { REG_DATA->SetRegClassSize(size); }

  // Set Tag values
  inline void SetValidRegTag() { REG_DATA->SetValidRegTag(); }
  inline void SetUnknownClassRegTag() { REG_DATA->SetUnknownClassRegTag(); }
  inline void SetIncorrectSizeRegTag() { REG_DATA->SetIncorrectSizeRegTag(); }

  // all Get functions
  inline char * GetRegClassName() { return REG_DATA->GetRegClassName(); }
  inline ushort GetRegClassSize() { return REG_DATA->GetRegClassSize(); }
  inline uint GetRegClassType() { return REG_DATA->GetRegClassType(); }
  inline TRRegTag GetRegTag() { return REG_DATA->GetRegTag(); }

 protected:
  // Redefine the set type, because this packet type always
  //  has bit type of 0.
  inline void SetType(uint val) {m_type = 0;};

  // do not want user to be able to set any tag
  inline void SetRegTag(TRRegTag tag) {REG_DATA->SetRegTag(tag);}

  // Fill the message buffer up with the registration data
  virtual int Fill(char *);

  // Read a registration packet from message buffer
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
