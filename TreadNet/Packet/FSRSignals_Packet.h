#ifndef FSRSIGNALS_PACKET_H
#define FSRSIGNALS_PACKET_H

#include "Packet.h"
#include "FSRSignals_Data.h"

// A quick definition to simplify the code below
#define FSRSIGNALS_DATA dynamic_cast<FSRSignals_Data*>(m_data)

class FSRSignals_Packet : public Packet {
 public:
  // Destructor
  ~FSRSignals_Packet();

  // Default Constructor
  FSRSignals_Packet();

  // Constructor w/ data passed in
  // Note: Data will be copied to a new FSRSignals_data class
  FSRSignals_Packet(const FSRSignals_Data&);

  /*************************************************************/
  /* Functions already contained (see Packet.h for details)    */
  /*                                                           */
  /* int GetSize();                                            */
  /* uint GetType();                                           */
  /* BOOL IsRegistered();                                      */
  /*************************************************************/

  // Returns Class Name, for identifying packet
  const char * ClassName() {return FSRSignals_Packet::m_class_name;}

  // Packet_Data functions - these affect the values of m_data, but
  //  does not delete/create m_data.
  int GetCopyOfData(Packet_Data * CopyDat);
  int SetData(Packet_Data * NewDat);

  // Access Functions - just pass along to the data class
  //  (both whole and particular)
  inline float GetValueRightHeel() const {return FSRSIGNALS_DATA->GetValueRightHeel();}
  inline float GetValueRightToe() const {return FSRSIGNALS_DATA->GetValueRightToe();}
  inline float GetValueLeftHeel() const {return FSRSIGNALS_DATA->GetValueLeftHeel();}
  inline float GetValueLeftToe() const {return FSRSIGNALS_DATA->GetValueLeftToe();}
  
  // Set functions - again pass along
  inline void SetValueRightHeel(const float rh) {FSRSIGNALS_DATA->SetValueRightHeel(rh);}
  inline void SetValueRightToe(const float rt) {FSRSIGNALS_DATA->SetValueRightToe(rt);}
  inline void SetValueLeftHeel(const float lh) {FSRSIGNALS_DATA->SetValueLeftHeel(lh);}
  inline void SetValueLeftToe(const float lt) {FSRSIGNALS_DATA->SetValueLeftToe(lt);}

protected:

  // Fill the message buffer up with the fsr data
  virtual int Fill(char *);

  // Read a fsr packet from message buffer
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
