#ifndef USERSTATE_PACKET_H
#define USERSTATE_PACKET_H

#include "Packet.h"
#include "UserState_Data.h"

// A quick definition to simplify the code below
#define USERSTATE_DATA dynamic_cast<UserState_Data*>(m_data)

class UserState_Packet : public Packet {
  friend class NewPos_Packet;
  
 public:
  // Destructor
  ~UserState_Packet();

  // Default Constructor
  UserState_Packet();

  // Constructor w/ data passed in
  // Note: Data will be copied to a new UserState_data class
  UserState_Packet(const UserState_Data&);

  /*************************************************************/
  /* Functions already contained (see Packet.h for details)    */
  /*                                                           */
  /* int GetSize();                                            */
  /* uint GetType();                                           */
  /* BOOL IsRegistered();                                      */
  /*************************************************************/

  // Returns Class Name, for identifying packet
  const char * ClassName() {return UserState_Packet::m_class_name;}

  // Packet_Data functions - these affect the values of m_data, but
  //  does not delete/create m_data.
  int GetCopyOfData(Packet_Data * CopyDat);
  int SetData(Packet_Data * NewDat);

  // Access Functions - just pass along to the data class
  //  (both whole and particular)
  inline TRPoint GetPos() const {return USERSTATE_DATA->GetPos();}
  inline TRPoint GetUserOffset() const {return USERSTATE_DATA->GetUserOffset();}
  inline TRPoint GetEyeOffset() const {return USERSTATE_DATA->GetEyeOffset();}
  inline TRVector GetFacing() const {return USERSTATE_DATA->GetFacing();}
  inline TRVector GetVelocity() const {return USERSTATE_DATA->GetVelocity();}
  inline TRStatus GetStatus() const {return USERSTATE_DATA->GetStatus();}
  inline char PosValid() const {return USERSTATE_DATA->PosValid();}
  
  inline double GetPosEasting() const     {return USERSTATE_DATA->GetPosEasting();}
  inline double GetPosNorthing() const    {return USERSTATE_DATA->GetPosNorthing();}
  inline double GetUserEasting() const    {return USERSTATE_DATA->GetUserEasting();}
  inline double GetUserNorthing() const   {return USERSTATE_DATA->GetUserNorthing();}
  inline double GetUserElevation() const  {return USERSTATE_DATA->GetUserElevation();}
  inline double GetEyeEasting() const     {return USERSTATE_DATA->GetEyeEasting();}
  inline double GetEyeNorthing() const    {return USERSTATE_DATA->GetEyeNorthing();}
  inline double GetEyeElevation() const   {return USERSTATE_DATA->GetEyeElevation();}
  inline double GetFacingEasting() const  {return USERSTATE_DATA->GetFacingEasting();}
  inline double GetFacingNorthing() const {return USERSTATE_DATA->GetFacingNorthing();}
  inline double GetVelEasting() const     {return USERSTATE_DATA->GetVelEasting();}
  inline double GetVelNorthing() const    {return USERSTATE_DATA->GetVelNorthing();}
  inline double GetVelVertical() const    {return USERSTATE_DATA->GetVelVertical();}

  // Set functions - again pass along
  inline void SetPos(const TRPoint& pos) {USERSTATE_DATA->SetPos(pos);}
  inline void SetUserOffset(const TRPoint& uoff) {USERSTATE_DATA->SetUserOffset(uoff);}
  inline void SetEyeOffset(const TRPoint& eoff) {USERSTATE_DATA->SetEyeOffset(eoff);}
  inline void SetFacing(const TRVector& facing) {USERSTATE_DATA->SetFacing(facing);}
  inline void SetVelocity(const TRVector& vel) {USERSTATE_DATA->SetVelocity(vel);}
  inline void SetStatus(TRStatus status) {USERSTATE_DATA->SetStatus(status);}
  inline void SetValidFlag() {USERSTATE_DATA->SetValidFlag();}

  inline void SetPosEasting(double e)     {USERSTATE_DATA->SetPosEasting(e);}
  inline void SetPosNorthing(double n)    {USERSTATE_DATA->SetPosNorthing(n);}
  inline void SetUserEasting(double e)    {USERSTATE_DATA->SetUserEasting(e);}
  inline void SetUserNorthing(double n)   {USERSTATE_DATA->SetUserNorthing(n);}
  inline void SetUserElevation(double e)  {USERSTATE_DATA->SetUserElevation(e);}
  inline void SetEyeEasting(double e)     {USERSTATE_DATA->SetEyeEasting(e);}
  inline void SetEyeNorthing(double n)    {USERSTATE_DATA->SetEyeNorthing(n);}
  inline void SetEyeElevation(double e)   {USERSTATE_DATA->SetEyeElevation(e);}
  inline void SetFacingEasting(double e)  {USERSTATE_DATA->SetFacingEasting(e);}
  inline void SetFacingNorthing(double n) {USERSTATE_DATA->SetFacingNorthing(n);}
  inline void SetVelEasting(double e)     {USERSTATE_DATA->SetVelEasting(e);}
  inline void SetVelNorthing(double n)    {USERSTATE_DATA->SetVelNorthing(n);}
  inline void SetVelVertical(double v)    {USERSTATE_DATA->SetVelVertical(v);}

 protected:
  inline void SetPosFlag(char flag) {USERSTATE_DATA->SetPosFlag(flag);}

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
