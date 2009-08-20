#ifndef USERSTATE_DATA_H
#define USERSTATE_DATA_H

#include "Packet_Data.h"

// This class is for UserState packet data

class UserState_Data : public Packet_Data {
 public:
  friend class UserState_Packet;

  // Destructor
  ~UserState_Data();
  
  // default constructor
  UserState_Data();

  // copy constructor
  UserState_Data(const UserState_Data&);

  // Assignment operator
  void operator = (const UserState_Data&);

  // Accessor functions
  
  // Whole values
  inline TRPoint GetPos() const {return m_pos;}
  inline TRPoint GetUserOffset() const {return m_user_offset;}
  inline TRPoint GetEyeOffset() const {return m_eye_offset;}
  inline TRVector GetFacing() const {return m_facing;}
  inline TRVector GetVelocity() const {return m_velocity;}
  inline TRStatus GetStatus() const {return m_status;}
  inline char PosValid() const {return m_valid;}

  // Local Whole Values (i.e. offsets in Treadport Local Coordinate System)
  TRPoint GetLocalUserOffset() const;
  TRPoint GetLocalEyeOffset() const;

  // Particular value functions
  inline double GetPosEasting() const     {return m_pos.X();}
  inline double GetPosNorthing() const    {return m_pos.Y();}
  inline double GetUserEasting() const    {return m_user_offset.X();}
  inline double GetUserNorthing() const   {return m_user_offset.Y();}
  inline double GetUserElevation() const  {return m_user_offset.Z();}
  inline double GetEyeEasting() const     {return m_eye_offset.X();}
  inline double GetEyeNorthing() const    {return m_eye_offset.Y();}
  inline double GetEyeElevation() const   {return m_eye_offset.Z();}
  inline double GetFacingEasting() const  {return m_facing.X();}
  inline double GetFacingNorthing() const {return m_facing.Y();}
  inline double GetVelEasting() const     {return m_velocity.X();}
  inline double GetVelNorthing() const    {return m_velocity.Y();}
  inline double GetVelVertical() const    {return m_velocity.Z();}

  // Set functions

  // Whole values
  inline void SetPos(const TRPoint& pos) {m_pos = pos;}
  inline void SetUserOffset(const TRPoint& uoff) {m_user_offset = uoff;}
  inline void SetEyeOffset(const TRPoint& eoff) {m_eye_offset = eoff;}
  inline void SetFacing(const TRVector& facing) {m_facing = facing;}
  inline void SetVelocity(const TRVector& vel) {m_velocity = vel;}
  inline void SetStatus(TRStatus status) {m_status = status;}
  inline void SetValidFlag() {m_valid = 1;}

  // Local Set Whole values (Make sure you have set facing correctly before
  //  using this function to ensure the correct conversion);
  void SetLocalUserOffset(const TRPoint& uoff);
  void SetLocalEyeOffset(const TRPoint& eoff);

  // Particular value
  inline void SetPosEasting(double e)     {m_pos.X(e);}
  inline void SetPosNorthing(double n)    {m_pos.Y(n);}
  inline void SetUserEasting(double e)    {m_user_offset.X(e);}
  inline void SetUserNorthing(double n)   {m_user_offset.Y(n);}
  inline void SetUserElevation(double e)  {m_user_offset.Z(e);}
  inline void SetEyeEasting(double e)     {m_eye_offset.X(e);}
  inline void SetEyeNorthing(double n)    {m_eye_offset.Y(n);}
  inline void SetEyeElevation(double e)   {m_eye_offset.Z(e);}
  inline void SetFacingEasting(double e)  {m_facing.X(e);}
  inline void SetFacingNorthing(double n) {m_facing.Y(n);}
  inline void SetVelEasting(double e)     {m_velocity.X(e);}
  inline void SetVelNorthing(double n)    {m_velocity.Y(n);}
  inline void SetVelVertical(double v)    {m_velocity.Z(v);}

  // Already provides uint GetSize();

 protected:
  // do not want user to set any value
  inline void SetPosFlag(char flag) {m_valid = flag;}

  // Inherits variable ushort m_size

  // Protected functions for my simplicity, convert between local treadport
  //  and virtual world coordinate systems using m_facing
  void ConvertLocalToVW(TRVector&) const;
  void ConvertVWToLocal(TRVector&) const;

 private:
  // For all Positions and Vectors the following holds true (X - East, Y - North)

  // Treadport center position in virtual world (Z - nothing)
  TRPoint  m_pos;

  // User Offset from treadport surface center in virtual world coordinates (Z - Up)
  //  (user is defined as tether attachment point - or end of tether, or user's back)
  TRPoint  m_user_offset;

  // Eye Offset from treadport surface center in virtual world coordinates (Z - Up)
  //  (eye is up to your eye - where torso height/width is defined in the Treadport app)
  TRPoint  m_eye_offset;

  // Facing in virtual world (Z - nothing)
  TRVector m_facing;

  // Treadport velocity in virtual world (Z - Up)
  TRVector m_velocity;

  // Treadport status (like active, enabled, ramping up/down, etc.)
  TRStatus m_status;

  // flag for determining if position is valid.  Essentially it becomes
  //  valid once a new Position is received by the server.
  char m_valid;
};

#endif
