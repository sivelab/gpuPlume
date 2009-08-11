#ifndef TRACKER_DATA_H
#define TRACKER_DATA_H

#include "Packet_Data.h"

// This class is for Tracker packet data

class Tracker_Data : public Packet_Data {
 public:
  friend class Tracker_Packet;

  // Destructor
  ~Tracker_Data();
  
  // default constructor
  Tracker_Data();

  // copy constructor
  Tracker_Data(const Tracker_Data&);

  // Assignment operator
  void operator = (const Tracker_Data&);

  // Accessor functions
  
  // Whole values
  inline TRPoint GetPos() const {return m_pos;}
  inline TRQuaternion GetOrientation() const {return m_orient;}

  // Particular value functions
  inline double   GetPosX() const         { return m_pos.X();}
  inline double   GetPosY() const         { return m_pos.Y();}
  inline double   GetPosZ() const         { return m_pos.Z();}
  inline double   GetOrientRot() const    { return m_orient.qW(); }
  inline TRVector GetOrientVector() const { return m_orient.qV(); }

  // Set functions

  // Whole values
  inline void SetPos(const TRPoint& pos) {m_pos = pos;}
  inline void SetOrientation(const TRQuaternion& orient) {m_orient = orient;}

  // Particular value
  inline void SetPosX(double x)                  {m_pos.X(x);}
  inline void SetPosY(double y)                  {m_pos.Y(y);}
  inline void SetPosZ(double z)                  {m_pos.Z(z);}
  inline void SetOrientRot(double r)             {m_orient.qW(r);}
  inline void SetOrientVector(const TRVector& v) {m_orient.qV(v);}

  // Already provides uint GetSize();

 protected:
  // Inherits variable ushort m_size

 private:
  // Tracker center position in virtual world
  TRPoint  m_pos;

  // Orientation of tracker in actual (tracker) world
  TRQuaternion m_orient;
};

#endif
