#ifndef COLLISION_DATA_H
#define COLLISION_DATA_H

#include "Packet_Data.h"

enum TRCollFlag {
  NO_COLLISION,
  COLLISION
};

// This class is for Collision packet data

class Collision_Data : public Packet_Data {
 public:
  friend class Collision_Packet;

  // Destructor
  ~Collision_Data();
  
  // default constructor
  Collision_Data();

  // copy constructor
  Collision_Data(const Collision_Data&);

  // Assignment operator
  void operator = (const Collision_Data&);

  // Accessor functions
  inline TRCollFlag GetCollisionFlag() const {return m_collision_flag;}
  inline TRVector GetCollisionNormal() const {return m_collision_normal;}
  inline TRPoint GetCollisionPos() const {return m_collision_pos;}

  // Set functions
  inline void SetNoCollision() {m_collision_flag = NO_COLLISION;}
  inline void SetCollision() {m_collision_flag = COLLISION;}

  inline void SetCollisionNormal(const TRVector& cn) {m_collision_normal = cn;}
  inline void SetCollisionPos(const TRVector& cp) {m_collision_pos = cp;}

  // Already provides uint GetSize();

 protected:
  // Do not want users to set any flag value
  inline void SetCollisionFlag(TRCollFlag f) {m_collision_flag = f;}

  // Inherits variable ushort m_size

 private:
  // For all Vectors/Points the following holds true (X - East, Y - North, Z - Up)

  // Collision flag
  TRCollFlag m_collision_flag;

  // Collision surface normal
  TRVector m_collision_normal;

  // Collision surface position
  TRPoint m_collision_pos;
};

#endif
