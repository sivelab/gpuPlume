#ifndef VIRTUALSTATE_DATA_H
#define VIRTUALSTATE_DATA_H

#include "Packet_Data.h"

// This class is for Virtual State packet data

class VirtualState_Data : public Packet_Data {
 public:
  // Destructor
  ~VirtualState_Data();
  
  // default constructor
  VirtualState_Data();

  // copy constructor
  VirtualState_Data(const VirtualState_Data&);

  // Assignment operator
  void operator = (const VirtualState_Data&);

  // Accessor functions
  inline TRVector GetSurfaceNormal() const {return m_surface_normal;}
  inline TRTerrain GetTerrain() const {return m_terrain;}
  inline TRVector GetVirtualForce() const {return m_virtual_force;}
  inline float GetMovementScale() const {return m_movement_scale;}
  inline float GetUserScale() const {return m_user_scale;}

  // Set functions
  inline void SetSurfaceNormal(const TRVector& sn) {m_surface_normal = sn;}
  inline void SetTerrain(TRTerrain ter) {m_terrain = ter;}
  inline void SetVirtualForce(const TRVector& vf) {m_virtual_force = vf;}
  inline void SetMovementScale(float ms) {m_movement_scale = ms;}
  inline void SetUserScale(float us) {m_user_scale = us;}

  // Already provides uint GetSize();

 protected:
  // Inherits variable ushort m_size

 private:
  // For all Vectors the following holds true (X - East, Y - North, Z - Up)

  // Surface Normal for virtual world
  TRVector m_surface_normal;

  // Terrain value
  TRTerrain m_terrain;

  // Virtual force to apply
  TRVector m_virtual_force;

  // Movement scale
  float m_movement_scale;

  // User scale
  float m_user_scale;
};

#endif
