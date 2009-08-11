#ifndef FSRSIGNALS_DATA_H
#define FSRSIGNALS_DATA_H

#include "Packet_Data.h"

// This class is for FSR Signal packet data

class FSRSignals_Data : public Packet_Data {
 public:
  friend class FSRSignals_Packet;

  // Destructor
  ~FSRSignals_Data();
  
  // default constructor
  FSRSignals_Data();

  // copy constructor
  FSRSignals_Data(const FSRSignals_Data&);

  // Assignment operator
  void operator = (const FSRSignals_Data&);

  // Accessor functions

  // Particular value functions
  inline float GetValueRightHeel() const {return m_right_heel;}
  inline float GetValueRightToe() const {return m_right_toe;}
  inline float GetValueLeftHeel() const {return m_left_heel;}
  inline float GetValueLeftToe() const {return m_left_toe;}

  // Set functions

  // Particular value
  inline void SetValueRightHeel(float rh) {m_right_heel = rh;}
  inline void SetValueRightToe(float rt) {m_right_toe = rt;}
  inline void SetValueLeftHeel(float lh) {m_left_heel = lh;}
  inline void SetValueLeftToe(float lt) {m_left_toe = lt;}

  // Already provides uint GetSize();

  // Inherits variable ushort m_size

private:

  // values for the each individual FSR sensor
  float m_right_heel;
  float m_right_toe;
  float m_left_heel;
  float m_left_toe;
};

#endif
