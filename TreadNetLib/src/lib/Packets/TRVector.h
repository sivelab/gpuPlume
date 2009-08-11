#ifndef TRVECTOR_H
#define TRVECTOR_H

#include "Treadport_Types.h"

#define TRPoint TRVector

class TRVector {
 public:
  // default constructor
  TRVector() { m_x = m_y = m_z = 0; }
  
  // copy constructor
  TRVector(const TRVector&);

  // indiv element constructor
  TRVector(double x, double y, double z) {
    m_x = x;
    m_y = y;
    m_z = z;
  }

  // single value constructor
  TRVector(double val) { m_x = m_y = m_z = val; }

  // Assignment operator
  void operator = (const TRVector&);

  // Access Functions
  inline double X() const { return m_x; }
  inline double Y() const { return m_y; }
  inline double Z() const { return m_z; }

  // Element Assignment Functions
  inline double& X() { return m_x; }
  inline double& Y() { return m_y; }
  inline double& Z() { return m_z; }

  // Set Functions
  inline void X(double new_x) { m_x = new_x; }
  inline void Y(double new_y) { m_y = new_y; }
  inline void Z(double new_z) { m_z = new_z; }

  // Return size of vector data
  inline int GetSize() {return 3*sizeof(m_x);}

  // Returns norm & squared norm of vector
  inline double SquaredNorm() const { return m_x*m_x + m_y*m_y + m_z*m_z; }
  inline double Norm() const { return sqrt(SquaredNorm()); }

  // Normalizes this vector
  void Normalize();

  // Other operators
  void operator += (const TRVector&);
  void operator -= (const TRVector&);
  void operator *= (double);
  void operator /= (double);
  TRVector operator - ();
  BOOL operator==(const TRVector&);
  BOOL operator!=(const TRVector&);

 private:
  double m_x;
  double m_y;
  double m_z;
};

// other TRVector operators
TRVector operator + (const TRVector& left, const TRVector& right);
TRVector operator - (const TRVector& left, const TRVector& right);
TRVector operator * (const TRVector& left, double mult);
TRVector operator * (double mult, const TRVector& right);
TRVector operator / (const TRVector& left, double div);

TRVector Cross (const TRVector& left, const TRVector& right);
double   Dot (const TRVector& left, const TRVector& right);

#endif
