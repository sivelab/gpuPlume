#ifndef TRQUATERNION_H
#define TRQUATERNION_H

#include "TRVector.h"

/***************************************************************/
/*   Class Name : TRQuaternion                                 */
/*                                                             */
/*   Description - This class represents quaternions used for  */
/*    rotating points and vectors (specifically TRVector &     */
/*    TRPoint).  Thus these types of quaternions are always    */
/*    unit length and are also called Euler Parameters.        */
/***************************************************************/

class TRQuaternion {
 public:
  // default constructor, returns identity quaternion
  TRQuaternion() : m_qv(0), m_qw(1) {}
  
  // copy constructor
  TRQuaternion(const TRQuaternion&);

  // indiv element constructor
  TRQuaternion(double x, double y, double z, double w);

  // two parts constructors
  TRQuaternion(TRVector v, double w) : m_qv(v), m_qw(w) {Normalize();}
  TRQuaternion(double w, TRVector v) : m_qv(v), m_qw(w) {Normalize();}

  // Assignment operator
  void operator = (const TRQuaternion&);

  // Access Functions
  inline double qX() const   { return m_qv.X(); }
  inline double qY() const   { return m_qv.Y(); }
  inline double qZ() const   { return m_qv.Z(); }
  inline double qW() const   { return m_qw;     }
  inline TRVector qV() const { return m_qv;     }

  // Element Assignment Functions
  inline double& qX()   { return m_qv.X(); }
  inline double& qY()   { return m_qv.Y(); }
  inline double& qZ()   { return m_qv.Z(); }
  inline double& qW()   { return m_qw;     }
  inline TRVector& qV() { return m_qv;     }

  // Set Functions
  inline void qX(double new_x)   { m_qv.X() = new_x; }
  inline void qY(double new_y)   { m_qv.Y() = new_y; }
  inline void qZ(double new_z)   { m_qv.Z() = new_z; }
  inline void qW(double new_w)   { m_qw = new_w;     }
  inline void qV(TRVector new_v) { m_qv = new_v;     }

  // Returns Inverse of quaternion - same as conjugate, and represents a
  //  rotation in the opposite direction
  TRQuaternion Inverse();

  // Sets a rotation matrix which equals quaternion
  void RetEquivRotMatrix(TRVector& xcol, TRVector& ycol, TRVector& zcol);

  // Normalizes quaternion - doesn't need to be public, but might
  //  want to ensure normalization after many operations
  void Normalize();  

  // Return size of quaternion data
  inline int GetSize() {return m_qv.GetSize() + sizeof(m_qw);}

 private:

  TRVector m_qv;  // vector portion
  double   m_qw;  // scaler portion
};

// other TRQuaternion operators
TRQuaternion operator * (const TRQuaternion& left, const TRQuaternion& right);

// Rotates a vector v using quaternion q (q should be formed from angle 
//  theta/2 to get a full theta rotation).
TRVector Rotate (const TRQuaternion& q, const TRVector& v);

// Create quaternion to represent rotation theta (radians) about vector k
TRQuaternion CreateQuatFromRotVec(double theta, const TRVector& k);

// Create quaternion from a Rotation matrix
TRQuaternion CreateQuatFromRotMat(const TRVector& x, const TRVector& y, const TRVector& z);

#endif
