#include "TRQuaternion.h"

/*********************************************************/
/*  Name : TRQuaternion                                  */
/*                                                       */
/*  Description: Individual element constructor          */
/*                                                       */
/*  Input : double, double, double, double               */
/*  Return Value : None                                  */
/*********************************************************/

TRQuaternion::TRQuaternion(double x, double y, double z, double w) {
  m_qv.X() = x;
  m_qv.Y() = y;
  m_qv.Z() = z;
  m_qw = w;
  Normalize();
}

/*********************************************************/
/*  Name : TRQuaternion                                  */
/*                                                       */
/*  Description: Copy Constructor                        */
/*                                                       */
/*  Input : TRQuaternion&                                */
/*  Return Value : None                                  */
/*********************************************************/

TRQuaternion::TRQuaternion(const TRQuaternion& right) {
  m_qv = right.qV();
  m_qw = right.qW();
}

/*********************************************************/
/*  Name : TRQuaternion                                  */
/*                                                       */
/*  Description: Assigment operator                      */
/*                                                       */
/*  Input : TRQuaternion&                                */
/*  Return Value : None                                  */
/*********************************************************/

void TRQuaternion::operator = (const TRQuaternion& right) {
  m_qv = right.qV();
  m_qw = right.qW();
}

/*********************************************************/
/*  Name : Inverse                                       */
/*                                                       */
/*  Description: Returns inverse of quaternion (this is  */
/*     the same as conjugate and represents a rotation in*/
/*     opposite direction)                               */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : TRQuaternion                          */
/*********************************************************/

TRQuaternion TRQuaternion::Inverse() {
  return TRQuaternion(-m_qv, m_qw);
  
  // Below is if quaternion is not unit
  //  TRQuaternion inv;
  //  double ns = Dot(m_qv,m_qv) + m_qw*m_qw;
  //  inv = Conjugate() / ns;
  //  return inv;
}

/*********************************************************/
/*  Name : RetEquivRotMatrix                             */
/*                                                       */
/*  Description: Returns equivalent rotation matrix of   */
/*     this quaternion.                                  */
/*                                                       */
/*     Assumes unit quaternion                           */
/*                                                       */
/*  Input : Rot Matrix vector references (3 cols)        */
/*  Return Value : none                                  */
/*********************************************************/

void TRQuaternion::RetEquivRotMatrix(TRVector& xcol, TRVector& ycol, TRVector& zcol) {
  double q1 = m_qv.X();
  double q2 = m_qv.Y();
  double q3 = m_qv.Z();

  xcol.X() = 1 - 2*q2*q2 - 2*q3*q3;
  xcol.Y() = 2*(q1*q2 + m_qw*q3);
  xcol.Z() = 2*(q1*q3 - m_qw*q2);
  
  ycol.X() = 2*(q1*q2 - m_qw*q3);
  ycol.Y() = 1 - 2*q1*q1 - 2*q3*q3;
  ycol.Z() = 2*(q2*q3 + m_qw*q1);

  zcol.X() = 2*(q1*q3 + m_qw*q2);
  zcol.Y() = 2*(q2*q3 - m_qw*q1);
  zcol.Z() = 1 - 2*q1*q1 - 2*q2*q2;
}

/*********************************************************/
/*  Name : Normalize                                     */
/*                                                       */
/*  Description: Normalizes quaternion                   */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

void TRQuaternion::Normalize() {
  // get norm and then divide it out
  double nq = sqrt(Dot(m_qv,m_qv) + m_qw*m_qw);
  m_qv /= nq;
  m_qw /= nq;
}


/*********************************************************/
/*  Name : operator *                                    */
/*                                                       */
/*  Description: Combines two quaternions (i.e. combines */
/*    a rotation of quaternion right then left, into a   */
/*    single rotation                                    */
/*                                                       */
/*  Input : TRQuaternion, TRQuaternion                   */
/*  Return Value : TRQuaternion                          */
/*********************************************************/

TRQuaternion operator * (const TRQuaternion& left, const TRQuaternion& right) {
  double qw = left.qW() * right.qW() - Dot(left.qV(), right.qV());
  TRVector qv = left.qW() * right.qV() +
    right.qW() * left.qV() +
    Cross(left.qV(), right.qV());
  return TRQuaternion(qw, qv);
}

/*********************************************************/
/*  Name : Rotate                                        */
/*                                                       */
/*  Description: Rotates a vector based on quaternion    */
/*  Note: q should be formed from angle theta/2 to get   */
/*   a full theta rotation (use CreateQuatFromRotVect w/ */
/*   theta and it will handle this automatically).       */
/*                                                       */
/*  Input : TRQuaternion, TRVector                       */
/*  Return Value : rotated TRVector                      */
/*********************************************************/

TRVector Rotate (const TRQuaternion& q, const TRVector& v) {
  double qw = q.qW();
  TRVector qv = q.qV();
  return (qw*qw - Dot(qv,qv)) * v +
	  2*qw*Cross(qv, v) +
	  2*qv*Dot(qv,v);
}

/*********************************************************/
/*  Name : CreateQuatFromRotVec                          */
/*                                                       */
/*  Description: Creates a Quaternion which represents   */
/*     a rotation theta about vector k                   */
/*                                                       */
/*  Note: Then can use Rotate(q,v) to rotate the vector v*/
/*     about k, theta degrees.                           */
/*                                                       */
/*  Input : angle theta in radians and vector k          */
/*  Return Value : TRQuaternion                          */
/*********************************************************/

TRQuaternion CreateQuatFromRotVec(double theta, const TRVector& k) {
  double qw = cos(theta/2.0);

  // use normalized vector
  TRVector knorm = k;
  knorm.Normalize();

  TRVector qv = knorm * sin(theta/2.0);
  return TRQuaternion(qv, qw);
}

/*********************************************************/
/*  Name : CreateQuatFromRotMat                          */
/*                                                       */
/*  Description: Creates a Quaternion which represents   */
/*     inserted rotation matrix [xcol, ycol, zcol]       */
/*                                                       */
/*  Note: Then can use Rotate(q,v) to rotate the vector v*/
/*     just the same as if using the rotation matrix     */
/*                                                       */
/*  Input : rotation matrix in 3 vector columns          */
/*  Return Value : TRQuaternion                          */
/*********************************************************/

TRQuaternion CreateQuatFromRotMat(const TRVector& x, const TRVector& y, const TRVector& z) {
  // formula from john's class notes chapter 3
  double qw = .5 * sqrt(1 + x.X() + y.Y() + z.Z());
  TRVector qv;

  if (qw > .001) {
    // good enough to use as divisor
    qv.X() = (y.Z() - z.Y()) / (4.0 * qw);
    qv.Y() = (z.X() - x.Z()) / (4.0 * qw);
    qv.Z() = (x.Y() - y.X()) / (4.0 * qw);
  }
  else {
    // must use a different method
    double q1s = .25*(1 + x.X() - y.Y() - z.Z());
    double q2s = .25*(1 - x.X() + y.Y() - z.Z());
    double q3s = .25*(1 - x.X() - y.Y() + z.Z());
    if (q1s > q2s) {
      if (q1s > q3s) {
	// use q1s
	qv.X() = sqrt(q1s);
	qv.Y() = .25*(y.X() + x.Y()) / qv.X();
	qv.Z() = .25*(z.X() + x.Z()) / qv.X();
      }
      else {
	// use q3s
	qv.Z() = sqrt(q3s);
	qv.X() = .25*(z.X() + x.Z()) / qv.Z();
	qv.Y() = .25*(z.Y() + y.Z()) / qv.Z();
      }
    }
    else if (q2s > q3s) {
      // use q2s
      qv.Y() = sqrt(q2s);
      qv.X() = .25*(y.X() + x.Y()) / qv.Y();
      qv.Z() = .25*(z.Y() + y.Z()) / qv.Y();
    }
    else {
      // use q3s
	qv.Z() = sqrt(q3s);
	qv.X() = .25*(z.X() + x.Z()) / qv.Z();
	qv.Y() = .25*(z.Y() + y.Z()) / qv.Z();
    }
  }
  return TRQuaternion(qv, qw);
}

