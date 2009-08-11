#include "TRVector.h"

/*********************************************************/
/*  Name : TRVector                                      */
/*                                                       */
/*  Description: Copy Constructor                        */
/*                                                       */
/*  Input : TRVector&                                    */
/*  Return Value : None                                  */
/*********************************************************/

TRVector::TRVector(const TRVector& right) {
  m_x = right.X();
  m_y = right.Y();
  m_z = right.Z();
}

/*********************************************************/
/*  Name : TRVector                                      */
/*                                                       */
/*  Description: Assigment operator                      */
/*                                                       */
/*  Input : TRVector&                                    */
/*  Return Value : None                                  */
/*********************************************************/

void TRVector::operator = (const TRVector& right) {
  m_x = right.X();
  m_y = right.Y();
  m_z = right.Z();
}

/*********************************************************/
/*  Name : TRVector                                      */
/*                                                       */
/*  Description: Normalize vector                        */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

void TRVector::Normalize() {
  *this /= Norm();
}

/*********************************************************/
/*  Name : Equality operators ==, !=                     */
/*                                                       */
/*  Description: equality operators                      */
/*                                                       */
/*  Input : TRVector&                                    */
/*  Return Value : Bool                                  */
/*********************************************************/

BOOL TRVector::operator == (const TRVector& right) {
  if (m_x != right.X()) return FALSE;
  if (m_y != right.Y()) return FALSE;
  if (m_z != right.Z()) return FALSE;
  return TRUE;
}

BOOL TRVector::operator != (const TRVector& right) {
  return !((*this) == right);
}

/*********************************************************/
/*  Name : Member operators +=, -=, *=, /=, -            */
/*                                                       */
/*  Description: Standard operators                      */
/*                                                       */
/*  Input : TRVector&, double, or nothing                */
/*  Return Value : Changes the internal values, except /,- */
/*********************************************************/

void TRVector::operator += (const TRVector& right) {
  // += operator
  m_x += right.X();
  m_y += right.Y();
  m_z += right.Z();
}

void TRVector::operator -= (const TRVector& right) {
  // -= operator
  m_x -= right.X();
  m_y -= right.Y();
  m_z -= right.Z();
}

void TRVector::operator *= (double right) {
  // *= operator, only for a single value on right, not vector
  m_x *= right;
  m_y *= right;
  m_z *= right;
}

void TRVector::operator /= (double right) {
  // /= operator, only for a single value on right, not vector
  m_x /= right;
  m_y /= right;
  m_z /= right;
}

TRVector TRVector::operator - () {
  // - operator, negative operator
  return TRVector(-m_x, -m_y, -m_z);
}

/*********************************************************/
/*  Name : Binary operators +, -, *, /                   */
/*                                                       */
/*  Description: Standard operators                      */
/*                                                       */
/*  Input : TRVector& or double                          */
/*  Return Value : another TRVector                      */
/*********************************************************/

TRVector operator + (const TRVector& left, const TRVector& right) {
  // addition operator
  return TRVector(left.X() + right.X(), 
		  left.Y() + right.Y(),
		  left.Z() + right.Z());
}

TRVector operator - (const TRVector& left, const TRVector& right) {
  // substraction operator
  return TRVector(left.X() - right.X(), 
		  left.Y() - right.Y(),
		  left.Z() - right.Z());
}

TRVector operator * (double mult, const TRVector& right) {
  // multiplication operator, only for a double input
  return TRVector(mult*right.X(), mult*right.Y(), mult*right.Z());
}

TRVector operator * (const TRVector& left, double mult) {
  return TRVector(mult*left.X(), mult*left.Y(), mult*left.Z());
}

TRVector operator / (const TRVector& left, double div) {
  // / operator, only for a single value on right, not left
  //  i.e. TRVector / 4, not 4 / TRVector
  return TRVector(left.X() / div, left.Y() / div, left.Z() / div);
}


/*********************************************************/
/*  Name : Cross                                         */
/*                                                       */
/*  Description: Cross Product                           */
/*                                                       */
/*  Input : 2 TRVectors                                  */
/*  Return Value : another TRVector                      */
/*********************************************************/

TRVector Cross (const TRVector& left, const TRVector& right) {
  TRVector result;
  // cross product goes as
  //  | lx ly lz |
  //  | rx ry rz | == (ly*rz - lz*ry)i - (lx*rz - lz*rx)j + (lx*ry - ly*rx)k
  //  | i  j  k  |
  result.X() = left.Y()*right.Z() - left.Z()*right.Y();
  result.Y() = left.Z()*right.X() - left.X()*right.Z();
  result.Z() = left.X()*right.Y() - left.Y()*right.X();
  return result;
}

/*********************************************************/
/*  Name : Dot                                           */
/*                                                       */
/*  Description: Dot Product                             */
/*                                                       */
/*  Input : 2 TRVectors                                  */
/*  Return Value : double                                */
/*********************************************************/

double Dot (const TRVector& left, const TRVector& right) {
  return left.X() * right.X() +
    left.Y() * right.Y() +
    left.Z() * right.Z();
}
