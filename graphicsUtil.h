#ifndef GRAPHICS_UTIL_H
#define GRAPHICS_UTIL_H 1

#include <math.h>
#include <string>

#include <assert.h>
#include <iostream>
#include <math.h>

#include "ansiCPP.H"
#include "approxMath.H"

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <stdlib.h>
#include <GL/glut.h>
// #include <GL/gl.h>
#endif

// class QuatR4;			// Defined in quatR4.H

using namespace std;

class Texture 
{
public:
  Texture();
  Texture( std::string texture_filename );

  void load();
  void bind() const;

  static void enable() { 
    #ifndef __NO_OGL__
    glEnable( GL_TEXTURE_2D ); 
    #endif
  }
  static void disable() { 
    #ifndef __NO_OGL__
    glDisable( GL_TEXTURE_2D ); 
    #endif
  }

private:
  std::string _tex_filename;
  #ifndef __NO_OGL__
  GLuint _tex_id;
  #else
  uint _tex_id;
  #endif
};


/* Camera class for use in OpenGL.  Based on movements common in
   various games. */
class OpenGLCamera {
  
};

/* 
class VertexProgram
{
  VertexProgram();
};

class PixelShader
{
  PixelShader();
}; */


// class QuatR4;			// Defined in quatR4.H

class Vector2 {
  friend ostream& operator<<(ostream& os, const Vector2& v);
  friend istream& operator>>(istream& is, Vector2& v); 
  double d[2];
 public:
  Vector2() {}
  Vector2(const double a, const double b);
  virtual const double& operator[](const int i) const  // inspection version
  { assert(i >= 0 && i < 2); return d[i]; }
  virtual double& operator[](const int i)           // lvalue version 
  { assert(i >= 0 && i < 2); return d[i]; }
  virtual void set(const double a, const double b);
  virtual void set(const double v[2]);
  double norm(void) const;                   // compute L_2 norm
  double normL1(void) const;	             // compute L_1 norm 
  double normLInf(void) const;	             // compute L_\infty norm
  void normalize(void);                      // destructively normalize self
  double dot(const Vector2&) const;          // dot product
  double cross(const Vector2&) const;       // cross product
  Vector2& operator+=(const Vector2&);
  Vector2& operator-=(const Vector2&);
  Vector2& operator*=(const double c);
  Vector2& operator/=(const double c);
  bool operator==(const Vector2& rhs) const;
  bool operator!=(const Vector2& rhs) const;
  bool operator<(const Vector2& rhs) const;   // checks if v1[i] < v2[i] for all i, true if all are true
  bool operator>(const Vector2& rhs) const;   // checks if v1[i] > v2[i] for all i, true if all are true
};

inline const Vector2 operator+(const Vector2& lhs, const Vector2& rhs);
inline const Vector2 operator-(const Vector2& lhs, const Vector2& rhs);
inline const Vector2 operator*(const double, const Vector2&);// scalar multiply
inline const Vector2 operator*(const Vector2&, const double);
inline const Vector2 operator/(const Vector2& v, const double c);// vect/scalar

double anglebetween(const Vector2&, const Vector2&);    // radians
double signedanglebetween(const Vector2 a, const Vector2 b, const Vector2 c); 
                                             // from a to b, wrt c in radians 

class Vector3 {
  friend ostream& operator<<(ostream& os, const Vector3& v);
  friend istream& operator>>(istream& is, Vector3& v); 
  double d[3];
 public:
  Vector3() {}
  Vector3(const double a, const double b, const double c);
  virtual const double& operator[](const int i) const  // inspection version
  { assert(i >= 0 && i < 3); return d[i]; }
  virtual double& operator[](const int i)           // lvalue version 
  { assert(i >= 0 && i < 3); return d[i]; }
  virtual void set(const double a, const double b, const double c);
  virtual void set(const double v[3]);
  double norm(void) const;                   // compute L_2 norm
  double normL1(void) const;	             // compute L_1 norm 
  double normLInf(void) const;	             // compute L_\infty norm
  void normalize(void);                      // destructively normalize self
  double dot(const Vector3&) const;          // dot product
  Vector3 cross(const Vector3&) const;       // cross product
  Vector3& operator+=(const Vector3&);
  Vector3& operator-=(const Vector3&);
  Vector3& operator*=(const double c);
  Vector3& operator/=(const double c);
  bool operator==(const Vector3& rhs) const;
  bool operator!=(const Vector3& rhs) const;
  bool operator<(const Vector3& rhs) const;   // checks if v1[i] < v2[i] for all i, true if all are true
  bool operator>(const Vector3& rhs) const;   // checks if v1[i] > v2[i] for all i, true if all are true
  const Vector3 mul_skew(const Vector3&) const;    // transpose(v)*skew(v2)
};


inline const Vector3 operator+(const Vector3& lhs, const Vector3& rhs);
inline const Vector3 operator-(const Vector3& lhs, const Vector3& rhs);
inline const Vector3 operator*(const double, const Vector3&);	// scalar multiply
inline const Vector3 operator*(const Vector3&, const double);
inline const Vector3 operator/(const Vector3& v, const double c);    // vect/scalar

double anglebetween(const Vector3&, const Vector3&);    // radians
double signedanglebetween(const Vector3 a, const Vector3 b, const Vector3 c); 
                                             // from a to b, wrt c in radians 


class Vector4 {
  friend ostream& operator<<(ostream& os, const Vector4& v);
  friend istream& operator>>(istream& is, Vector4& v); 
  double d[4];
 public:
  Vector4() {}
  Vector4(const double a, const double b, const double c, const double e);
  virtual const double& operator[](const int i) const  // inspection version
  { assert(i >= 0 && i < 4); return d[i]; }
  virtual double& operator[](const int i)           // lvalue version 
  { assert(i >= 0 && i < 4); return d[i]; }
  virtual void set(const double a, const double b, const double c, const double e);
  virtual void set(const double v[4]);
  double norm(void) const;                   // compute L_2 norm
  double normL1(void) const;	             // compute L_1 norm 
  double normLInf(void) const;	             // compute L_\infty norm
  void normalize(void);                      // destructively normalize self
  double dot(const Vector4&) const;          // dot product
  Vector4 cross(const Vector4&) const;       // cross product
  Vector4& operator+=(const Vector4&);
  Vector4& operator-=(const Vector4&);
  Vector4& operator*=(const double c);
  Vector4& operator/=(const double c);
  bool operator==(const Vector4& rhs) const;
  bool operator!=(const Vector4& rhs) const;
  bool operator<(const Vector4& rhs) const;   // checks if v1[i] < v2[i] for all i, true if all are true
  bool operator>(const Vector4& rhs) const;   // checks if v1[i] > v2[i] for all i, true if all are true
  const Vector4 mul_skew(const Vector4&) const;    // transpose(v)*skew(v2)
};


class Matrix3x3 {
  friend ostream& operator<<(ostream&, const Matrix3x3&);
  friend istream& operator>>(istream&, Matrix3x3&);  // untested
  double d[3][3];
 public:
  Matrix3x3() {};
  Matrix3x3(const double, const double, const double,
	    const double, const double, const double,
	    const double, const double, const double);
  void set(const double, const double, const double,
	   const double, const double, const double,
	   const double, const double, const double);
  void set2Identity(void);
  void set2Zeros(void);
  void setColV(const Vector3& c1, const Vector3& c2, const Vector3& c3);
  void setRowV(const Vector3& r1, const Vector3& r2, const Vector3& r3);
  void swapRows(const int r1, const int r2);
  void swapCols(const int c1, const int c2);
  const double* operator[](const int c) const;  // inspection value
  double* operator[](const int c);      // lvalue version 
  const Vector3 col(const int c) const;         // return column c as a vector
  const Vector3 row(const int r) const;         // return row r as a vector
  double trace(void) const;              // return the trace of the matrix
  double det(void) const;                // return the determinant of self
  const Matrix3x3 operator+(const Matrix3x3&) const;
  const Matrix3x3 operator-(const Matrix3x3&) const;
  const Matrix3x3 operator*(const Matrix3x3&) const;
  const Matrix3x3 mul_skew(const Vector3& v) const; // M*skew(v)
  const Matrix3x3 neg_mul_skew(const Vector3& v) const; // -M*skew(v)
  const Matrix3x3 transpose(void) const;       // return the transpose of self
  const Vector3 mul(const Vector3& v) const;   // M*v
  const Vector3 transpose_mul(const Vector3& v) const;    // transpose(M)*v
  const Matrix3x3 transpose_mul(const Matrix3x3&) const;  // transpose(M)*A
  void dumpto3x3(double (*dest)[3][3]) const;
  double normFrobenius(const int skipSqrt = 0) const;
  double normOff(const int skipSqrt = 0) const;
  unsigned int symEigen(Matrix3x3& eigVec, Vector3& eigVal, const double tol) const;
  void makeQuat(double& q0, double& q1, double& q2, double& q3) const;
  void makeQuat(double q[4]) const;
  void setFromQuat(const double q0, const double q1, const double q2,
		   const double q3);
  void setFromQuat(const double q[4]);
  // Matrix3x3& operator=(const QuatR4& q);
  bool operator==(const Matrix3x3& rhs) const;
  bool operator!=(const Matrix3x3& rhs) const;
};

inline const Matrix3x3 operator*(const Matrix3x3& M, const double c);
int Ax_b(Matrix3x3 M, Vector3& x, Vector3 b);
inline const Vector3 operator*(const Vector3&, const Matrix3x3&); // transpose(v)*M
inline const Vector3 operator*(const Matrix3x3&, const Vector3&); // M*v
inline void vEmXvPv(Vector3& v, const Matrix3x3& A, const Vector3& x, const
		    Vector3& b); // v = A * x + b

//
// Matrix4x4
//
class Matrix4x4 {
  friend ostream& operator<<(ostream&, const Matrix4x4&);
  friend istream& operator>>(istream&, Matrix4x4&);  // untested
  double d[4][4];
public:
  Matrix4x4() {};
  Matrix4x4(const double, const double, const double, const double,
	    const double, const double, const double, const double,
	    const double, const double, const double, const double,
	    const double, const double, const double, const double);
  void set(const double, const double, const double, const double,
	   const double, const double, const double, const double,
	   const double, const double, const double, const double,
	   const double, const double, const double, const double);
  void set2Identity(void);
  void set2Zeros(void);
  void setColV(const Vector4& c1, const Vector4& c2, const Vector4& c3, const Vector4& c4);
  void setRowV(const Vector4& r1, const Vector4& r2, const Vector4& r3, const Vector4& r4);
  void swapRows(const int r1, const int r2);
  void swapCols(const int c1, const int c2);
  const double* operator[](const int c) const;  // inspection value
  double* operator[](const int c);      // lvalue version 
  const Vector4 col(const int c) const;         // return column c as a vector
  const Vector4 row(const int r) const;         // return row r as a vector
  // double trace(void) const;              // return the trace of the matrix
  // double det(void) const;                // return the determinant of self
  const Matrix4x4 operator+(const Matrix4x4&) const;
  const Matrix4x4 operator-(const Matrix4x4&) const;
  const Matrix4x4 operator*(const Matrix4x4&) const;
  // const Matrix4x4 mul_skew(const Vector4& v) const; // M*skew(v)
  // const Matrix4x4 neg_mul_skew(const Vector4& v) const; // -M*skew(v)
  const Matrix4x4 transpose(void) const;       // return the transpose of self
  // const Vector4 transpose_mul(const Vector4& v) const;    // transpose(M)*v
  // const Matrix4x4 transpose_mul(const Matrix4x4&) const;  // transpose(M)*A
  void dumpto4x4(double (*dest)[4][4]) const;
  double normFrobenius(const int skipSqrt = 0) const;
  // double normOff(const int skipSqrt = 0) const;
  // unsigned int symEigen(Matrix4x4& eigVec, Vector4& eigVal, const double tol) const;
  // void makeQuat(double& q0, double& q1, double& q2, double& q4) const;
  // void makeQuat(double q[4]) const;
  // void setFromQuat(const double q0, const double q1, const double q2, const double q4);
  // void setFromQuat(const double q[4]);
  // Matrix4x4& operator=(const QuatR4& q);
  bool operator==(const Matrix4x4& rhs) const;
  bool operator!=(const Matrix4x4& rhs) const;

  void setTranslationV( const Vector3& t );  // when M4x4 is used as homogenous transform matrix
  void setRotationM( const Matrix3x3& m );  // when M4x4 is used as homogenous transform matrix

  void setRotationAboutX( const double angle_in_radians );  // when M4x4 is used as homogenous transform matrix
  void setRotationAboutY( const double angle_in_radians );  // when M4x4 is used as homogenous transform matrix
  void setRotationAboutZ( const double angle_in_radians );  // when M4x4 is used as homogenous transform matrix
};

inline const Matrix4x4 operator*(const Matrix4x4& M, const double c);
// int Ax_b(Matrix4x4 M, Vector4& x, Vector4 b);
inline const Vector4 operator*(const Vector4&, const Matrix4x4&); // transpose(v)*M
inline const Vector4 operator*(const Matrix4x4&, const Vector4&); // M*v
inline void vEmXvPv(Vector4& v, const Matrix4x4& A, const Vector4& x, const
		    Vector4& b); // v = A * x + b


// Implementation - Vector2

inline Vector2::Vector2(const double a, const double b)
{ d[0] = a; d[1] = b; }

inline bool Vector2::operator==(const Vector2& rhs) const
{ 
  return (d[0] == rhs.d[0] && d[1] == rhs.d[1]);
}

inline bool Vector2::operator!=(const Vector2& rhs) const
{
  return !(operator==(rhs));
}

inline bool Vector2::operator<(const Vector2& rhs) const
{
  return (d[0] < rhs.d[0] && d[1] < rhs.d[1]);
}

inline bool Vector2::operator>(const Vector2& rhs) const
{
  return (d[0] > rhs.d[0] && d[1] > rhs.d[1]);
}

inline void 
Vector2::set(const double a, const double b)
{
  d[0]=a;
  d[1]=b;
}

inline void 
Vector2::set(const double v[2])
{
  d[0]=v[0];
  d[1]=v[1];
}

inline double 
Vector2::norm(void) const
{
  return sqrt(d[0]*d[0] + d[1]*d[1]);
}

inline double
Vector2::normL1(void) const
{
  return myfabs(d[0]) + myfabs(d[1]);
}

inline double
Vector2::normLInf(void) const
{
  double vmax = myfabs(d[0]);
  register double tmp = myfabs(d[1]);
  if(tmp > vmax) vmax = tmp;
  return vmax;
}

inline void 
Vector2::normalize(void)
{
  const double n = norm();
  if(n > 0.0)
    {
      d[0] /= n;
      d[1] /= n;
    }
}

inline double 
Vector2::dot(const Vector2& v) const	// dot product
{
  return double(d[0]*v[0] + d[1]*v[1]);
}

inline double
Vector2::cross(const Vector2& v) const        // cross product
{
  return double(d[0] * v[1] - d[1] * v[0]);
}

inline Vector2&
Vector2::operator+=(const Vector2& v)
{
  d[0] += v[0];
  d[1] += v[1];
  return *this;
}

inline const Vector2
operator+(const Vector2& lhs, const Vector2& rhs)
{
  return Vector2(lhs) += rhs;
}

inline Vector2&
Vector2::operator-=(const Vector2& v)
{
  d[0] -= v[0];
  d[1] -= v[1];
  return *this;
}

inline const Vector2
operator-(const Vector2& lhs, const Vector2& rhs)
{
  return Vector2(lhs) -= rhs;
}


inline Vector3::Vector3(const double a, const double b, const double c)
{ d[0] = a; d[1] = b; d[2] = c; }

inline Vector4::Vector4(const double a, const double b, const double c, const double e)
{ d[0] = a; d[1] = b; d[2] = c; d[3] = e; }



// matrix 3x3
inline const double* Matrix3x3::operator[](const int c) const  
{ assert(c >= 0 && c < 3); return d[c]; }
inline double* Matrix3x3::operator[](const int c) 
{ assert(c >= 0 && c < 3); return d[c]; }

inline const Vector3 Matrix3x3::col(const int c) const
{ return Vector3(d[0][c], d[1][c], d[2][c]); }

inline const Vector3 Matrix3x3::row(const int r) const
{ return Vector3(d[r][0], d[r][1], d[r][2]); }


// matrix 4x4
inline const double* Matrix4x4::operator[](const int c) const  
{ assert(c >= 0 && c < 4); return d[c]; }
inline double* Matrix4x4::operator[](const int c) 
{ assert(c >= 0 && c < 4); return d[c]; }

inline const Vector4 Matrix4x4::col(const int c) const
{ return Vector4(d[0][c], d[1][c], d[2][c], d[3][c]); }

inline const Vector4 Matrix4x4::row(const int r) const
{ return Vector4(d[r][0], d[r][1], d[r][2], d[r][3]); }


// Implementation - Vector3

inline bool Vector3::operator==(const Vector3& rhs) const
{ 
  return (d[0] == rhs.d[0] && d[1] == rhs.d[1] && d[2] == rhs.d[2]);
}

inline bool Vector3::operator!=(const Vector3& rhs) const
{
  return !(operator==(rhs));
}

inline bool Vector3::operator<(const Vector3& rhs) const
{
  return (d[0] < rhs.d[0] && d[1] < rhs.d[1] && d[2] < rhs.d[2]);
}

inline bool Vector3::operator>(const Vector3& rhs) const
{
  return (d[0] > rhs.d[0] && d[1] > rhs.d[1] && d[2] > rhs.d[2]);
}

inline void 
Vector3::set(const double a, const double b, const double c)
{
  d[0]=a;
  d[1]=b;
  d[2]=c;
}

inline void 
Vector3::set(const double v[3])
{
  d[0]=v[0];
  d[1]=v[1];
  d[2]=v[2];
}

inline double 
Vector3::norm(void) const
{
  return sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
}

inline double
Vector3::normL1(void) const
{
  return myfabs(d[0]) + myfabs(d[1]) + myfabs(d[2]);
}

inline double
Vector3::normLInf(void) const
{
  double vmax = myfabs(d[0]);
  register double tmp = myfabs(d[1]);
  if(tmp > vmax) vmax = tmp;
  tmp = myfabs(d[2]);
  if(tmp > vmax) vmax = tmp;
  return vmax;
}

inline void 
Vector3::normalize(void)
{
  const double n = norm();
  if(n > 0.0)
    {
      d[0] /= n;
      d[1] /= n;
      d[2] /= n;
    }
}


inline double 
Vector3::dot(const Vector3& v) const	// dot product
{
  return double(d[0]*v[0] + d[1]*v[1] + d[2]*v[2]);
}


inline Vector3 
Vector3::cross(const Vector3& v) const        // cross product
{
  return Vector3(d[1] * v[2] - d[2] * v[1],
		 d[2] * v[0] - d[0] * v[2],
		 d[0] * v[1] - d[1] * v[0]);
}


inline const Vector3 
Vector3::mul_skew(const Vector3& v) const        // transpose(d)*skew(v)
{
  return Vector3(d[1] * v[2] - d[2] * v[1],
		 d[2] * v[0] - d[0] * v[2],
		 d[0] * v[1] - d[1] * v[0]);
}


inline Vector3&
Vector3::operator+=(const Vector3& v)
{
  d[0] += v[0];
  d[1] += v[1];
  d[2] += v[2];
  return *this;
}



inline const Vector3 
operator+(const Vector3& lhs, const Vector3& rhs)
{
  return Vector3(lhs) += rhs;
}



inline Vector3&
Vector3::operator-=(const Vector3& v)
{
  d[0] -= v[0];
  d[1] -= v[1];
  d[2] -= v[2];
  return *this;
}


inline const Vector3
operator-(const Vector3& lhs, const Vector3& rhs)
{
  return Vector3(lhs) -= rhs;
}


// Implementation - Vector4

inline bool Vector4::operator==(const Vector4& rhs) const
{ 
  return (d[0] == rhs.d[0] && d[1] == rhs.d[1] && d[2] == rhs.d[2] && d[3] == rhs.d[3]);
}

inline bool Vector4::operator!=(const Vector4& rhs) const
{
  return !(operator==(rhs));
}

inline bool Vector4::operator<(const Vector4& rhs) const
{
  return (d[0] < rhs.d[0] && d[1] < rhs.d[1] && d[2] < rhs.d[2] && d[3] < rhs.d[3]);
}

inline bool Vector4::operator>(const Vector4& rhs) const
{
  return (d[0] > rhs.d[0] && d[1] > rhs.d[1] && d[2] > rhs.d[2] && d[3] > rhs.d[3]);
}

inline void 
Vector4::set(const double a, const double b, const double c, const double e)
{
  d[0]=a;
  d[1]=b;
  d[2]=c;
  d[3]=e;
}

inline void 
Vector4::set(const double v[4])
{
  d[0]=v[0];
  d[1]=v[1];
  d[2]=v[2];
  d[3]=v[3];
}

inline double 
Vector4::norm(void) const
{
  return sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2] + d[3]*d[3]);
}

inline double
Vector4::normL1(void) const
{
  return myfabs(d[0]) + myfabs(d[1]) + myfabs(d[2]) + myfabs(d[3]);
}

inline double
Vector4::normLInf(void) const
{
  double vmax = myfabs(d[0]);
  register double tmp = myfabs(d[1]);
  if(tmp > vmax) vmax = tmp;
  tmp = myfabs(d[2]);
  if(tmp > vmax) vmax = tmp;
  tmp = myfabs(d[3]);
  if(tmp > vmax) vmax = tmp;
  return vmax;
}

inline void 
Vector4::normalize(void)
{
  const double n = norm();
  if(n > 0.0)
    {
      d[0] /= n;
      d[1] /= n;
      d[2] /= n;
      d[3] /= n;
    }
}

inline double 
Vector4::dot(const Vector4& v) const	// dot product
{
  return double(d[0]*v[0] + d[1]*v[1] + d[2]*v[2] + d[3]*v[3]);
}

inline Vector4
Vector4::cross(const Vector4& v) const        // cross product
{
  std::cout << "Vector4::cross(Vector4) NOT IMPLEMENTED!" << std::endl;
  assert( true == false );
  return Vector4(0,0,0,0);
}

inline Vector4&
Vector4::operator+=(const Vector4& v)
{
  d[0] += v[0];
  d[1] += v[1];
  d[2] += v[2];
  d[3] += v[3];
  return *this;
}

inline const Vector4
operator+(const Vector4& lhs, const Vector4& rhs)
{
  return Vector4(lhs) += rhs;
}

inline Vector4&
Vector4::operator-=(const Vector4& v)
{
  d[0] -= v[0];
  d[1] -= v[1];
  d[2] -= v[2];
  d[3] -= v[3];
  return *this;
}

inline const Vector4
operator-(const Vector4& lhs, const Vector4& rhs)
{
  return Vector4(lhs) -= rhs;
}


// Implementation - Matrix3x3

inline const Vector3
operator*(const Vector3& v, const Matrix3x3& m) 
{
  return Vector3(v[0]*m[0][0] + v[1]*m[1][0] + v[2]*m[2][0],
                 v[0]*m[0][1] + v[1]*m[1][1] + v[2]*m[2][1],
                 v[0]*m[0][2] + v[1]*m[1][2] + v[2]*m[2][2]);
}

inline const Vector3 
operator*(const Matrix3x3& m, const Vector3& v) 
{
  return Vector3(m[0][0]*v[0]+m[0][1]*v[1]+m[0][2]*v[2],
		 m[1][0]*v[0]+m[1][1]*v[1]+m[1][2]*v[2],
		 m[2][0]*v[0]+m[2][1]*v[1]+m[2][2]*v[2]);
}


inline const Vector3 
Matrix3x3::transpose_mul(const Vector3& v) const
{
  return Vector3(d[0][0]*v[0] + d[1][0]*v[1] + d[2][0]*v[2], 
		 d[0][1]*v[0] + d[1][1]*v[1] + d[2][1]*v[2], 
		 d[0][2]*v[0] + d[1][2]*v[1] + d[2][2]*v[2]);
}

inline const Vector3
Matrix3x3::mul(const Vector3& v) const
{
  return Vector3(d[0][0]*v[0] + d[0][1]*v[1] + d[0][2]*v[2],
                 d[1][0]*v[0] + d[1][1]*v[1] + d[1][2]*v[2],
                 d[2][0]*v[0] + d[2][1]*v[1] + d[2][2]*v[2]);
}


// Implementation - Matrix4x4

inline const Vector4
operator*(const Vector4& v, const Matrix4x4& m) 
{
  return Vector4(v[0]*m[0][0] + v[1]*m[1][0] + v[2]*m[2][0] + v[3]*m[2][0],
                 v[0]*m[0][1] + v[1]*m[1][1] + v[2]*m[2][1] + v[3]*m[2][1],
                 v[0]*m[0][2] + v[1]*m[1][2] + v[2]*m[2][2] + v[3]*m[2][2],
		 v[0]*m[0][3] + v[1]*m[1][3] + v[2]*m[2][3] + v[3]*m[2][3]);
}

inline const Vector4
operator*(const Matrix4x4& m, const Vector4& v) 
{
  return Vector4(m[0][0]*v[0]+m[0][1]*v[1]+m[0][2]*v[2]+m[0][3]*v[3],
		 m[1][0]*v[0]+m[1][1]*v[1]+m[1][2]*v[2]+m[1][3]*v[3],
		 m[2][0]*v[0]+m[2][1]*v[1]+m[2][2]*v[2]+m[2][3]*v[3],
		 m[3][0]*v[0]+m[3][1]*v[1]+m[3][2]*v[2]+m[3][3]*v[3]);
}



inline Vector3&
Vector3::operator*=(const double c)
{
  d[0] *= c; d[1] *= c; d[2] *= c;
  return *this;
}


inline Vector3&
Vector3::operator/=(const double c)
{
  d[0] /= c; d[1] /= c; d[2] /= c;
  return *this;
}


inline const Vector3 
operator*(const double c, const Vector3& v)
{
  return Vector3(v[0]*c, v[1]*c, v[2]*c);
}


inline const Vector3 
operator*(const Vector3& v, const double c)
{
  return Vector3(v[0]*c, v[1]*c, v[2]*c);
}


inline const Vector3 
operator/(const Vector3& v, const double c)
{
  return Vector3(v[0]/c, v[1]/c, v[2]/c);
}




inline
Matrix3x3::Matrix3x3(const double a11, const double a12, const double a13,
		     const double a21, const double a22, const double a23,
		     const double a31, const double a32, const double a33)
{
  d[0][0]=a11;  d[0][1]=a12;  d[0][2]=a13;
  d[1][0]=a21;  d[1][1]=a22;  d[1][2]=a23;
  d[2][0]=a31;  d[2][1]=a32;  d[2][2]=a33;
}


inline void
Matrix3x3::set(const double a11, const double a12, const double a13,
	       const double a21, const double a22, const double a23,
	       const double a31, const double a32, const double a33)
{
  d[0][0]=a11;  d[0][1]=a12;  d[0][2]=a13;
  d[1][0]=a21;  d[1][1]=a22;  d[1][2]=a23;
  d[2][0]=a31;  d[2][1]=a32;  d[2][2]=a33;
}




inline void Matrix3x3::set2Identity(void)
{
  d[0][0] = ONE;  d[0][1] = ZERO; d[0][2] = ZERO;
  d[1][0] = ZERO; d[1][1] = ONE;  d[1][2] = ZERO;
  d[2][0] = ZERO; d[2][1] = ZERO; d[2][2] = ONE;
}


inline const Matrix3x3 Matrix3x3::transpose(void) const
{
  return Matrix3x3(
	   d[0][0], d[1][0], d[2][0],
	   d[0][1], d[1][1], d[2][1],
	   d[0][2], d[1][2], d[2][2]);
}

inline double Matrix3x3::det(void) const
{
// Cofactor expansion along the first column.
  return d[0][0]*(d[1][1]*d[2][2]-d[1][2]*d[2][1]) +
         d[1][0]*(d[0][2]*d[2][1]-d[0][1]*d[2][2]) +
	 d[2][0]*(d[0][1]*d[1][2]-d[0][2]*d[1][1]);
}


inline double 
Matrix3x3::trace(void) const
{
  return d[0][0]+d[1][1]+d[2][2];
}


inline const Matrix3x3 
Matrix3x3::operator+(const Matrix3x3& m) const
{
  return Matrix3x3(
	   d[0][0]+m[0][0],  d[0][1]+m[0][1],   d[0][2]+m[0][2],
	   d[1][0]+m[1][0],  d[1][1]+m[1][1],   d[1][2]+m[1][2],
	   d[2][0]+m[2][0],  d[2][1]+m[2][1],   d[2][2]+m[2][2]);
}


inline const Matrix3x3 
Matrix3x3::operator-(const Matrix3x3& m) const
{
  return Matrix3x3(
	   d[0][0]-m[0][0],  d[0][1]-m[0][1],   d[0][2]-m[0][2],
	   d[1][0]-m[1][0],  d[1][1]-m[1][1],   d[1][2]-m[1][2],
	   d[2][0]-m[2][0],  d[2][1]-m[2][1],   d[2][2]-m[2][2]);
}

inline const Matrix3x3
operator*(const Matrix3x3& M, const double c)
{
  return Matrix3x3(M[0][0]*c, M[0][1]*c, M[0][2]*c,
		   M[1][0]*c, M[1][1]*c, M[1][2]*c,
		   M[2][0]*c, M[2][1]*c, M[2][2]*c);
}


inline void 
Matrix3x3::setColV(const Vector3& c1, const Vector3& c2, const Vector3& c3)
{
  d[0][0] = c1[0];  d[0][1] = c2[0];  d[0][2] = c3[0];
  d[1][0] = c1[1];  d[1][1] = c2[1];  d[1][2] = c3[1];
  d[2][0] = c1[2];  d[2][1] = c2[2];  d[2][2] = c3[2];
}


inline void 
Matrix3x3::setRowV(const Vector3& r1, const Vector3& r2, const Vector3& r3)
{
  d[0][0] = r1[0];  d[0][1] = r1[1];  d[0][2] = r1[2];
  d[1][0] = r2[0];  d[1][1] = r2[1];  d[1][2] = r2[2];
  d[2][0] = r3[0];  d[2][1] = r3[1];  d[2][2] = r3[2];
}

inline void
Matrix3x3::makeQuat(double q[4]) const
{
  makeQuat(q[0], q[1], q[2], q[3]);
}

inline void
Matrix3x3::setFromQuat(const double q[4])
{
  setFromQuat(q[0], q[1], q[2], q[3]);
}


// This is ugly, but it is the fastest, most compiler friendly way I 
// know of doing this.  No temporaries
inline void vEmXvPv(Vector3& v, const Matrix3x3& A, const Vector3& x, 
		    const Vector3& b)
{
  v[0] = A[0][0] * x[0] + A[0][1] * x[1] + A[0][2] * x[2] + b[0];
  v[1] = A[1][0] * x[0] + A[1][1] * x[1] + A[1][2] * x[2] + b[1];
  v[2] = A[2][0] * x[0] + A[2][1] * x[1] + A[2][2] * x[2] + b[2];
}


inline bool Matrix3x3::operator==(const Matrix3x3& rhs) const
{

  return (d[0][0] == rhs.d[0][0] && d[0][1] == rhs.d[0][1] && 
	  d[0][2] == rhs.d[0][2] && d[1][0] == rhs.d[1][0] && 
	  d[1][1] == rhs.d[1][1] && d[1][2] == rhs.d[1][2] && 
	  d[2][0] == rhs.d[2][0] && d[2][1] == rhs.d[2][1] && 
	  d[2][2] == rhs.d[2][2]); 
}

inline bool Matrix3x3::operator!=(const Matrix3x3& rhs) const
{
  return !(operator==(rhs));
}

//
// Matrix4x4
//

inline
Matrix4x4::Matrix4x4(const double a11, const double a12, const double a13, const double a14,
		     const double a21, const double a22, const double a23, const double a24,
		     const double a31, const double a32, const double a33, const double a34,
		     const double a41, const double a42, const double a43, const double a44)
{
  d[0][0]=a11;  d[0][1]=a12;  d[0][2]=a13; d[0][3]=a14;
  d[1][0]=a21;  d[1][1]=a22;  d[1][2]=a23; d[1][3]=a24;
  d[2][0]=a31;  d[2][1]=a32;  d[2][2]=a33; d[2][3]=a34;
  d[3][0]=a41;  d[3][1]=a42;  d[3][2]=a43; d[3][3]=a44;
}


inline void
Matrix4x4::set(const double a11, const double a12, const double a13, const double a14,
	       const double a21, const double a22, const double a23, const double a24,
	       const double a31, const double a32, const double a33, const double a34,
	       const double a41, const double a42, const double a43, const double a44)
{
  d[0][0]=a11;  d[0][1]=a12;  d[0][2]=a13;  d[0][3]=a14;
  d[1][0]=a21;  d[1][1]=a22;  d[1][2]=a23;  d[1][3]=a24;
  d[2][0]=a31;  d[2][1]=a32;  d[2][2]=a33;  d[2][3]=a34;
  d[3][0]=a41;  d[3][1]=a42;  d[3][2]=a43;  d[3][3]=a44;
}

inline void Matrix4x4::set2Identity(void)
{
  d[0][0] = ONE;  d[0][1] = ZERO; d[0][2] = ZERO; d[0][3] = ZERO;
  d[1][0] = ZERO; d[1][1] = ONE;  d[1][2] = ZERO; d[1][3] = ZERO;
  d[2][0] = ZERO; d[2][1] = ZERO; d[2][2] = ONE;  d[2][3] = ZERO;
  d[3][0] = ZERO; d[3][1] = ZERO; d[3][2] = ZERO; d[3][3] = ONE;
}

inline const Matrix4x4 Matrix4x4::transpose(void) const
{
  return Matrix4x4(d[0][0], d[1][0], d[2][0], d[3][0],
		   d[0][1], d[1][1], d[2][1], d[3][1],
		   d[0][2], d[1][2], d[2][2], d[3][2],
		   d[0][3], d[1][3], d[2][3], d[3][3]);
}

inline const Matrix4x4
Matrix4x4::operator+(const Matrix4x4& m) const
{
  return Matrix4x4( d[0][0]+m[0][0],  d[0][1]+m[0][1],   d[0][2]+m[0][2],  d[0][3]+m[0][3],
		    d[1][0]+m[1][0],  d[1][1]+m[1][1],   d[1][2]+m[1][2],  d[1][3]+m[1][3],
		    d[2][0]+m[2][0],  d[2][1]+m[2][1],   d[2][2]+m[2][2],  d[2][3]+m[2][3],
		    d[3][0]+m[3][0],  d[3][1]+m[3][1],   d[3][2]+m[3][2],  d[3][3]+m[3][3] );
}

inline const Matrix4x4
Matrix4x4::operator-(const Matrix4x4& m) const
{
  return Matrix4x4( d[0][0]-m[0][0],  d[0][1]-m[0][1],   d[0][2]-m[0][2],   d[0][3]-m[0][3],
		    d[1][0]-m[1][0],  d[1][1]-m[1][1],   d[1][2]-m[1][2],   d[1][3]-m[1][3],
		    d[2][0]-m[2][0],  d[2][1]-m[2][1],   d[2][2]-m[2][2],   d[2][3]-m[2][3],
		    d[3][0]-m[3][0],  d[3][1]-m[3][1],   d[3][2]-m[3][2],   d[3][3]-m[3][3] );
}

inline const Matrix4x4
operator*(const Matrix4x4& M, const double c)
{
  return Matrix4x4(M[0][0]*c, M[0][1]*c, M[0][2]*c, M[0][3]*c,
		   M[1][0]*c, M[1][1]*c, M[1][2]*c, M[1][3]*c,
		   M[2][0]*c, M[2][1]*c, M[2][2]*c, M[2][3]*c,
		   M[3][0]*c, M[3][1]*c, M[3][2]*c, M[3][3]*c);
}


inline void 
Matrix4x4::setColV(const Vector4& c1, const Vector4& c2, const Vector4& c3, const Vector4& c4)
{
  d[0][0] = c1[0];  d[0][1] = c2[0];  d[0][2] = c3[0];  d[0][3] = c3[0];
  d[1][0] = c1[1];  d[1][1] = c2[1];  d[1][2] = c3[1];  d[1][3] = c3[1];
  d[2][0] = c1[2];  d[2][1] = c2[2];  d[2][2] = c3[2];  d[2][3] = c3[2];
  d[3][0] = c1[3];  d[3][1] = c2[3];  d[3][2] = c3[3];  d[3][3] = c3[3];
}

inline void 
Matrix4x4::setRowV(const Vector4& r1, const Vector4& r2, const Vector4& r3, const Vector4& r4)
{
  d[0][0] = r1[0];  d[0][1] = r1[1];  d[0][2] = r1[2];  d[0][3] = r1[3];
  d[1][0] = r2[0];  d[1][1] = r2[1];  d[1][2] = r2[2];  d[1][3] = r2[3];
  d[2][0] = r3[0];  d[2][1] = r3[1];  d[2][2] = r3[2];  d[2][3] = r3[3];
  d[3][0] = r4[0];  d[3][1] = r4[1];  d[3][2] = r4[2];  d[3][3] = r4[3];
}

// This is ugly, but it is the fastest, most compiler friendly way I 
// know of doing this.  No temporaries
inline void vEmXvPv(Vector4& v, const Matrix4x4& A, const Vector4& x, 
		    const Vector4& b)
{
  v[0] = A[0][0] * x[0] + A[0][1] * x[1] + A[0][2] * x[2] + A[0][2] * x[3] + b[0];
  v[1] = A[1][0] * x[0] + A[1][1] * x[1] + A[1][2] * x[2] + A[1][2] * x[3] + b[1];
  v[2] = A[2][0] * x[0] + A[2][1] * x[1] + A[2][2] * x[2] + A[2][2] * x[3] + b[2];
  v[3] = A[3][0] * x[0] + A[3][1] * x[1] + A[3][2] * x[2] + A[3][2] * x[3] + b[3];
}

inline bool Matrix4x4::operator==(const Matrix4x4& rhs) const
{
  return (d[0][0] == rhs.d[0][0] && d[0][1] == rhs.d[0][1] && d[0][2] == rhs.d[0][2] && d[0][3] == rhs.d[0][3] && 
	  d[1][0] == rhs.d[1][0] && d[1][1] == rhs.d[1][1] && d[1][2] == rhs.d[1][2] && d[1][3] == rhs.d[1][3] && 
	  d[2][0] == rhs.d[2][0] && d[2][1] == rhs.d[2][1] && d[2][2] == rhs.d[2][2] && d[2][3] == rhs.d[2][3] && 
	  d[3][0] == rhs.d[3][0] && d[3][1] == rhs.d[3][1] && d[3][2] == rhs.d[3][2] && d[3][3] == rhs.d[3][3]);
}

inline bool Matrix4x4::operator!=(const Matrix4x4& rhs) const
{
  return !(operator==(rhs));
}

// when M4x4 is used as homogenous transform matrix
inline void Matrix4x4::setTranslationV( const Vector3& t ) 
{
  d[0][3] = t[0];
  d[1][3] = t[1];
  d[2][3] = t[2];
}

// when M4x4 is used as homogenous transform matrix
inline void Matrix4x4::setRotationM( const Matrix3x3& m )
{
  d[0][0] = m[0][0];  d[0][1] = m[0][1];  d[0][2] = m[0][2];
  d[1][0] = m[1][0];  d[1][1] = m[1][1];  d[1][2] = m[1][2];
  d[2][0] = m[2][0];  d[2][1] = m[2][1];  d[2][2] = m[2][2];
}

// when M4x4 is used as homogenous transform matrix
inline void Matrix4x4::setRotationAboutX( const double rad )
{
  d[0][0] = 1;  d[0][1] = 0;          d[0][2] = 0;          d[0][3] = 0;
  d[1][0] = 0;  d[1][1] = cos(rad);   d[1][2] = -sin(rad);  d[1][3] = 0;
  d[2][0] = 0;  d[2][1] = sin(rad);   d[2][2] = cos(rad);   d[2][3] = 0;
  d[3][0] = 0;  d[3][1] = 0;          d[3][2] = 0;          d[3][3] = 1;
}

// when M4x4 is used as homogenous transform matrix
inline void Matrix4x4::setRotationAboutY( const double rad )
{
  d[0][0] = cos(rad);   d[0][1] = 0;  d[0][2] = sin(rad);  d[0][3] = 0;
  d[1][0] = 0;          d[1][1] = 1;  d[1][2] = 0;         d[1][3] = 0;
  d[2][0] = -sin(rad);  d[2][1] = 0;  d[2][2] = cos(rad);  d[2][3] = 0;
  d[3][0] = 0;          d[3][1] = 0;  d[3][2] = 0;         d[3][3] = 1;
}

// when M4x4 is used as homogenous transform matrix
inline void Matrix4x4::setRotationAboutZ( const double rad )
{
  d[0][0] = cos(rad);  d[0][1] = -sin(rad);  d[0][2] = 0;  d[0][3] = 0;
  d[1][0] = sin(rad);  d[1][1] = cos(rad);   d[1][2] = 0;  d[1][3] = 0;
  d[2][0] = 0;         d[2][1] = 0;          d[2][2] = 1;  d[2][3] = 0;
  d[3][0] = 0;         d[3][1] = 0;          d[3][2] = 0;  d[3][3] = 1;
}


void glErrorCheck( const char* msg );

inline double normalize( double *v )
{
  double mag = sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
  v[0] /= mag;
  v[1] /= mag;
  v[2] /= mag;
  return mag;
}

#endif // GRAPHICS_UTIL_H
