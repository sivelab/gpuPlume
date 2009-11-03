// #include <iostream>

#include "graphicsUtil.h"
#include "rgb.h"


void glErrorCheck( const char* msg )
{
#ifndef NDEBUG
  GLenum err_code;
  err_code = glGetError();
  while (err_code != GL_NO_ERROR) {
    // std::cerr << "OpenGL Error: " << gluErrorString(err_code) << ", Context[" << msg << "]" << std::std::endl;
    err_code = glGetError();
  }
#endif
}

Texture::Texture()
{
}

Texture::Texture( std::string texture_filename )
  : _tex_filename( texture_filename ), _tex_id(0)
{
}

void Texture::load()
{
  RGBImageRec *tex_data;
  tex_data = rgbImageLoad( _tex_filename.c_str() );
  if (tex_data) {
    // std::std::cout << "Loaded \"" << _tex_filename << "\"..." << std::std::endl;

    glGenTextures(1, &_tex_id);
    glErrorCheck("glGenTextures");

    glBindTexture(GL_TEXTURE_2D, _tex_id);   // 2d texture (x and y size)
    glErrorCheck("glBindTexture");

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, tex_data->sizeX, tex_data->sizeY,
		 0, GL_RGB, GL_UNSIGNED_BYTE, tex_data->data);
    glErrorCheck("glTexImage2D");

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glErrorCheck("glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT)");
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glErrorCheck("glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT)");

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glErrorCheck("glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)");

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glErrorCheck("glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)");

    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
    glErrorCheck("glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL)");
  }
}

void Texture::bind() const
{
  glBindTexture(GL_TEXTURE_2D, _tex_id);   // 2d texture (x and y size)
  glErrorCheck("glBindTexture");
}



// ------------------------------------------------------------------------
// Vector Stuff
// ------------------------------------------------------------------------

ostream &
operator<<(ostream& os, const Vector2& v)
{
  os << '[' << v[0] << ' ' << v[1] << ']';
  return os;
}

istream &
operator>>(istream& is, Vector2& v)
{
  is >> v[0] >> v[1];
  return is;
}


ostream &
operator<<(ostream& os, const Vector3& v)
{
  os << '[' << v[0] << ' ' << v[1] << ' ' << v[2] << ']';
  return os;
}

istream &
operator>>(istream& is, Vector3& v)
{
  is >> v[0] >> v[1] >> v[2];
  return is;
}

ostream &
operator<<(ostream& os, const Vector4& v)
{
  os << '[' << v[0] << ' ' << v[1] << ' ' << v[2] << ' ' << v[3] << ']';
  return os;
}

istream &
operator>>(istream& is, Vector4& v)
{
  is >> v[0] >> v[1] >> v[2] >> v[3];
  return is;
}

double
anglebetween(const Vector2& u, const Vector2& v)
{
  std::cout << "anglebetween(Vector2, Vector2) in linearR3.C - NOT COMPLETE" << std::endl;
  assert( true == false );
  return 0.0;
// Obvious function:
//    acos(u.dot(v)/(u.norm() * v.norm()))
// I think this may be faster though.  Traded a sqrt for a double and a mul.
}


double
signedanglebetween(const Vector2 a, const Vector2 b, const Vector3 c)
{
  std::cout << "signedanglebetween(Vector2, Vector2, Vector2) in linearR3.C - NOT COMPLETE" << std::endl;
  assert( true == false );
  return 0.0;
}


double
anglebetween(const Vector3& u, const Vector3& v)
{
// Obvious function:
//    acos(u.dot(v)/(u.norm() * v.norm()))
// I think this may be faster though.  Traded a sqrt for a double and a mul.

  double temp=(u[0]*u[0]+u[1]*u[1]+u[2]*u[2])*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  return acos(((u[0]*v[0]+u[1]*v[1]+u[2]*v[2])*sqrt(temp))/temp);
}


double
signedanglebetween(const Vector3 a, const Vector3 b, const Vector3 c)
{
  // Faden wrote this
  Vector3 axb= a.cross(b);        // |a x b| = sin(theta) * |a| * |b|
  double sinAlpha= (axb).norm();  // I'm assuming a, b are normalized.
  if (sinAlpha > ONE)
    sinAlpha = ONE;
  if (sinAlpha < -ONE)
    sinAlpha = -ONE;
  double Alpha= asin( sinAlpha );
  if (axb.dot(c) < 0)
    Alpha= -Alpha; // check if antiparallel
//  cerr << "signed angle between=" << Alpha*180/3.14 << std::endl;
  return Alpha;
}



// ------------------------------------------------------------------------
// Matrix3x3 Stuff
// ------------------------------------------------------------------------


ostream
&operator <<(ostream& os, const Matrix3x3& m)
{
  os << "[[" << m[0][0] << ' ' << m[0][1] << ' ' << m[0][2] << "];...\n"
     << " [" << m[1][0] << ' ' << m[1][1] << ' ' << m[1][2] << "];...\n"
     << " [" << m[2][0] << ' ' << m[2][1] << ' ' << m[2][2] << "]]";
  return os;
}


istream
&operator >>(istream& is, Matrix3x3& v)
{
  is >> v[0][0] >> v[0][1] >> v[0][2]
     >> v[1][0] >> v[1][1] >> v[1][2]
     >> v[2][0] >> v[2][1] >> v[2][2];
  return is;
}


const Matrix3x3
Matrix3x3::operator*(const Matrix3x3& m) const
{
  return Matrix3x3(
	   d[0][0]*m[0][0]+d[0][1]*m[1][0]+d[0][2]*m[2][0],
	   d[0][0]*m[0][1]+d[0][1]*m[1][1]+d[0][2]*m[2][1],
	   d[0][0]*m[0][2]+d[0][1]*m[1][2]+d[0][2]*m[2][2],
	   d[1][0]*m[0][0]+d[1][1]*m[1][0]+d[1][2]*m[2][0],
	   d[1][0]*m[0][1]+d[1][1]*m[1][1]+d[1][2]*m[2][1],
	   d[1][0]*m[0][2]+d[1][1]*m[1][2]+d[1][2]*m[2][2],
	   d[2][0]*m[0][0]+d[2][1]*m[1][0]+d[2][2]*m[2][0],
	   d[2][0]*m[0][1]+d[2][1]*m[1][1]+d[2][2]*m[2][1],
	   d[2][0]*m[0][2]+d[2][1]*m[1][2]+d[2][2]*m[2][2]);
}



const Matrix3x3
Matrix3x3::mul_skew(const Vector3& v) const
{
// A*skew(v)
  return Matrix3x3(
       d[0][1]*v[2] - d[0][2]*v[1],
      -d[0][0]*v[2] + d[0][2]*v[0],
       d[0][0]*v[1] - d[0][1]*v[0],
       d[1][1]*v[2] - d[1][2]*v[1],
      -d[1][0]*v[2] + d[1][2]*v[0],
       d[1][0]*v[1] - d[1][1]*v[0],
       d[2][1]*v[2] - d[2][2]*v[1],
      -d[2][0]*v[2] + d[2][2]*v[0],
       d[2][0]*v[1] - d[2][1]*v[0]);
}


const Matrix3x3
Matrix3x3::neg_mul_skew(const Vector3& v) const
{
// -A*skew(v)

  return Matrix3x3(
      -d[0][1]*v[2] + d[0][2]*v[1],
       d[0][0]*v[2] - d[0][2]*v[0],
      -d[0][0]*v[1] + d[0][1]*v[0],
      -d[1][1]*v[2] + d[1][2]*v[1],
       d[1][0]*v[2] - d[1][2]*v[0],
      -d[1][0]*v[1] + d[1][1]*v[0],
      -d[2][1]*v[2] + d[2][2]*v[1],
       d[2][0]*v[2] - d[2][2]*v[0],
      -d[2][0]*v[1] + d[2][1]*v[0]);
}



const Matrix3x3
Matrix3x3::transpose_mul(const Matrix3x3& m) const
{
  return Matrix3x3(
           d[0][0]*m[0][0] + d[1][0]*m[1][0] + d[2][0]*m[2][0],
           d[0][0]*m[0][1] + d[1][0]*m[1][1] + d[2][0]*m[2][1],
           d[0][0]*m[0][2] + d[1][0]*m[1][2] + d[2][0]*m[2][2],
           d[0][1]*m[0][0] + d[1][1]*m[1][0] + d[2][1]*m[2][0],
           d[0][1]*m[0][1] + d[1][1]*m[1][1] + d[2][1]*m[2][1],
           d[0][1]*m[0][2] + d[1][1]*m[1][2] + d[2][1]*m[2][2],
           d[0][2]*m[0][0] + d[1][2]*m[1][0] + d[2][2]*m[2][0],
           d[0][2]*m[0][1] + d[1][2]*m[1][1] + d[2][2]*m[2][1],
           d[0][2]*m[0][2] + d[1][2]*m[1][2] + d[2][2]*m[2][2]);
}


void
Matrix3x3::dumpto3x3(double (*dest)[3][3]) const
{

  (*dest)[0][0] = d[0][0];
  (*dest)[0][1] = d[0][1];
  (*dest)[0][2] = d[0][2];
  (*dest)[1][0] = d[1][0];
  (*dest)[1][1] = d[1][1];
  (*dest)[1][2] = d[1][2];
  (*dest)[2][0] = d[2][0];
  (*dest)[2][1] = d[2][1];
  (*dest)[2][2] = d[2][2];
}



void
Matrix3x3::swapCols(const int c1, const int c2)
{
  double t;

  t = d[0][c1];
  d[0][c1] = d[0][c2];
  d[0][c2] = t;
  t = d[1][c1];
  d[1][c1] = d[1][c2];
  d[1][c2] = t;
  t = d[2][c1];
  d[2][c1] = d[2][c2];
  d[2][c2] = t;
}


void
Matrix3x3::swapRows(const int r1, const int r2)
{
  double t;

  t = d[r1][0];
  d[r1][0] = d[r2][0];
  d[r2][0] = t;
  t = d[r1][1];
  d[r1][1] = d[r2][1];
  d[r2][1] = t;
  t = d[r1][2];
  d[r1][2] = d[r2][2];
  d[r2][2] = t;
}


double
Matrix3x3::normFrobenius(const int skipSqrt) const
{
  // Frobenius norm
  double norm = ZERO;
  norm += d[0][0] * d[0][0];
  norm += d[0][1] * d[0][1];
  norm += d[0][2] * d[0][2];
  norm += d[1][0] * d[1][0];
  norm += d[1][1] * d[1][1];
  norm += d[1][2] * d[1][2];
  norm += d[2][0] * d[2][0];
  norm += d[2][1] * d[2][1];
  norm += d[2][2] * d[2][2];
  if(skipSqrt)
    return norm;
  return sqrt(norm);
}


double
Matrix3x3::normOff(const int skipSqrt) const
{
  // the "norm" of the off diagonal elements.
  double norm = ZERO;
  norm += d[0][1] * d[0][1];
  norm += d[0][2] * d[0][2];
  norm += d[1][0] * d[1][0];
  norm += d[1][2] * d[1][2];
  norm += d[2][0] * d[2][0];
  norm += d[2][1] * d[2][1];

  if(skipSqrt)
    return norm;
  return sqrt(norm);
}





inline int
symSchur(const double d[3][3], const int p, const int q,
	 double &c, double &s)
{
  if(d[p][q] != ZERO)
    {
      double tau = (d[q][q] - d[p][p])/(TWO * d[p][q]);
      double t = ((double)signOf(tau))/(myfabs(tau)+sqrt(ONE + tau*tau));
      c = ONE/sqrt(1+t*t); s = t*c;
      return 1;
    }
  else
    {
      c = ONE; s = ZERO;
      return 0;
    }
}


unsigned int
Matrix3x3::symEigen(Matrix3x3 &eigVec, Vector3 &eigVal, const double tol) const
{
  // Returns the eigenvalues and eigenvectors of a symmetric matrix

  // symmetry test
  assert(d[0][1] == d[1][0]);
  assert(d[0][2] == d[2][0]);
  assert(d[1][2] == d[2][1]);

  const double eps = tol * normFrobenius();
  unsigned int iterCount = 0;

  // We are using a classical Jacobi method (Golub & van Loan, sec
  // 8.5) Yeah, I know that is overkill for a 3x3, but it is fairly
  // straight-forward to code.  Note that we may do more rotations
  // than are absolutely necessary because for each pass of the while
  // loop, we do three rotations, though the off(V'AV)<=eps condition
  // may have been satisfied by only one or two of them.  This only
  // improves accuracy.


  Matrix3x3 A = *this;
  Matrix3x3 J;			// Jacobi/Givens rotation

  eigVec.set2Identity();
  J.set2Identity();
  double c, s;
  int p,q;
  while(A.normOff() > eps)
    {
      for(p = 0; p < 2; ++p)
	for(q = p+1; q < 3; ++q)
	  {
	    if(symSchur(A.d, p, q, c, s))
	      {
		J[p][p] = c;
		J[q][q] = c;
		J[p][q] = s;
		J[q][p] = -s;
		A = J.transpose_mul(A*J);
		eigVec = eigVec*J;
		J[p][p] = ONE;
		J[q][q] = ONE;
		J[p][q] = ZERO;
		J[q][p] = ZERO;
		iterCount += 1;
	      }
	  }
    }
  eigVal[0] = A[0][0];
  eigVal[1] = A[1][1];
  eigVal[2] = A[2][2];
  return iterCount;
}


int Ax_b(Matrix3x3 A, Vector3& x, Vector3 b)
{
  if(approxZero(A.det()))
    {
      // Numerically singular, implies that the q1q2 segment is
      // orthogonal to the norm of the triangle
      return 1;
    }

  // Find pivot for first row
  int q = 0;
  if(myfabs(A[1][0]) > myfabs(A[0][0]))
    q = 1;
  if(myfabs(A[2][0]) > myfabs(A[q][0]))
    q = 2;

  // Swap rows, if necessary
  if(q != 0)
    {
      A.swapRows(0,q);
      const double t = b[0];
      b[0] = b[q];
      b[q] = t;
    }

  // Scale top row
  A[0][1] /= A[0][0];
  A[0][2] /= A[0][0];
  b[0] /= A[0][0];
  // A[0][0] = ONE;

  // Eliminate lower entries
  A[1][1] -= A[1][0] * A[0][1];
  A[1][2] -= A[1][0] * A[0][2];
  b[1] -= A[1][0] * b[0];
  // A[1][0] = ZERO;
  A[2][1] -= A[2][0] * A[0][1];
  A[2][2] -= A[2][0] * A[0][2];
  b[2] -= A[2][0] * b[0];
  // A[2][0] = ZERO;

  // Find pivot for second row
  if(myfabs(A[2][1]) > myfabs(A[1][1]))
    {
      A.swapRows(1,2);
      const double t = b[1];
      b[1] = b[2];
      b[2] = t;
    }

  // Scale middle row
  A[1][2] /= A[1][1];
  b[1] /= A[1][1];
  // A[1][1] = ONE;

  // Eliminate lower entries
  A[2][2] -= A[2][1] * A[1][2];
  b[2] -= A[2][1] * b[1];
  // A[2][1] = ZERO;

  // Backsub
  x[2] = b[2]/A[2][2];
  x[1] = b[1] - A[1][2] * x[2];
  x[0] = b[0] - A[0][2] * x[2] - A[0][1] * x[1];
  return 0;
}


void Matrix3x3::makeQuat(double &q0, double &q1, double &q2, double &q3) const
{
  const double tr = trace();
  double a;
  bool ok = false;
  q0 = sqrt((tr + ONE)/FOUR);

  if(!approxZero(q0))
    {
      // q0 is not zero, do the easy case;
      a = FOUR * q0;
      q1 = (d[2][1] - d[1][2])/a;
      q2 = (d[0][2] - d[2][0])/a;
      q3 = (d[1][0] - d[0][1])/a;
      ok = true;
    }
  else
    {
      q0 = ZERO;
      a = (ONE - tr)/FOUR;
      q1 = sqrt(a + d[0][0]/TWO);
      if(!approxZero(q1))
	{
	  q2 = (d[1][0] + d[0][1])/(FOUR * q1);
	  q3 = (d[2][0] + d[0][2])/(FOUR * q1);
	  ok = true;
	}
      else
	{
	  q2 = sqrt(a + d[1][1]/TWO);
	  if(!approxZero(q2))
	    {
	      q1 = (d[1][0] + d[0][1])/(FOUR * q2);
	      q3 = (d[2][1] + d[1][2])/(FOUR * q2);
	      ok = true;
	    }
	  else
	    {
	      q3 = sqrt(a + d[2][2]/TWO);
	      if(!approxZero(q3))
		{
		  q1 = (d[2][0] + d[0][2])/(FOUR * q3);
		  q2 = (d[2][1] + d[1][2])/(FOUR * q3);
		  ok = true;
		}
	    }

	}
    }
  assert(ok == true);
  a = sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
  q0 /= a;
  q1 /= a;
  q2 /= a;
  q3 /= a;
}


void Matrix3x3::setFromQuat(const double q0, const double q1,
			    const double q2, const double q3)
{
  double q02, q0sqr2, q0q12, q0q22, q0q32;

  if(q0 == ONE)
    {
      this->set2Identity();
      return;
    }

  if(q0 == ZERO)
    {
      q02 = ZERO;
      q0sqr2 = ZERO;
      q0q12 = ZERO;
      q0q22 = ZERO;
      q0q32 = ZERO;
    }
  else
    {
      q02 = TWO*q0;
      q0sqr2 = q02*q0;
      q0q12 = (q1 == ZERO) ? ZERO : q02*q1;
      q0q22 = (q2 == ZERO) ? ZERO : q02*q2;
      q0q32 = (q3 == ZERO) ? ZERO : q02*q3;
    }

  const double q1q22 = TWO*q1*q2;
  const double q1q32 = TWO*q1*q3;
  const double q2q32 = TWO*q2*q3;


  d[0][0] = q0sqr2 + q1*q1*TWO - ONE;
  d[0][1] = q1q22 - q0q32;
  d[0][2] = q1q32 + q0q22;

  d[1][0] = q1q22 + q0q32;
  d[1][1] = q0sqr2 + q2*q2*TWO - ONE;
  d[1][2] = q2q32 - q0q12;

  d[2][0] = q1q32 - q0q22;
  d[2][1] = q2q32 + q0q12;
  d[2][2] = q0sqr2 + q3*q3*TWO - ONE;

}

#if 0
Matrix3x3&
Matrix3x3::operator=(const QuatR4& q)
{
  double q02, q0sqr2, q0q12, q0q22, q0q32;

  if(q[QuatR4::QT] == ONE)
    {
      this->set2Identity();
      return *this;
    }

  if(q[QuatR4::QT] == ZERO)
    {
      q02 = ZERO;
      q0sqr2 = ZERO;
      q0q12 = ZERO;
      q0q22 = ZERO;
      q0q32 = ZERO;
    }
  else
    {
      q02 = TWO*q[QuatR4::QT];
      q0sqr2 = q02*q[QuatR4::QT];
      q0q12 = (q[QuatR4::QX] == ZERO) ? ZERO : q02*q[QuatR4::QX];
      q0q22 = (q[QuatR4::QY] == ZERO) ? ZERO : q02*q[QuatR4::QY];
      q0q32 = (q[QuatR4::QZ] == ZERO) ? ZERO : q02*q[QuatR4::QZ];
    }

  const double q1q22 = TWO*q[QuatR4::QX]*q[QuatR4::QY];
  const double q1q32 = TWO*q[QuatR4::QX]*q[QuatR4::QZ];
  const double q2q32 = TWO*q[QuatR4::QY]*q[QuatR4::QZ];


  d[0][0] = q0sqr2 + q[QuatR4::QX]*q[QuatR4::QX]*TWO - ONE;
  d[0][1] = q1q22 - q0q32;
  d[0][2] = q1q32 + q0q22;

  d[1][0] = q1q22 + q0q32;
  d[1][1] = q0sqr2 + q[QuatR4::QY]*q[QuatR4::QY]*TWO - ONE;
  d[1][2] = q2q32 - q0q12;

  d[2][0] = q1q32 - q0q22;
  d[2][1] = q2q32 + q0q12;
  d[2][2] = q0sqr2 + q[QuatR4::QZ]*q[QuatR4::QZ]*TWO - ONE;

  return *this;

}
#endif

// ------------------------------------------------------------------------
// Matrix4x4 Stuff
// ------------------------------------------------------------------------

ostream
&operator <<(ostream& os, const Matrix4x4& m)
{
  os << "[[" << m[0][0] << ' ' << m[0][1] << ' ' << m[0][2] << ' ' << m[0][3] << "];...\n"
     << " [" << m[1][0] << ' ' << m[1][1] << ' ' << m[1][2] << ' ' << m[1][3] << "];...\n"
     << " [" << m[2][0] << ' ' << m[2][1] << ' ' << m[2][2] << ' ' << m[2][3] << "];...\n"
     << " [" << m[3][0] << ' ' << m[3][1] << ' ' << m[3][2] << ' ' << m[3][3] << "]]";
  return os;
}

istream
&operator >>(istream& is, Matrix4x4& v)
{
  is >> v[0][0] >> v[0][1] >> v[0][2] >> v[0][3]
     >> v[1][0] >> v[1][1] >> v[1][2] >> v[1][3]
     >> v[2][0] >> v[2][1] >> v[2][2] >> v[2][3]
     >> v[3][0] >> v[3][1] >> v[3][2] >> v[3][3];
  return is;
}

const Matrix4x4
Matrix4x4::operator*(const Matrix4x4& m) const
{
  return Matrix4x4(
	   d[0][0]*m[0][0]+d[0][1]*m[1][0]+d[0][2]*m[2][0]+d[0][3]*m[3][0],
	   d[0][0]*m[0][1]+d[0][1]*m[1][1]+d[0][2]*m[2][1]+d[0][3]*m[3][1],
	   d[0][0]*m[0][2]+d[0][1]*m[1][2]+d[0][2]*m[2][2]+d[0][3]*m[3][2],
	   d[0][0]*m[0][3]+d[0][1]*m[1][3]+d[0][2]*m[2][3]+d[0][3]*m[3][3],

	   d[1][0]*m[0][0]+d[1][1]*m[1][0]+d[1][2]*m[2][0]+d[1][3]*m[3][0],
	   d[1][0]*m[0][1]+d[1][1]*m[1][1]+d[1][2]*m[2][1]+d[1][3]*m[3][1],
	   d[1][0]*m[0][2]+d[1][1]*m[1][2]+d[1][2]*m[2][2]+d[1][3]*m[3][2],
	   d[1][0]*m[0][3]+d[1][1]*m[1][3]+d[1][2]*m[2][3]+d[1][3]*m[3][3],

	   d[2][0]*m[0][0]+d[2][1]*m[1][0]+d[2][2]*m[2][0]+d[2][3]*m[3][0],
	   d[2][0]*m[0][1]+d[2][1]*m[1][1]+d[2][2]*m[2][1]+d[2][3]*m[3][1],
	   d[2][0]*m[0][2]+d[2][1]*m[1][2]+d[2][2]*m[2][2]+d[2][3]*m[3][2],
	   d[2][0]*m[0][3]+d[2][1]*m[1][3]+d[2][2]*m[2][3]+d[2][3]*m[3][3],

	   d[3][0]*m[0][0]+d[3][1]*m[1][0]+d[3][2]*m[2][0]+d[3][3]*m[3][0],
	   d[3][0]*m[0][1]+d[3][1]*m[1][1]+d[3][2]*m[2][1]+d[3][3]*m[3][1],
	   d[3][0]*m[0][2]+d[3][1]*m[1][2]+d[3][2]*m[2][2]+d[3][3]*m[3][2],
	   d[3][0]*m[0][3]+d[3][1]*m[1][3]+d[3][2]*m[2][3]+d[3][3]*m[3][3]);
}

void
Matrix4x4::dumpto4x4(double (*dest)[4][4]) const
{
  (*dest)[0][0] = d[0][0];
  (*dest)[0][1] = d[0][1];
  (*dest)[0][2] = d[0][2];
  (*dest)[0][3] = d[0][3];
  (*dest)[1][0] = d[1][0];
  (*dest)[1][1] = d[1][1];
  (*dest)[1][2] = d[1][2];
  (*dest)[1][3] = d[1][3];
  (*dest)[2][0] = d[2][0];
  (*dest)[2][1] = d[2][1];
  (*dest)[2][2] = d[2][2];
  (*dest)[2][3] = d[2][3];
  (*dest)[3][0] = d[3][0];
  (*dest)[3][1] = d[3][1];
  (*dest)[3][2] = d[3][2];
  (*dest)[3][3] = d[3][3];
}

void
Matrix4x4::swapCols(const int c1, const int c2)
{
  double t;

  t = d[0][c1];
  d[0][c1] = d[0][c2];
  d[0][c2] = t;
  t = d[1][c1];
  d[1][c1] = d[1][c2];
  d[1][c2] = t;
  t = d[2][c1];
  d[2][c1] = d[2][c2];
  d[2][c2] = t;
  t = d[3][c1];
  d[3][c1] = d[3][c2];
  d[3][c2] = t;
}


void
Matrix4x4::swapRows(const int r1, const int r2)
{
  double t;

  t = d[r1][0];
  d[r1][0] = d[r2][0];
  d[r2][0] = t;
  t = d[r1][1];
  d[r1][1] = d[r2][1];
  d[r2][1] = t;
  t = d[r1][2];
  d[r1][2] = d[r2][2];
  d[r2][2] = t;
  t = d[r1][3];
  d[r1][3] = d[r2][3];
  d[r2][3] = t;
}


double
Matrix4x4::normFrobenius(const int skipSqrt) const
{
  // Frobenius norm
  double norm = ZERO;
  norm += d[0][0] * d[0][0];
  norm += d[0][1] * d[0][1];
  norm += d[0][2] * d[0][2];
  norm += d[0][3] * d[0][3];
  norm += d[1][0] * d[1][0];
  norm += d[1][1] * d[1][1];
  norm += d[1][2] * d[1][2];
  norm += d[1][3] * d[1][3];
  norm += d[2][0] * d[2][0];
  norm += d[2][1] * d[2][1];
  norm += d[2][2] * d[2][2];
  norm += d[2][3] * d[2][3];
  norm += d[3][0] * d[3][0];
  norm += d[3][1] * d[3][1];
  norm += d[3][2] * d[3][2];
  norm += d[3][3] * d[3][3];
  if(skipSqrt)
    return norm;
  return sqrt(norm);
}

#if 0
int Ax_b(Matrix4x4 A, Vector4& x, Vector4 b)
{
  if(approxZero(A.det()))
    {
      // Numerically singular, implies that the q1q2 segment is
      // orthogonal to the norm of the triangle
      return 1;
    }

  // Find pivot for first row
  int q = 0;
  if(myfabs(A[1][0]) > myfabs(A[0][0]))
    q = 1;
  if(myfabs(A[2][0]) > myfabs(A[q][0]))
    q = 2;

  // Swap rows, if necessary
  if(q != 0)
    {
      A.swapRows(0,q);
      const double t = b[0];
      b[0] = b[q];
      b[q] = t;
    }

  // Scale top row
  A[0][1] /= A[0][0];
  A[0][2] /= A[0][0];
  b[0] /= A[0][0];
  // A[0][0] = ONE;

  // Eliminate lower entries
  A[1][1] -= A[1][0] * A[0][1];
  A[1][2] -= A[1][0] * A[0][2];
  b[1] -= A[1][0] * b[0];
  // A[1][0] = ZERO;
  A[2][1] -= A[2][0] * A[0][1];
  A[2][2] -= A[2][0] * A[0][2];
  b[2] -= A[2][0] * b[0];
  // A[2][0] = ZERO;

  // Find pivot for second row
  if(myfabs(A[2][1]) > myfabs(A[1][1]))
    {
      A.swapRows(1,2);
      const double t = b[1];
      b[1] = b[2];
      b[2] = t;
    }

  // Scale middle row
  A[1][2] /= A[1][1];
  b[1] /= A[1][1];
  // A[1][1] = ONE;

  // Eliminate lower entries
  A[2][2] -= A[2][1] * A[1][2];
  b[2] -= A[2][1] * b[1];
  // A[2][1] = ZERO;

  // Backsub
  x[2] = b[2]/A[2][2];
  x[1] = b[1] - A[1][2] * x[2];
  x[0] = b[0] - A[0][2] * x[2] - A[0][1] * x[1];
  return 0;
}
#endif
