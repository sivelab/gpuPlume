/****************************************************************************
 * Approximate Floating Point Math Functions
 *
 * Paul Henning
 * The University of Iowa Computer Science Department
 * phenning@cs.uiowa.edu
 *
 * $Id: approxMath.H,v 1.3 2000/06/28 16:11:16 willemsn Exp $
 ****************************************************************************/


#ifndef APPROXMATH_HEADER
#define APPROXMATH_HEADER

#include <math.h>

const double ZERO = 0.00000000000000000;
const double ONE = 1.00000000000000000;
const double TWO = 2.00000000000000000;
const double THREE = 3.00000000000000000;
const double FOUR = 4.00000000000000000;

const double pPI = 3.1415926535897932385;
const double p2PI = 6.2831853071795864770;
const double rad2deg = 57.295779513082320875;
const double deg2rad = 0.017453292519943295770;


// How close is near?
const double approxTol = 1.0e-8;

const double approx0 = -approxTol;
const double approx1 = ONE + approxTol;

namespace ApproxMath {
  double fabs(const double a);
}

inline double ApproxMath::fabs(const double a) 
{
  return (a < ZERO)?(-a):a;
}

inline double myfabs(const double a) 
{
  return (a < ZERO)?(-a):a;
}

inline int signOf(const double a)
{
  return (a < ZERO)?(-1):1;
}

inline bool approxEq(const double a, const double b)
{
  if(myfabs(a-b) < approxTol)
    return true;
  return false;
}

// template<class T>
// inline bool approxEq(const T& a, const T& b)
// {
//  if(myfabs(a-b) < approxTol)
//    return true;
//  return false;
// }


inline bool approxGt(const double a, const double b)
{
  if(a > (b - approxTol))
    return true;
  return false;
}

inline bool approxLt(const double a, const double b)
{
  if(a < (b + approxTol))
    return true;
  return false;
}

inline bool approxZero(const double a)
{
  if(a < approxTol && a > -approxTol)
    return true;
  return false;
}


inline bool approxIn01(const double a)
{
  if(a > approx0 && a < approx1)
    return true;
  return false;
}


inline bool approxIn01(const double a, const double b)
{
  if(a > approx0 && a < approx1 && b > approx0 && b < approx1)
    return true;
  return false;
}


// Truncates the parameter to 10 decimal places
inline double truncate10(double a)
{
  double temp = floor(a * 1.0e10);
  return temp/1.0e10;
}


#endif

