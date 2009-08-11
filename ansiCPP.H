/****************************************************************************
 * ISO C++ Approximation functions
 *
 * Paul Henning
 * The University of Iowa Computer Science Department
 * phenning@cs.uiowa.edu
 *
 * $Id: ansiCPP.H,v 1.1.1.1 1999/06/10 18:48:05 willemsn Exp $
 ****************************************************************************/

#ifndef ANSICPP_HEADER
#define ANSICPP_HEADER 1



#ifdef NO_CPP_STATIC_CAST
#define STATIC_CAST(T,V) (T)(V)
#else
#define STATIC_CAST(T,V) static_cast<T>(V)
#endif

#ifdef NO_CPP_CONST_CAST
#define CONST_CAST(T,V) (T)(V)
#else
#define CONST_CAST(T,V) const_cast<T>(V)
#endif


#ifdef NO_CPP_BOOL
typedef int bool;
const int true = 1;
const int false = 0;
#endif

#endif
