//
// This code is from the rtfs library and The Timer class in
// OpenSceneGraph, both written by Robert Osfield of OpenSceneGraph.
// "rtfs" is a real-time frame scheduler.
//
// It has been modified to provide functionality for our applications.
// -Pete Willemsen
// 
//C++ header - Open Scene Graph - Copyright (C) 1998-2001 Robert Osfield
//Distributed under the terms of the GNU Library General Public License (LGPL)
//as published by the Free Software Foundation.

#ifndef TIMER_H
#define TIMER_H 1

#if defined(_MSC_VER)
typedef __int64 Timer_t;
#else
typedef unsigned long long Timer_t;
#endif

class Timer {
  
 public:
    
  Timer();
  ~Timer() {}
    
  Timer_t tic() const;
    
  inline double deltas( Timer_t t1, Timer_t t2 ) const { return (double)(t2 - t1)*_secsPerTic; }
  inline Timer_t deltam( Timer_t t1, Timer_t t2 ) const { return Timer_t(deltas(t1,t2)*1e3); }
  inline Timer_t deltau( Timer_t t1, Timer_t t2 ) const { return Timer_t(deltas(t1,t2)*1e6); }
  inline Timer_t deltan( Timer_t t1, Timer_t t2 ) const { return Timer_t(deltas(t1,t2)*1e9); }
  
  double getSecsPerClick() { return _secsPerTic; }
  
  private :
      double                  _secsPerTic;
      bool                    _useStandardClock;
};

#endif
