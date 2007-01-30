#include <stdio.h>
#include <string.h>

#include "Timer.h"

#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>

#ifdef WIN32

#include <sys/types.h>
#include <fcntl.h>
#include <windows.h>
#include <winbase.h>

Timer::Timer()
{
  LARGE_INTEGER frequency;
  if(QueryPerformanceFrequency(&frequency))
    {
      _secsPerTic = 1.0/(double)frequency.QuadPart;
    }
  else
    {
      _secsPerTic = 1.0;
      notify(NOTICE)<<"Error: Timer::Timer() unable to use QueryPerformanceFrequency, "<<std::endl;
      notify(NOTICE)<<"timing code will be wrong, Windows error code: "<<GetLastError()<<std::endl;
    }
}

Timer_t Timer::tic() const
{
  LARGE_INTEGER qpc;
  if (QueryPerformanceCounter(&qpc))
    {
      return qpc.QuadPart;
    }
  else
    {
      notify(NOTICE)<<"Error: Timer::Timer() unable to use QueryPerformanceCounter, "<<std::endl;
      notify(NOTICE)<<"timing code will be wrong, Windows error code: "<<GetLastError()<<std::endl;
      return 0;
    }
}

#else

#include <sys/time.h>

Timer::Timer( void )
{
  _secsPerTic = (1.0 / (double) 1000000);
}

Timer_t Timer::tic() const
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return ((Timer_t)tv.tv_sec)*1000000+(Timer_t)tv.tv_usec;
}

#endif
