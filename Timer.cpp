#include <iostream>

#include "Timer.h"

#if 0
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#endif

#ifdef WIN32

#include <sys/types.h>
#include <fcntl.h>
#include <windows.h>
#include <winbase.h>

Timer::Timer( bool enable_high_res_timer )
  : _use_high_res_timer(false)
{
  LARGE_INTEGER frequency;
  if(QueryPerformanceFrequency(&frequency))
    {
      _secsPerTic = 1.0/(double)frequency.QuadPart;
    }
  else
    {
      _secsPerTic = 1.0;
      std::cerr << "Error: Timer::Timer() unable to use QueryPerformanceFrequency, " << std::endl;
      std::cerr << "timing code will be wrong, Windows error code: "<< GetLastError() << std::endl;
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
      std::cerr << "Error: Timer::Timer() unable to use QueryPerformanceCounter, " << std::endl;
      std::cerr << "timing code will be wrong, Windows error code: "<< GetLastError() << std::endl;
      return 0;
    }
}

#else

#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/time.h>

// This is the hardware assembly code for the high resolution, low
// latency timer on x86 machines.  This timer will only work on x86
// hardware running Linux and must be enabled as a default argument to
// the constructor.
#define CLK(x)      __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x))

Timer::Timer( bool enable_high_res_timer )
  : _use_high_res_timer(enable_high_res_timer)
{
  if (_use_high_res_timer)
    {
      char buff[128];
      FILE *fp = fopen( "/proc/cpuinfo", "r" );
      
      double cpu_mhz=0.0f;
      while( fgets( buff, sizeof( buff ), fp ) > 0 )
	{
	  if( !strncmp( buff, "cpu MHz", strlen( "cpu MHz" )))
	    {
	      char *ptr = buff;
	      
	      while( ptr && *ptr != ':' ) ptr++;
	      if( ptr ) 
		{
		  ptr++;
		  sscanf( ptr, "%lf", &cpu_mhz );
		}
	      break;
	    }
	}
      fclose( fp );
      
      if (cpu_mhz==0.0f)
	{
	  // error - no cpu_mhz found, guess the secs per tic...
	  Timer_t start_time = tic();
	  sleep (1);
	  Timer_t end_time = tic();
	  _secsPerTic = 1.0/(double)(end_time-start_time);
	}
      else
	{
	  _secsPerTic = 1e-6/cpu_mhz;
	}
    }
  else 
    {
      // use standard, gettimeofday timing mechanism
      _secsPerTic = (1.0 / (double) 1000000);
    }
}

Timer_t Timer::tic() const
{
  if (_use_high_res_timer)
    {
      Timer_t x;CLK(x);return x;
    }
  else
    {
      struct timeval tv;
      gettimeofday(&tv, NULL);
      return ((Timer_t)tv.tv_sec)*1000000+(Timer_t)tv.tv_usec;
    }
}

#endif // else clause for ifdef WIN32
