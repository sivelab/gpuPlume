#include <math.h>
#include "Random.h"

static float randVal()
{
#ifndef WIN32
  return (float)(rand()/(float)RAND_MAX); 
#else
  return drand48(); 
#endif
}

Random::Random()
{
  // Seed the random number generator with a somewhat random (or at
  // least changing) value

#ifndef WIN32
  // On Unix systems, this should do it...
  srand48( time(0) % getpid() );
#else
  // Otherwise, just give a value until we figure out how to do it on
  // Windows.
  srand(2);
#endif
}

Random::Random(long s)
{
  // Seed the random number generator with the value passed as an
  // argument. 

#ifndef WIN32
  // On Unix systems, this should do it...
  srand48( s );
#else
  // Otherwise, just give a value until we figure out how to do it on
  // Windows.
  srand( s );
#endif
}

float Random::uniform()
{
  return randVal();
}

float Random::normal()
{
  float fac, rsq, v1, v2;

  if (m_normal_value == false)
    {
      do 
	{
	  v1 = 2.0 * randVal() - 1.0;
	  v2 = 2.0 * randVal() - 1.0;
	  rsq = v1*v1 + v2*v2;
	} while (rsq >= 1.0);
      
      rsq = sqrt( (-2.0 * log(rsq) ) / rsq );
      
      m_remaining_value = v2 * rsq;
      m_normal_value = true;

      return v1*rsq;
    }
  else
    {
      m_normal_value = false;
      return m_remaining_value;
    }
}

