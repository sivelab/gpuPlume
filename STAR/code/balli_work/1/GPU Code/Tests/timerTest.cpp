#include <iostream>
#include <unistd.h>
#include "Timer.h"

int main(int argc, char*argv[])
{
  Timer t1;
  Timer t2(true);

  Timer_t start_t1, end_t1;
  Timer_t start_t2, end_t2;

  std::cout << "Using a usleep to wait for 0.5 seconds" << std::endl;
  for (int i=0; i<2; i++)
    {
      start_t1 = t1.tic();
      
      // wait for at least N microseconds 
      // 1 second = 1000000
      // 0.5 seconds = 500000
      // 1 ms = 1000
      usleep(500000);
      
      end_t1 = t1.tic();
      std::cout << "Low Res Timer: " << t1.deltas(start_t1, end_t1) << " sec" << std::endl;
      std::cout << "\t" << t1.deltam(start_t1, end_t1) << " msec" << std::endl;
      std::cout << "\t" << t1.deltau(start_t1, end_t1) << " usec" << std::endl;
    }

  // HIGH Res Timer

  for (int i=0; i<2; i++)
    {
      start_t2 = t2.tic();

      // wait for at least N microseconds 
      // 1 second = 1000000
      // 1 ms = 1000
      usleep(500000);
      
      end_t2 = t2.tic();
      std::cout << "High Res Timer: " << t2.deltas(start_t2, end_t2) << " sec" << std::endl;
      std::cout << "\t" << t2.deltam(start_t2, end_t2) << " msec" << std::endl;
      std::cout << "\t" << t2.deltau(start_t2, end_t2) << " usec" << std::endl;
    }

}
