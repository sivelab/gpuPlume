#include <iostream>
#include <math.h>
#include <assert.h>
#include <vector>
#include "Random.h"

int main(int argc, char *argv[])
{
  if (argc != 2) exit(1);
  int N = atoi(argv[1]);

  Random random(2);

  std::cout << "UNIFORM RANDOM" << std::endl;
  for (int i=0; i<N; i++)
    std::cout << "Random: " << Random::uniform() << std::endl;

  std::vector<float> random_values(N);
  std::cout << "\nNORMAL RANDOM" << std::endl;
  float v, avg = 0;
  for (int i=0; i<N; i++)
    {
      random_values[i] = Random::normal();
      avg += random_values[i];
      std::cout << "Random: " << random_values[i] << std::endl;
    }
  float mean = avg/(float)N;
  std::cout << "Mean: " << mean << std::endl;

  float isum = 0.0;
  for (int i=0; i<N; i++)
    {
      isum += ((random_values[i] - mean) * (random_values[i] - mean));
    }
  float sigma = sqrt(1.0/(float)N * isum);
  std::cout << "Standard Deviation: " << sigma << ", Variance = " << sigma*sigma << std::endl;
  
  std::cout << "-----------------------------" << std::endl;
  
#if 0
  int twidth = 256, theight = 256;
  float f1 = Random::uniform() * twidth;
  float f2 = Random::uniform() * theight;

  std::cout << "Setting random offset: " << f1 << ", " << f2 << std::endl;

  while (1) {
  float r1, r2;
  for (int i=0; i<20; i++)
    {
      r1 = Random::uniform() * twidth;
      r2 = Random::uniform() * theight;
      
      r1 = r1 + f1;
      r2 = r2 + f2;

      // std::cout << "Offset value: " << r1 << ", " << r2 << std::endl;

      if (r1 > twidth) r1 -= twidth;
      if (r2 > theight) r2 -= theight;
      assert( r1 <= twidth );
      assert( r2 <= theight );
      assert( r1 >= 0.0 && r2 >= 0.0 );
    }
  }
#endif
}
