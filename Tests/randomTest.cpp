#include <iostream>
#include "Random.h"

int main(int argc, char *argv[])
{
  if (argc != 2) exit(1);
  int N = atoi(argv[1]);

  Random random(2);

  std::cout << "UNIFORM RANDOM" << std::endl;
  for (int i=0; i<N; i++)
    std::cout << "Random: " << Random::uniform() << std::endl;

  std::cout << "\nNORMAL RANDOM" << std::endl;
  float v, avg = 0;
  for (int i=0; i<N; i++)
    {
      v = Random::normal();
      avg += v;
      std::cout << "Random: " << v << std::endl;
    }
  std::cout << "Mean: " << avg/N << std::endl;
}
