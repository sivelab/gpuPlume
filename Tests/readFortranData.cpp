#include <iostream>
#include <fstream>

int main(int argc, char *argv[])
{
  std::ifstream fortranInput("QU_velocity.bin", std::ios::binary);

  fortranInput.seekg(0, std::ios::end);
  int numBytes = fortranInput.tellg();
  fortranInput.seekg(0, std::ios::beg);

  // allocate memory for whole file and read everything in at once
  char *fileBuffer = new char[numBytes];
  fortranInput.read(fileBuffer, numBytes);
  fortranInput.close();

  int cd = 0;

  // skip over 4 or 8 bytes since that's how fortran dumps binary
  // writes (depends on 32-bit versus 64-bit, respectively and I'm not
  // sure how to determine that yet)

  cd += sizeof(int);

  int nx, ny, nz;
  nx = ny = nz = 2;

  double *u = new double[nx * ny * nz];
  double *v = new double[nx * ny * nz];
  double *w = new double[nx * ny * nz];

  double *dval;
  int dataN = nx*ny*nz;
  int dataSz = dataN * sizeof(double);

  for (int i=0; i<(dataN); i++)
    {
      dval = reinterpret_cast<double*>(&fileBuffer[cd]);
      u[i] = *dval;

      dval = reinterpret_cast<double*>(&fileBuffer[cd + dataSz]);
      v[i] = *dval;

      dval = reinterpret_cast<double*>(&fileBuffer[cd + dataSz*2]);
      w[i] = *dval;

      cd += sizeof(double);
    }

  for (int i=0; i<(nx * ny * nz); i++)
    std::cout << "uvw(" << i << ") = " << u[i] << "  " << v[i] << "  " << w[i] << std::endl;

  delete [] fileBuffer;
}

