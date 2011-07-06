#include <iostream>
#include <fstream>

int main(int argc, char *argv[])
{
  std::ifstream fortranInput("QP_turbfield.bin", std::ios::binary);

  fortranInput.seekg(0, std::ios::end);
  int numBytes = fortranInput.tellg();
  fortranInput.seekg(0, std::ios::beg);

  // allocate memory for whole file and read everything in at once
  char *fileBuffer = new char[numBytes];
  fortranInput.read(fileBuffer, numBytes);
  fortranInput.close();

  int cd = 0;


  int nx, ny, nz;
  nx = ny = nz = 2;

  // double *u = new double[nx * ny * nz];
  // double *v = new double[nx * ny * nz];
  // double *w = new double[nx * ny * nz];

  short *ival;
  float *dval;
  int dataN = nx*ny*nz;
  int dataSz = dataN * sizeof(double);



  while (cd < numBytes)
    {
      // write(63)sigu_c,sigv_c,sigw_c,ustarij(i,j,k),eps,elz(i,j,k)

      // skip over 4 or 8 bytes since that's how fortran dumps binary
      // writes (depends on 32-bit versus 64-bit, respectively and I'm not
      // sure how to determine that yet)

      // HEADER: so, extract the header values ... two 4 byte values...
      ival = reinterpret_cast<short*>(&fileBuffer[cd]);
      cd += sizeof(short);
      std::cout << *ival << ' ';

      ival = reinterpret_cast<short*>(&fileBuffer[cd]);
      cd += sizeof(short);
      std::cout << *ival << ' ';

      // SIGU_C
      dval = reinterpret_cast<float*>(&fileBuffer[cd]);
      cd += sizeof(float);
      std::cout << *dval << ' ';

      // SIGV_C
      dval = reinterpret_cast<float*>(&fileBuffer[cd]);
      cd += sizeof(float);
      std::cout << *dval << ' ';

      // SIGW_C
      dval = reinterpret_cast<float*>(&fileBuffer[cd]);
      cd += sizeof(float);
      std::cout << *dval << ' ';

      // USTAR_{ij}(i,j,k)
      dval = reinterpret_cast<float*>(&fileBuffer[cd]);
      cd += sizeof(float);
      std::cout << *dval << ' ';

      // EPS
      dval = reinterpret_cast<float*>(&fileBuffer[cd]);
      cd += sizeof(float);
      std::cout << *dval << ' ';

      // ELZ(i,j,k)
      dval = reinterpret_cast<float*>(&fileBuffer[cd]);
      cd += sizeof(float);
      std::cout << *dval << ' ';

      // TRAILER
      ival = reinterpret_cast<short*>(&fileBuffer[cd]);
      cd += sizeof(short);
      std::cout << *ival << ' ';

      ival = reinterpret_cast<short*>(&fileBuffer[cd]);
      cd += sizeof(short);
      std::cout << *ival << std::endl;
    }

  delete [] fileBuffer;
}

