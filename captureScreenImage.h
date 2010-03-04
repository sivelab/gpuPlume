#ifdef __CAPTURESCREENIMAGE_H__
#define __CAPTURESCREENIMAGE_H__ 1

#include <cstdio>
#include <string>

#include <setjmp.h>     /* for jmpbuf declaration in writepng.h */
#include "writepng.h"

class CaptureScreenImage
  {
  public:
    PNGImage() {};

    bool writeFileData(const std::string& filename, const int width, const int height, const float *data);	
    
  private:
    void cleanup(void);

    mainprog_info m_png_fileinfo;
  };

}

#endif // __PNGIMAGE_H__


#endif // #define __CAPTURESCREENIMAGE_H__
