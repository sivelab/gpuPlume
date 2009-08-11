#ifndef RGB_H
#define RGB_H

/*
** RGB Image Structure
*/

#include "GL/gl.h"

typedef struct _RGBImageRec {
    GLint sizeX, sizeY;
    unsigned char *data;
} RGBImageRec;

extern RGBImageRec *rgbImageLoad(const char *);

#endif // RGB_H
