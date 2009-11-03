#ifndef RGB_H
#define RGB_H

/*
** RGB Image Structure
*/

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <stdlib.h>
#include <GL/glut.h>
// #include <GL/gl.h>
#endif

typedef struct _RGBImageRec {
    GLint sizeX, sizeY;
    unsigned char *data;
} RGBImageRec;

extern RGBImageRec *rgbImageLoad(const char *);

#endif // RGB_H
