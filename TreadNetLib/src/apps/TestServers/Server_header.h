#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Server_Controller.h"
#include "Command_Packet.h"
#include "NewPos_Packet.h"
#include "NewPosResponse_Packet.h"
#include "ServerRequest_Packet.h"
#include "UserState_Packet.h"
#include "VirtualState_Packet.h"
#include "Collision_Packet.h"
#include "FSRSignals_Packet.h"
#include "Stairs_Packet.h"
#include "Tracker_Packet.h"
#include "ByteArray_Packet.h"

/* Added VxWorks support, though the print out is bad */

#if defined(RTI_VXWORKS)

#include "vxWorks.h"
#include <time.h>
#include <ioLib.h>

#define printw1 printf
#define printw2 printf
#define printw3 printf
#define printw4 printf

/* need to define different for vxworks, because the key presses are
   not registering.  Use t,g,f,h (increase/decrease speed,turn left/right) */
#define KEY_UP        116
#define KEY_DOWN      103
#define KEY_LEFT      102
#define KEY_RIGHT     104

/* include correct clock */
#include "vxSync.h"
#define rtcSync vxSync

#define getch getchar

#elif defined(_WIN32)
// ****** WINDOWS *****
#include <windows.h>
#include <winbase.h>
#include <conio.h>
#include <wincon.h>

// Include correct clock sync
#include "winSync.h"
#define rtcSync winSync

#ifndef WIN_NOVARINCLUDE
// messy, but it avoids double including these variables which causes
//  errors in windows (looks nicer then adding each individually to all
//  .cpp files
HANDLE visiblebuf;
unsigned long numwritten;

char blankstring[128] = 
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

char tmpstring[128];
#endif

#define sleep(x) Sleep(1000*x)
#define getch _getch

#define printw1(a) { \
		strncpy(tmpstring,blankstring,128); \
		sprintf(tmpstring,a); \
		WriteConsole(visiblebuf,tmpstring,strlen(tmpstring),&numwritten,NULL); \
}

#define printw2(a,b) { \
		strncpy(tmpstring,blankstring,128); \
		sprintf(tmpstring,a,b); \
		WriteConsole(visiblebuf,tmpstring,strlen(tmpstring),&numwritten,NULL); \
}
#define printw3(a,b,c) { \
		strncpy(tmpstring,blankstring,128); \
		sprintf(tmpstring,a,b,c); \
		WriteConsole(visiblebuf,tmpstring,strlen(tmpstring),&numwritten,NULL); \
}
#define printw4(a,b,c,d) { \
		strncpy(tmpstring,blankstring,128); \
		sprintf(tmpstring,a,b,c,d); \
		WriteConsole(visiblebuf,tmpstring,strlen(tmpstring),&numwritten,NULL); \
}

// Adding 255 for these keys so don't conflict w/ standard ascii
#define KEY_UP        72+255
#define KEY_DOWN      80+255
#define KEY_LEFT      75+255
#define KEY_RIGHT     77+255

#else
// linux/unix
#include <curses.h>
#include <time.h>
#include <unistd.h>
#include "rtcSync.h"

#define printw1 printw
#define printw2 printw
#define printw3 printw
#define printw4 printw
#endif


#define EPS .0001
#define MAX_SPEED        10.0

/* trig for 5 degree rotation */
#define COS5             0.9961947
#define SIN5             0.087155743

#define SIGN(x)     ((x>0)?(1.0):(-1.0))

/* function declarations */
void RegisterPackets();
TRVector convert2dheadingto3d(TRVector surfnormal, TRVector v);


// Functions for handling input, windows, etc. - works for any
//  OS - linux, Win, VxWorks
int GetCharacter();
void InitWindow();
void EraseWindow();
void StopWindow();
double mycurs_getdouble(void);
