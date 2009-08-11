#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Client_Controller.h"
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
#define printw5 printf
#define printwvis1 printf

/* include correct clock */
#include "vxSync.h"
#define rtcSync vxSync

#define sleep(a) { \
  struct timespec vx_r; \
  vx_r.tv_sec = a; \
  vx_r.tv_nsec = 0; \
}

//  nanosleep(&vx_r, NULL); \ - need to add timerLib into kernel

/* not supporting the read double correctly, not sure what will happen if I try */
#define getch getchar
#define addch(a) NULL

#elif defined(WIN32)
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

#define printwvis1(a) \
	WriteConsole(visiblebuf,a,strlen(a),&numwritten,NULL)

#define addch(x) { \
	char tmponecharstr[1] = {x}; \
	unsigned long numwritten; \
    WriteConsole(visiblebuf,tmponecharstr,1,&numwritten,NULL); \
}
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

#define printw5(a,b,c,d,e) { \
		strncpy(tmpstring,blankstring,128); \
		sprintf(tmpstring,a,b,c,d,e); \
		WriteConsole(visiblebuf,tmpstring,strlen(tmpstring),&numwritten,NULL); \
}


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
#define printw5 printw
#define printwvis1 printw
#endif


#define EPS .0001

/* function declarations */
int TRGetUserState();
void RegisterPackets(int type = 0);
int TRDisable();
int TREnable();
int TRSetPosition(const UserState_Data& usd, const VirtualState_Data& vsd);
int TRSetVirtualState(VirtualState_Data& vsd);
int TRGetAllServerInfo();
int TRSetAllVirtInfo(VirtualState_Data& vsd, Collision_Data& cd, Stairs_Data& sd);
int TRSendByteArray(ByteArray_Data& cad);

// Functions for handling input, windows, etc. - works for any
//  OS - linux, Win, VxWorks
int GetCharacter();
void InitWindow();
void EraseWindow();
void StopWindow();
double mycurs_getdouble(void);
