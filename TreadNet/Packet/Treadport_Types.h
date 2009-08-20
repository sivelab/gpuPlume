#ifndef TREADPORT_TYPES_H
#define TREADPORT_TYPES_H

#include <math.h>
#include <string.h>

#if defined(RTI_VXWORKS)
// Definitions/includes needed for vxworks
#include <unistd.h>
#define uint unsigned int
#define ushort unsigned short
#define BOOL bool
#undef m_data
#undef m_type
char * strdup(const char *str);

#elif defined(WIN32)
// Definitions/includes for windows
#include <windows.h>
#define uint UINT
#define ushort USHORT

#else
// Definitions/includes for linux/solaris
#include <unistd.h>
#define BOOL bool

#ifndef FALSE
#define FALSE false
#endif

#ifndef TRUE
#define TRUE true
#endif

#endif

enum TRStatusBit {
  ENABLED_BIT = 1,
  DISABLED_BIT = 2,
  ACTIVE_BIT = 4,
  INACTIVE_BIT = 8,
  RAMPING_UP_BIT = 16,
  RAMPING_DOWN_BIT = 32,
  MAN_TILT_BIT = 64
};

enum TRStatus {
  OFF = 0,
  ACTIVE_ENABLED = ACTIVE_BIT | ENABLED_BIT,
  ACTIVE_DISABLED = ACTIVE_BIT | DISABLED_BIT,
  INACTIVE_ENABLED = INACTIVE_BIT | ENABLED_BIT,
  INACTIVE_DISABLED = INACTIVE_BIT | DISABLED_BIT,
  ACTIVE_ENABLED_RAMPING_UP = ACTIVE_ENABLED | RAMPING_UP_BIT,
  ACTIVE_ENABLED_RAMPING_DOWN = ACTIVE_ENABLED | RAMPING_DOWN_BIT,
  MANUAL_TILT_MODE = MAN_TILT_BIT
};

enum TRTerrain {
  NORMAL_TERRAIN,
  STAIRS_TERRAIN,
  ROUGH_TERRAIN,
  ICY_TERRAIN
};

TRTerrain operator ++ (TRTerrain& t, int);
TRTerrain operator ++ (TRTerrain& t);

enum TRCommand {
  NO_CMD,
  ENABLE_CMD,
  DISABLE_CMD
};


#endif
