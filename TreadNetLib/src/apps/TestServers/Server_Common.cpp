#define WIN_NOVARINCLUDE
#include "Server_header.h"

/* Variable declarations */
char *statstr[64] = {
"INVALID!!!!!0",
"-- ENABLED BIT",
"-- DISABLED BIT",
"INVALID!!!!!3",
"-- ACTIVE_BIT",
"ACTIVE_ENABLED",
"ACTIVE_DISABLED",
"INVALID!!!!!7",
"-- INACTIVE BIT",
"INACTIVE_ENABLED",
"INACTIVE_DISABLED",
"INVALID!!!!!11",
"INVALID!!!!!12",
"INVALID!!!!!13",
"INVALID!!!!!14",
"INVALID!!!!!15",
"-- RAMPING_UP_BIT",
"INVALID!!!!!17",
"INVALID!!!!!18",
"INVALID!!!!!19",
"INVALID!!!!!20",
"ACTIVE_ENABLED RAMPING UP",
"INVALID!!!!!22",
"INVALID!!!!!23",
"INVALID!!!!!24",
"INVALID!!!!!25",
"INVALID!!!!!26",
"INVALID!!!!!27",
"INVALID!!!!!28",
"INVALID!!!!!29",
"INVALID!!!!!30",
"INVALID!!!!!31",
"INVALID!!!!!32",
"INVALID!!!!!33",
"INVALID!!!!!34",
"INVALID!!!!!35",
"INVALID!!!!!36",
"ACTIVE_ENABLED RAMPING DOWN",
"INVALID!!!!!38",
"INVALID!!!!!39",
"INVALID!!!!!40",
"INVALID!!!!!41",
"INVALID!!!!!42",
"INVALID!!!!!43",
"INVALID!!!!!44",
"INVALID!!!!!45",
"INVALID!!!!!46",
"INVALID!!!!!47",
"INVALID!!!!!48",
"INVALID!!!!!49",
"INVALID!!!!!50",
"INVALID!!!!!51",
"INVALID!!!!!52",
"INVALID!!!!!53",
"INVALID!!!!!54",
"INVALID!!!!!55",
"INVALID!!!!!56",
"INVALID!!!!!57",
"INVALID!!!!!58",
"INVALID!!!!!59",
"INVALID!!!!!60",
"INVALID!!!!!61",
"INVALID!!!!!62",
"INVALID!!!!!63"
};

char *terrstr[4] = {"NORMAL","STAIRS","ROUGH","ICY"};

#ifdef RTI_VXWORKS
char *speedkey_string = "t/g          ";
char *turnkey_string  = "f/h          ";
#else
char *speedkey_string = "Arrow Up/Down";
char *turnkey_string  = "Arrow R/L    ";
#endif

char *stairs_dir_string[3] = {"IGNORE", "UP", "DOWN"};
char *stairs_cmd_string[4] = {"NONE", "SWITCH", "INC FACT", "DEC FACT"};

// server controller
Server_Controller *server;

// user state packets
UserState_Packet * UsrState_pkt;

// virtual state packet
VirtualState_Packet * VirtState_pkt;
int VirtState_type;

// Command packet
Command_Packet * Cmd_pkt;
int Cmd_type;

// New position packets
NewPos_Packet * NewPos_pkt;
NewPosResponse_Packet * NewPosResp_pkt;
int NewPos_type;

// Collision packet
Collision_Packet * Coll_pkt;
int Coll_type;

// FSRSignals packet
FSRSignals_Packet * Fsr_pkt;

// Stairs packet
Stairs_Packet * Stairs_pkt;
int Stairs_type;

// Tracker packet
Tracker_Packet * Tracker_pkt;

// Character array packet
ByteArray_Packet * ByteArray_pkt;
int ByteArray_type;

// Server Request packet
ServerRequest_Packet * ServerRequest_pkt;
int ServerRequest_type;

#if defined(WIN32)
// include windows variable
extern HANDLE visiblebuf;
#endif

void RegisterPackets() {
  // Register packets we are going to send/accept
  // Create all packets
  UsrState_pkt = new UserState_Packet();
  VirtState_pkt = new VirtualState_Packet();
  Cmd_pkt = new Command_Packet();
  NewPos_pkt = new NewPos_Packet();
  NewPosResp_pkt = new NewPosResponse_Packet();
  Coll_pkt = new Collision_Packet();
  Fsr_pkt = new FSRSignals_Packet();
  Stairs_pkt = new Stairs_Packet();
  Tracker_pkt = new Tracker_Packet();
  ServerRequest_pkt = new ServerRequest_Packet();
  ByteArray_pkt = new ByteArray_Packet();
  
  // Request and response for User State
  int result = server->Register(UsrState_pkt);
  if (result < 1) {
    // error in Register
    server->PrintLastError();
    exit(1);
  }

  ServerRequest_type = server->Register(ServerRequest_pkt);
  if (ServerRequest_type < 1) {
    // error in Register
    server->PrintLastError();
    exit(1);
  }

  // Virtual state
  VirtState_type = server->Register(VirtState_pkt);
  if (VirtState_type < 1) {
    server->PrintLastError();
    exit(1);
  }

  // Command packet
  Cmd_type = server->Register(Cmd_pkt);
  if (Cmd_type < 1) {
    server->PrintLastError();
    exit(1);
  }

  // New Position request and response
  NewPos_type = server->Register(NewPos_pkt);
  if (NewPos_type < 1) {
    server->PrintLastError();
    exit(1);
  }

  result = server->Register(NewPosResp_pkt);
  if (result < 1) {
    server->PrintLastError();
    exit(1);
  }

  // Collision
  Coll_type = server->Register(Coll_pkt);
  if (Coll_type < 1) {
    server->PrintLastError();
    exit(1);
  }

  // Stairs
  Stairs_type = server->Register(Stairs_pkt);
  if (Stairs_type < 1) {
    server->PrintLastError();
    exit(1);
  }

  // Tracker
  result = server->Register(Tracker_pkt);
  if (result < 1) {
    server->PrintLastError();
    exit(1);
  }

  // FSR Signals
  result = server->Register(Fsr_pkt);
  if (result < 1) {
    server->PrintLastError();
    exit(1);
  }

  // Byte array
  ByteArray_type = server->Register(ByteArray_pkt);
  if (result < 1) {
    server->PrintLastError();
    exit(1);
  }
}

int GetCharacter() {
#if defined(RTI_VXWORKS)
  fd_set read_fds;
  struct timeval poll;
  poll.tv_sec = 0;
  poll.tv_usec = 0;

  // setup stdin correctly - not wait for line, read char by char
  int options = ioctl(fileno(stdin), FIOGETOPTIONS,0);
  options &= ~(OPT_LINE);
  ioctl(fileno(stdin), FIOSETOPTIONS, options);

  // Need non-blocking, this was the only way I could figure to do this
  FD_ZERO(&read_fds);
  FD_SET(0,&read_fds);
  int ch = 0;

  select(1, &read_fds, (fd_set *) 0, (fd_set *) 0, &poll);
  if (FD_ISSET(0, &read_fds)) {
    // stdin is ready to read
    ch = getchar();
  }
  return ch;
  
#elif defined(_WIN32)
  int ch = 0;
  if (_kbhit()) {
    ch = getch();
    // some keys are handled by returning 0 and then the key code
    if (ch == 0 || ch == 0xE0)
      ch = getch() + 255;
  }
  return ch;
 
#else // Linux
  return getch();
#endif
}

void InitWindow() {
#if defined(RTI_VXWORKS)
  // nothing to do since no curses, or windows
#elif defined(_WIN32)
  visiblebuf = CreateConsoleScreenBuffer( GENERIC_READ | GENERIC_WRITE,
					  FILE_SHARE_READ | FILE_SHARE_WRITE,
					  NULL, CONSOLE_TEXTMODE_BUFFER, NULL);
  SetConsoleActiveScreenBuffer(visiblebuf);
  //  invisiblebuf = CreateConsoleScreenBuffer( GENERIC_READ | GENERIC_WRITE,
  //					    FILE_SHARE_READ | FILE_SHARE_WRITE,
  //					    NULL, CONSOLE_TEXTMODE_BUFFER, NULL);  
#else // Linux
  /* start curses up */
  initscr();
  cbreak();
  noecho();

  nonl();
  intrflush(stdscr,FALSE);
  keypad(stdscr,TRUE);
  nodelay(stdscr,TRUE);
#endif
}

void EraseWindow() {
#if defined(RTI_VXWORKS)
  // nothing
#elif defined(_WIN32)
  COORD  cursorpos;
  unsigned long numwritten;
  cursorpos.X = 0;
  cursorpos.Y = 0;
  
  FillConsoleOutputCharacter( 
			     visiblebuf,  
			     ' ',
			     80*25,
			     cursorpos,
			     &numwritten); 
  SetConsoleCursorPosition(visiblebuf,cursorpos); 
  
#else // Linux
  erase();
#endif
}

void StopWindow() {
#if defined(RTI_VXWORKS)
  // nothing
#elif defined(_WIN32)
  // nothing either
#else // Linux
  endwin();
#endif
}



/* These functions convert from a 2d unit
   heading vector to a 3d unit heading vector

   I've decided that the treadport simply keeps
   track of the user's two-dimensional heading
   (i.e. facing) and has to use the slope info 
   to calculate the 3d velocity vector and the
   actual change in the east/north position. 

   This function describes the conversion */

/* WARNING: this one normalizes the vector,
   and so only works for unit vectors. */
TRVector convert2dheadingto3d(TRVector surfnormal, TRVector v) {
  v.Z() = -(surfnormal.X()*v.X() + surfnormal.Y()*v.Y())/surfnormal.Z();
  v /= v.Norm();
  return v;
}





/* These are a few depricated functions.
   I was trying to set up the treadport simulator
   so that it would keep track of the user's heading
   on some flat surface, and then just map that flat
   surface onto the current virtual slope to determine
   his true east/north heading. 

   This leads to strange effects when the user crosses
   a slope discontinuity, so I've decided on the simpler
   approach described above */

/*  
TRVector surf_to_vw(TRVector n, TRVector sf) {
  TRVector vf;

  double nx = n.X;
  double ny = n.Y;
  double nz = n.Z;
  double recip = 1.0/(1.0+nz);

  vf.X = (ny*ny*recip+nz)*sf.X + (-nx*ny*recip)*sf.Y + nx*sf.Z;
  vf.Y = (-nx*ny*recip)*sf.X + (nx*nx*recip+nz)*sf.Y + ny*sf.Z;
  vf.Z = -nx*sf.X +            -ny*sf.Y        +        nz*sf.Z;

  return vf;
}


  
TRVector vw_to_surf(TRVector n, TRVector vf) {
  TRVector sf;

  double nx = n.X;
  double ny = n.Y;
  double nz = n.Z;
  double recip = 1.0/(1.0+nz);

  sf.X = (ny*ny*recip+nz)*vf.X + (-nx*ny*recip)*vf.Y - nx*vf.Z;
  sf.Y = (-nx*ny*recip)*vf.X + (nx*nx*recip+nz)*vf.Y - ny*vf.Z;
  sf.Z = nx*vf.X +             ny*vf.Y        +       nz*vf.Z;

  return sf;
}
*/

