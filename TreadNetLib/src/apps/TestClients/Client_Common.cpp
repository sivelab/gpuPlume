#define WIN_NOVARINCLUDE
#include "Client_header.h"

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

// client controller
Client_Controller *client;

// user state packets
UserState_Packet * UsrState_pkt;
int UsrState_type;

// virtual state packet
VirtualState_Packet * VirtState_pkt;

// Command packet
Command_Packet * Cmd_pkt;

// New position packets
NewPos_Packet * NewPos_pkt;
NewPosResponse_Packet * NewPosResp_pkt;
int NewPosResp_type;

// Collision packet
Collision_Packet * Coll_pkt;

// FSRSignals packet
FSRSignals_Packet * Fsr_pkt;
int FSR_type;

// Stairs packet
Stairs_Packet * Stairs_pkt;

// Tracker packet
Tracker_Packet * Tracker_pkt;
int Tracker_type;

// Server Request packet for all server packets
ServerRequest_Packet * ServerRequestAll_pkt;

// Server Request packet for just user state
ServerRequest_Packet * RequestUsrState_pkt;

// Character array packet for whatever
ByteArray_Packet * ByteArray_pkt;

#if defined(WIN32)
// include windows variable
extern HANDLE visiblebuf;
#endif

// returns 0 if fail, 1 if success
int TRDisable() {
  Cmd_pkt->SetDisableCmd();
  if (client->SendSinglePacket(Cmd_pkt)) {
    return 1;
  }
  // else error
  client->PrintLastError();
  return 0;
}

// returns 0 if fail, 1 if success
int TREnable() {
  Cmd_pkt->SetEnableCmd();
  if (client->SendSinglePacket(Cmd_pkt)) {
    return 1;
  }
  // else error
  client->PrintLastError();
  return 0;
}

// returns 0 if fail, 1 if success
int TRSetPosition(const UserState_Data& usd, const VirtualState_Data& vsd) {
  // Put these user and virtual state into new position and set it
  NewPos_pkt->SetUserState(usd);
  NewPos_pkt->SetVirtualState(vsd);

  if (client->SendSinglePacket(NewPos_pkt) &&
      client->ReceiveSinglePacket(NewPosResp_pkt)) {
    // new pos was sent and a response returned, check it
    if (NewPosResp_pkt->GetFlag() == NEW_POS_SUCCESS) {
      // success
      return 1;
    }
    else {
      // failure returned
      fprintf(stderr, "TRSetPosition::Error - Server returned failure in set position\n");
      return 0;
    }
  }
  // if got here there was some client error
  client->PrintLastError();
  return 0;
}

// returns 0 if fail, 1 if success
int TRSetVirtualState(VirtualState_Data& vsd) {
  VirtState_pkt->SetData(&vsd);
  if (client->SendSinglePacket(VirtState_pkt)) {
    return 1;
  }
  // else error
  client->PrintLastError();
  return 0;
}

// returns 0 if fail, 1 if success - leaves values in UsrState_pkt
int TRGetUserState() {
  if (client->SendSinglePacket(RequestUsrState_pkt) &&
      client->ReceiveSinglePacket(UsrState_pkt)) {
    // both worked correctly
    return 1;
  }
  // otherwise some error, print
  client->PrintLastError();
  return 0;
}

// returns 0 if fail, 1 if success - leaves values in UsrState_pkt, Tracker_pkt, & Fsr_pkt
int TRGetAllServerInfo() {
  if (client->SendSinglePacket(ServerRequestAll_pkt)) {
    // wait for a response
    if (client->Receive(TRUE) == 1) {
      int packet_type = client->GetPacketType();
      int pkts_received = 0;
      int cont_loop = 1;
      while (packet_type != -1 && cont_loop) {
	if (packet_type == Tracker_type) {
	  pkts_received |= 1;
	  if (client->GetPacketFromMsg(Tracker_pkt) != 1) {
	    cont_loop = 0;
	  }
	}
	else if (packet_type == FSR_type) {
	  pkts_received |= 2;
	  if (client->GetPacketFromMsg(Fsr_pkt) != 1) {
	    cont_loop = 0;
	  }
	}
	else if (packet_type == UsrState_type) {
	  pkts_received |= 4;
	  if (client->GetPacketFromMsg(UsrState_pkt) != 1) {
	    cont_loop = 0;
	  }
	}
	else {
	  // something was sent back wrong!!!
	  fprintf(stderr, "Error::GetAllServerInfo - sent back a type of %i\n", packet_type);
	  return 0;
	}
	packet_type = client->GetPacketType();
      }
      if (cont_loop) {
	// if finished loop regularly (otherwise want the client->PrintLastError())
	if (pkts_received == 7) {
	  // received all
	  return 1;
	}
	else {
	  fprintf(stderr, "Error::GetAllServerInfo - did not send back all packets - sent %i\n", 
		  pkts_received);
	  return 0;
	}	
      }
    }
  }
  // otherwise some error, print
  client->PrintLastError();
  return 0;
}

// returns 0 if fail, 1 if success
int TRSetAllVirtInfo(VirtualState_Data& vsd, Collision_Data& cd, Stairs_Data& sd) {
  VirtState_pkt->SetData(&vsd);
  if (client->AddPacket(VirtState_pkt) == 1) {
    Coll_pkt->SetData(&cd);
    if (client->AddPacket(Coll_pkt) == 1) {
      Stairs_pkt->SetData(&sd);
      if (client->AddPacket(Stairs_pkt) == 1) {
	if (client->Send() == 1) {
	  // successful all the way
	  return 1;
	}
      }
    }
  }
  // else error
  client->PrintLastError();
  return 0;
}

int TRSendByteArray(ByteArray_Data& cad) {
  ByteArray_pkt->SetData(&cad);
  if (client->AddPacket(ByteArray_pkt) == 1) {
    if (client->Send() == 1) {
      // successfull all the way
      return 1;
    }
  }
  // else error
  client->PrintLastError();
  return 0;
}

void RegisterPackets(int type) {
  // Register packets we are going to send/accept
  // Type = 0 is normal, 1 - is just for sound client
  //   2 - is full

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
  ServerRequestAll_pkt = new ServerRequest_Packet();
  RequestUsrState_pkt = new ServerRequest_Packet();
  ByteArray_pkt = new ByteArray_Packet();
  
  // Request and response for User State
  UsrState_type = client->Register(UsrState_pkt);
  if (UsrState_type < 1) {
    // error in Register
    client->PrintLastError();
    exit(1);
  }
  fprintf(stderr, "User state packet type is %i\n", UsrState_type);

  // not saving this type, because this is something we send
  //  and thus we don't care about the type
  int result = client->Register(RequestUsrState_pkt);
  if (result < 1) {
    // error in Register
    client->PrintLastError();
    exit(1);
  }
  // set up request packet for user states
  if (RequestUsrState_pkt->AddPacketRequest(UsrState_pkt) == -1) {
    fprintf(stderr, "Failure in setting up user state request packet\n");
    exit(1);
  }

  if (type != 1) {
    // Virtual state
    result = client->Register(VirtState_pkt);
    if (result < 1) {
      client->PrintLastError();
      exit(1);
    }
    
    // Command packet
    result = client->Register(Cmd_pkt);
    if (result < 1) {
      client->PrintLastError();
      exit(1);
    }
    
    // New Position request and response
    result = client->Register(NewPos_pkt);
    if (result < 1) {
      client->PrintLastError();
      exit(1);
    }
    
    NewPosResp_type = client->Register(NewPosResp_pkt);
    if (NewPosResp_type < 1) {
      client->PrintLastError();
      exit(1);
    }
    fprintf(stderr, "New Position Response packet type is %i\n", NewPosResp_type);
    
    // Collision
    result = client->Register(Coll_pkt);
    if (result < 1) {
      client->PrintLastError();
      exit(1);
    }

    if (type == 2) {
      // add other packets
      // FSRSignals packet
      FSR_type = client->Register(Fsr_pkt);
      if (FSR_type < 1) {
	client->PrintLastError();
	exit(1);
      }
      fprintf(stderr, "FSRSignals packet type is %i\n", FSR_type);
     
      // Stairs packet
      result = client->Register(Stairs_pkt);
      if (result < 1) {
	client->PrintLastError();
	exit(1);
      }

      // Tracker packet
      Tracker_type = client->Register(Tracker_pkt);
      if (Tracker_type < 1) {
	client->PrintLastError();
	exit(1);
      }
      fprintf(stderr, "Tracker packet type is %i\n", Tracker_type);

      // Char array, let's send it so it's a write here
      result = client->Register(ByteArray_pkt);
      if (result < 1) {
	client->PrintLastError();
	exit(1);
      }
      
      // All server request packet (already registered this type of packet,
      //  but need to get the registration #/info into the packet - could
      //  have copied from RequestUsrState_pkt - this is easier)
      result = client->Register(ServerRequestAll_pkt);
      if (result < 1) {
	client->PrintLastError();
	exit(1);
      }      
      // set up this packet
      if (ServerRequestAll_pkt->AddPacketRequest(UsrState_pkt) == -1) {
	fprintf(stderr, "Failure in setting up server request all packet for user state\n");
	exit(1);
      }
      if (ServerRequestAll_pkt->AddPacketRequest(Fsr_pkt) == -1) {
	fprintf(stderr, "Failure in setting up server request all packet for user state\n");
	exit(1);
      }
      if (ServerRequestAll_pkt->AddPacketRequest(Tracker_pkt) == -1) {
	fprintf(stderr, "Failure in setting up server request all packet for user state\n");
	exit(1);
      }
    }      
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
      ch = getch();
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

double mycurs_getdouble(void) {
  int key;
  char str[8];
  int length = 0;
  double retval = 0.0;
  int decflag = 0;
  

  while ((key=getch())!=13) {

    switch(key) {
    case '-': 
      if (length==0) {
	str[length++] = key;
	addch(key);
      }      
      break;
      
    case '.':
      
      if (decflag) break;
      decflag = 1;
      
    case '0':
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9':

      if (length<5) {
	str[length++] = key;
	addch(key);
      }      
      break;

    case 8:
      
      addch(key);
      addch(' ');
      addch(key);
      length--;
      if (str[length] == '.') decflag = 0;
      str[length] = 0;
      break;

    }
  }
  retval = atof(str);
  return retval;
}

