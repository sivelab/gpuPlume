#include "Client_header.h"

/* variables to use from client_common */
extern Client_Controller *client;
extern UserState_Packet * UsrState_pkt;
extern char* statstr[64];

#ifndef RTI_VXWORKS
int main(int argc, char *argv[]) {
#else
int vxmain(int arg1, int arg2, int arg3) {
#endif
    
  char hostname[128];
  short port;
  double hz_request = 64;

  // process command line
#if defined (RTI_VXWORKS)
  // vxworks allows only 9 int args, so take first arg to define which machine, and
  //  the second arg as the port number - I know it's messy
  if (arg1 > 3 || arg2 == 0) {
    printf(" \n usage: sp <vxmain address>, server num, port, <hertz>\n\n");
    printf(" Current servers available are \n");
    printf(" 0 - roboworks\n");
    printf(" 1 - roboshell\n");
    printf(" 2 - buzzworm\n");
    printf(" 3 - hovenweep\n");
    exit(1);
  }

  switch (arg1) {
  case 0:
    strcpy(hostname, "roboworks");
    break;
  case 1:
    strcpy(hostname, "roboshell");
    break;
  case 2:
    strcpy(hostname, "buzzworm");
    break;
  case 3:
    strcpy(hostname, "hovenweep");
    break;
  }

  port = arg2;
  if (arg3 != 0) {
    hz_request = (double)pow(2.0, (int)(log(arg3)/log(2) + 0.5));
  }
#else
  if (argc < 3) {
    fprintf(stderr," \n usage:  sc server port <hertz>\n");
    printf(" \n and the best way to run this in unix is: \n");
    printf("  (sc host portnum <hertz> > /dev/tty) >& /dev/console \n");
    printf("  This way, all stderr messages go to your console window.\n");
    printf("\n  In Windows, simply run from a console:\n");
    printf("  sc host portnum <hertz>\n");
    printf(" Note: hertz should be a power of 2 from 1 - 8192\n"); 
    exit(1);
  }
  strcpy(hostname,argv[1]);
  port = atoi(argv[2]);
  if (argc > 3) {
    hz_request = atof(argv[3]);
    // convert hertz to power of 2
    hz_request = (double)pow(2.0, (int)(log(hz_request)/log(2) + 0.5));
  }
#endif
  
  rtcSync my_sync((int)hz_request);

  // Create client controller and connect to server
  client = new Client_Controller();
  //  client->SetVerbosity(1);

  switch (client->ConnectToServer(hostname,port)) {
  case SOCKET_ERR:
    // socket problem
    client->PrintLastError();
    exit(1);
    break;
  case NO_SERVER:
    // server does not have a socket open
    client->PrintLastError();
    exit(1);
    break;
  default:
    // connection successful
    fprintf(stderr, "Connection successful to %s:%d\n", hostname, port);
    break;
  }

  // Register all packets (it creates the packets now)
  RegisterPackets(1);
  printf("Registered all Packets\n");

  // Set some initial values - use data classes
  UserState_Data us_data;
  
  // User state first, Only can change the position and facing
  //  as all other values are determined by server
  us_data.SetPosEasting(0.0);
  us_data.SetPosNorthing(0.0);

  // Set facing North
  us_data.SetFacingEasting(0.0);
  us_data.SetFacingNorthing(1.0);
  
  TRGetUserState();

  // Initialize window
  InitWindow();
    
  int quit=0;
  int counter = 0;
  
  /* ok, lets do something interesting */
  UserState_Data tmp_usd = us_data;

  while(!quit) {
    
    /* request the treadport position and the status */
    TRGetUserState();

    // new data received, get copy of it
    UsrState_pkt->GetCopyOfData(&us_data);

    /* check to see if the user wants to
       command start/stop or change the terrain. */
    switch (GetCharacter()) {
	
    case 'q':
      client->Disconnect();
      sleep(1);
      quit = 1;
      break;

    }
    
    TRVector facing   = us_data.GetFacing();
    TRVector velocity = us_data.GetVelocity();

    /* display this stuff on the text screen*/
    EraseWindow();
#ifdef RTI_VXWORKS
    // only printing at 1 hz since lack curses, or windows
    if (counter%100 == 0) {
#endif
    printw2(" Counter: %d\n",counter);
    printw3(" Treadport Client connected to %s, port %d        \n\n",
	   hostname, port);
    printw1(" -------------------------------------      \n");
    printw1(" Keys: \n");
    printw1("    q:  quit      \n");
    printw1(" -------------------------------------      \n");
    printw2(" Treadport Status: %s\n",statstr[us_data.GetStatus()]);
    printw1(" -------------------------------------      \n");
    printw3(" Treadport Position:  %7.2lf East %7.2lf North     \n",
	   us_data.GetPosEasting(),
	   us_data.GetPosNorthing());
    printw4(" User Position:       %7.2lf E  %7.2lf N %7.2lf E        \n",
	    us_data.GetPosEasting() + us_data.GetUserEasting(),
	    us_data.GetPosNorthing() + us_data.GetUserNorthing(),
	    us_data.GetUserElevation());

    printw4(" Eye Position:       %7.2lf E  %7.2lf N %7.2lf E        \n",
	    us_data.GetPosEasting() + us_data.GetEyeEasting(),
	    us_data.GetPosNorthing() + us_data.GetEyeNorthing(),
	    us_data.GetEyeElevation());
	   
    printw4(" Treadport Facing:    %7.2lf X  %7.2lf Y  %7.2lf Z        \n",
	   facing.X(),facing.Y(),facing.Z());
    printw4(" Treadport Velocity:  %7.2lf X  %7.2lf Y  %7.2lf Z        \n",
	   velocity.X(),velocity.Y(),velocity.Z());

#ifdef RTI_VXWORKS
    }
#endif
    
    /* sleep a small while */

    //#if   defined(sparc_solaris)
    //      struct timespec r;
    //      r.tv_sec = 0;
    //      r.tv_nsec = 2000000000; /* ??? 30 hz? */
    //    usleep(r.tv_nsec/1000);
    //#else
    my_sync.block();
    //#endif
    
    counter++;
  }/* end while */

  StopWindow();
  
  return 1;  
}


