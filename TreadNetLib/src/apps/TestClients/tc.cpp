#include "Client_header.h"

/* function declarations */
double UpdateVW(const UserState_Data &u, VirtualState_Data &v);
double mycurs_getdouble(void);

/* Variables to use from client_common */
extern char* statstr[64];
extern char* terrstr[4];
extern Client_Controller *client;
extern UserState_Packet * UsrState_pkt;


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
    fprintf(stderr," \n usage:  tc server port <hertz>\n\n");
    printf(" \n\n and the best way to run this in unix is: \n");
    printf("  (tc host portnum > /dev/tty) >& /dev/console \n");
    printf("  This way, all stderr messages go to your console window.\n");
    printf("\n  In Windows, simply run from a console:\n");
    printf("  tc host portnum <hertz> \n");
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
  RegisterPackets();
  printf("Registered all Packets\n");

  // Set some initial values - use data classes
  UserState_Data us_data;
  VirtualState_Data vs_data;
  
  // User state first, Only can change the position and facing
  //  as all other values are determined by server
  us_data.SetPosEasting(0.0);
  us_data.SetPosNorthing(0.0);

  // Set facing North
  us_data.SetFacingEasting(0.0);
  us_data.SetFacingNorthing(1.0);
  //  TRSetUserElevation(&user_state,1);
  //  TRSetEyeElevation(&user_state,1.8);
  
  // Virtual state
  TRVector surf_normal(0,0,1);
  TRTerrain terrain = NORMAL_TERRAIN;
  TRVector virtual_force(0);
  float movement_scale = 1;
  float user_scale = 1;

  vs_data.SetSurfaceNormal(surf_normal);
  vs_data.SetTerrain(terrain);
  vs_data.SetVirtualForce(virtual_force);
  vs_data.SetMovementScale(movement_scale);
  vs_data.SetUserScale(user_scale);

  // Treadport must be disabled for New Position to work
  TRDisable();
  sleep(1);
     
  // Need to check that status is disabled
  TRGetUserState();
  while (UsrState_pkt->GetStatus() & ENABLED_BIT) {
    printf(" awaiting disabled mode on the treadport...\n");
    sleep(2);
    TRGetUserState();
  }

  /* start the user out on the top of the hill */
  TRSetPosition(us_data, vs_data);
 
  double elevation;
  elevation = UpdateVW(us_data,vs_data);

  // Initialize window
  InitWindow();
  
  int quit=0;
  int counter = 0;

  /* ok, lets do something interesting */
  TRVector tmpfacing;
  UserState_Data tmp_usd = us_data;
  VirtualState_Data tmp_vsd = vs_data;


  while(!quit) {
    
    /* request the treadport position and the status */
    TRGetUserState();

    // check if valid position returned
    if (!UsrState_pkt->PosValid()) {
      // invalid position returned, which means we lost the server and
      //  when it rebooted it no longer has a valid position, set it
      TRDisable();
      TRSetPosition(us_data, vs_data);
    }
    else {
      // new data received, get copy of it
      UsrState_pkt->GetCopyOfData(&us_data);
    }

    /* check to see if the user wants to
       command start/stop or change the terrain. */
    switch (GetCharacter()) {
	
    case 'D':
      /* Disable */
      TRDisable();
      break;
      
    case 'E':
      /*  */
      TREnable();
      break;
      
    case 'M': 
      movement_scale *= 1.1;
      //      fprintf(stderr," increased to %f \n", movement_scale);
      vs_data.SetMovementScale(movement_scale);
      break;

    case 'm':
      movement_scale *= .9;
      vs_data.SetMovementScale(movement_scale);
      break;
      
    case 'U':
      user_scale *= 1.1;
      vs_data.SetUserScale(user_scale);
      break;

    case 'u':
      user_scale *= .9;
      vs_data.SetUserScale(user_scale);
      break;

    case 'V':
      // increase virtual force in direction of facing by 10 N
      virtual_force.X() += 10*us_data.GetFacingEasting();
      virtual_force.Y() += 10*us_data.GetFacingNorthing();
      virtual_force.Z() += 1;
      vs_data.SetVirtualForce(virtual_force);
      break;

    case 'v':
      // decrease virtual force in direction of facing by 10 N
      virtual_force.X() -= 10*us_data.GetFacingEasting();
      virtual_force.Y() -= 10*us_data.GetFacingNorthing();
      virtual_force.Z() -= 1;
      vs_data.SetVirtualForce(virtual_force);
      break;

    case 't':
      /* change the terrain */
      terrain++;
      if (terrain>ICY_TERRAIN) 
	terrain = NORMAL_TERRAIN;
      vs_data.SetTerrain(terrain);
      break;
      
    case 'j':
      /* jump to new location */
      if (us_data.GetStatus() & ENABLED_BIT ) 
	break;

      double easting;
      double northing;
      double tmp_elev;
      double invnorm;

      printwvis1(" \n INPUT: \n");
      printwvis1("   new location: \n");
      printwvis1("    Easting:  ");
      easting=mycurs_getdouble();
      printwvis1("\n    Northing: ");
      northing=mycurs_getdouble();
      printwvis1("\n    Facing X: ");
      tmpfacing.X(mycurs_getdouble());
      printwvis1("\n    Facing Y: ");
      tmpfacing.Y(mycurs_getdouble());
      tmpfacing.Z(0);

      invnorm = tmpfacing.Norm();

      if (invnorm < EPS) {
	tmpfacing.X(1.0);
	tmpfacing.Y(0.0);
      }
      else {
	tmpfacing /= invnorm;
      }
      tmp_usd.SetPosEasting(easting);
      tmp_usd.SetPosNorthing(northing);
      tmp_usd.SetFacing(tmpfacing);
      tmp_vsd = vs_data;
      tmp_elev = UpdateVW(tmp_usd,tmp_vsd);
      
      if (TRSetPosition(tmp_usd,tmp_vsd)) {
	/* we were successful, so we should
	   update ourselves the same way. */
	us_data = tmp_usd;
	vs_data = tmp_vsd;
      }
      // otherwise an error occurred, so don't save new user state
      //  or virtual state data
      
      break;
      
    case 'q':
      client->Disconnect();
      sleep(1);
      quit = 1;
      break;

    }
    
    // Get new elevation, virtual state data
    elevation = UpdateVW(us_data,vs_data);
    // Send virtual state to server
    TRSetVirtualState(vs_data);

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
    printw1("    E:  Enable Treadport      \n");
    printw1("    D:  Disable  Treadport      \n");
    printw1("    t:  change the terrain      \n");
    printw1("    j:  jump user to new position (when stopped)      \n");
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
	    elevation+us_data.GetUserElevation());

    printw4(" Eye Position:       %7.2lf E  %7.2lf N %7.2lf E        \n",
	    us_data.GetPosEasting() + us_data.GetEyeEasting(),
	    us_data.GetPosNorthing() + us_data.GetEyeNorthing(),
	    elevation+us_data.GetEyeElevation());
	   
    printw4(" Treadport Facing:    %7.2lf X  %7.2lf Y  %7.2lf Z        \n",
	   facing.X(),facing.Y(),facing.Z());
    printw4(" Treadport Velocity:  %7.2lf X  %7.2lf Y  %7.2lf Z        \n",
	   velocity.X(),velocity.Y(),velocity.Z());
    printw1(" -------------------------------------        \n ");
    printw2(" Current Terrain Setting: %s          \n",
	   terrstr[terrain]);

#ifdef RTI_VXWORKS
    }
#endif

    /* sleep a small while */
    
    //#if   defined(sparc_solaris)
    //      struct timespec r;
    //      r.tv_sec = 0;
    //      r.tv_nsec = 2000000000; /* ??? 30 hz? */
    //      usleep(r.tv_nsec/1000);
    //#else
    my_sync.block();
    //#endif
    
    counter++;
  }/* end while */

  StopWindow();
  
  return 1;  
}

/******************************* Virtual World Info ***********************************/

#define MAXX +10.0
#define MINX -10.0
#define MAXY +10.0
#define MINY -10.0

/* this defines sqrt(2)/2 */
#define R2O2 0.70710678

double UpdateVW(const UserState_Data &u, VirtualState_Data &v) {
  
  /* here we assume the world is a gausian
     hill defined on x = [MINX,MAXX], 
     y = [MINY,MAXY], with z = 5*exp(-(x^2+y^2)/100).
     
     The square world is surrounded by walls that
     produce collisions. 
     
     There is no virtual force yet. */

  double x,y,z;
  TRVector normal,collisionn,collisionp;

  x = u.GetPosEasting();
  y = u.GetPosNorthing();

  double elevation;
  elevation = z = 5.0*exp(-(x*x+y*y)/100);

  normal.X() = 2*z*x/100;
  normal.Y() = 2*z*y/100;
  normal.Z() = 1;
  normal /= normal.Norm();

  v.SetSurfaceNormal(normal);
  
  /* worry about collisions with the walls
     This is really ugly code, but it tries
     to do the tom-style basic blending of normals
     at the contact between two surfaces. */
  

  /*  if (x>MAXX) {
    if (y>MAXY) {  
      collisionn.X = -R2O2;
      collisionn.Y = -R2O2;
      collisionn.Z = 0;
      collisionp.X = MAXX;
      collisionp.Y = MAXY;
    }
    else if ( y<MINY) {
      collisionn.X = -R2O2;
      collisionn.Y = R2O2;
      collisionn.Z = 0;
      collisionp.X = MAXX;
      collisionp.Y = MINY;
    }
    else {
      collisionn.X = -1;
      collisionn.Y = 0;
      collisionn.Z = 0;
      collisionp.X = MAXX;
      collisionp.Y = y;

    }
    //    TRSetCollision(v,collisionn,collisionp);    
  }
  else if (x<MINX) {
    if (y>MAXY) {  
      collisionn.X = +R2O2;
      collisionn.Y = -R2O2;
      collisionn.Z = 0;
      collisionp.X = MINX;
      collisionp.Y = MAXY;
    }
    else if ( y<MINY) {
      collisionn.X = +R2O2;
      collisionn.Y = +R2O2;
      collisionn.Z = 0;
      collisionp.X = MINX;
      collisionp.Y = MINY;
    }
    else {
      collisionn.X = 1;
      collisionn.Y = 0;
      collisionn.Z = 0;
      collisionp.X = MINX;
      collisionp.Y = y;
     }
    //    TRSetCollision(v,collisionn,collisionp);
  }
  else {
    if (y>MAXY) {  
      collisionn.X = 0;
      collisionn.Y = -1;
      collisionn.Z = 0;
      collisionp.X = x;
      collisionp.Y = MAXY;
      //      TRSetCollision(v,collisionn,collisionp);
    }
    else if (y<MINY) {
      collisionn.X = 0;
      collisionn.Y = 1;   
      collisionn.Z = 0;
      collisionp.X = x;
      collisionp.Y = MINY;
      //      TRSetCollision(v,collisionn,collisionp);
    }
    else {
      //      TRClearCollision(v);
    }
  }
  */
  return elevation;
}  


