#include "Client_header.h"

/* function declarations */
double UpdateVW(const UserState_Data &u, VirtualState_Data &v, int slope_ang, Collision_Data &c, Stairs_Data &s, TRTerrain surf_terrain);
double mycurs_getdouble(void);

/* Variables to use from client_common */
extern char* statstr[64];
extern char* terrstr[4];
extern Client_Controller *client;
// only need packets which I'm receiving, the ones I'm sending I'll just
//  create the data portions and pass them to the function
extern UserState_Packet * UsrState_pkt;
extern Tracker_Packet * Tracker_pkt;
extern FSRSignals_Packet * Fsr_pkt;

#ifndef RTI_VXWORKS
int main(int argc, char *argv[]) {
#else
int vxmain(int arg1, int arg2, int arg3, int arg4) {
#endif
  char hostname[128];
  short port;
  int slope_angle = 0;
  double hz_request = 64;

  // process command line
#if defined (RTI_VXWORKS)
  // vxworks allows only 9 int args, so take first arg to define which machine, and
  //  the second arg as the port number - I know it's messy
  if (arg1 > 3 || arg2 == 0) {
    printf(" \n usage: sp <vxmain address>, server num, port, <slope angle>, <hertz>\n\n");
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
  slope_angle = arg3;
  if (arg4 != 0) {
    hz_request = (double)pow(2.0, (int)(log(arg4)/log(2) + 0.5));
  }
#else  
  if (argc < 3) {
    fprintf(stderr," \n usage:  tc_single_slope server port <slope angle> <hertz>\n\n");
    printf(" \n\n and the best way to run this in unix is: \n");
    printf("  (tc_single_slope host portnum <angle> <hertz> > /dev/tty) >& /dev/console \n");
    printf("  This way, all stderr messages go to your console window.\n");
    printf("\n  In Windows, simply run from a console:\n");
    printf("  tc_single_slope host portnum <angle> <hertz>\n");
    exit(1);
  }
  strcpy(hostname,argv[1]);
  port = atoi(argv[2]);
  if (argc > 3) {
    slope_angle = atoi(argv[3]);
    if (argc > 4) {
      hz_request = atof(argv[4]);
      // convert hertz to power of 2
      hz_request = (double)pow(2.0, (int)(log(hz_request)/log(2) + 0.5));
    }
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
  RegisterPackets(2);
  printf("Registered all Packets\n");

  // Set some initial values - use data classes
  UserState_Data us_data;
  VirtualState_Data vs_data;
  Collision_Data cd_data;
  Stairs_Data sd_data;
  ByteArray_Data ba_data;
  
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
  TRTerrain surface_terrain = NORMAL_TERRAIN;
  TRVector virtual_force(0);
  float movement_scale = 1;
  float user_scale = 1;

  vs_data.SetSurfaceNormal(surf_normal);
  vs_data.SetTerrain(surface_terrain);
  vs_data.SetVirtualForce(virtual_force);
  vs_data.SetMovementScale(movement_scale);
  vs_data.SetUserScale(user_scale);

  // Collision Data (just setting zero, though already is so can see
  //  the function calls)
  cd_data.SetNoCollision();
  cd_data.SetCollisionNormal(TRVector(0));
  cd_data.SetCollisionPos(TRPoint(0));

  // Stairs Data - set no commands or not on stairs
  sd_data.SetNoStairsCmd();
  sd_data.SetIgnoreStairsDir();

  // Character array data
  char * message = strdup("Hello server");
  int i;
  for (i = 0; i < strlen(message); i++) {
    ba_data.SetByte(i, message[i]);
  }
  ba_data.SetByte(i, 0);
  ba_data.SetByte(i+1,'5');
  ba_data.SetLastElementInd(i+1);

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
  elevation = UpdateVW(us_data,vs_data,slope_angle, cd_data, sd_data, surface_terrain);

  // Initialize window
  InitWindow();
      
  int quit=0;
  int counter = 0;

  /* ok, lets do something interesting */
  TRVector tmpfacing;
  UserState_Data tmp_usd = us_data;
  VirtualState_Data tmp_vsd = vs_data;
  Collision_Data tmp_cd = cd_data;
  Stairs_Data tmp_sd = sd_data;

  while(!quit) {
    
    /* request the treadport position, status, fsrs, and tracker */
    TRGetAllServerInfo();

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

    // reset stairs command since in current implementation I send a stair packet
    //  every cycle even if not necessary (could just send when there is a new command,
    //  or different state of the stairs - this is how I would implement really)
    sd_data.SetNoStairsCmd();

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
      
    case 'c':
      sd_data.SetDecConstFactCmd();
      break;

    case 'C':
      sd_data.SetIncConstFactCmd();
      break;

    case 'S':
      sd_data.SetSwitchStairsCmd();
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
      /* change the surface terrain */
      surface_terrain++;
      // don't allow stairs terrain here, that is set elsewhere
      if (surface_terrain == STAIRS_TERRAIN)
	surface_terrain++;

      if (surface_terrain>ICY_TERRAIN) 
	surface_terrain = NORMAL_TERRAIN;

      if (vs_data.GetTerrain() != STAIRS_TERRAIN) {
	// don't override stairs terrain
	vs_data.SetTerrain(surface_terrain);
      }
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
      tmp_elev = UpdateVW(tmp_usd,tmp_vsd,slope_angle,tmp_cd,tmp_sd,surface_terrain);
      
      if (TRSetPosition(tmp_usd,tmp_vsd)) {
	/* we were successful, so we should
	   update ourselves the same way. */
	us_data = tmp_usd;
	vs_data = tmp_vsd;
	cd_data = tmp_cd;
	sd_data = tmp_sd;
      }
      // otherwise an error occurred, so don't save new user state
      //  or virtual state data
      
      break;
      
    case 'q':
      client->Disconnect();
      sleep(3);
      quit = 1;
      break;

    }
    
    // Get new elevation, virtual state data
    elevation = UpdateVW(us_data,vs_data,slope_angle,cd_data,sd_data,surface_terrain);
    // Send virtual state to server
    TRSetAllVirtInfo(vs_data,cd_data,sd_data);
    // send character array separately for now
    TRSendByteArray(ba_data);
    
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
	   terrstr[vs_data.GetTerrain()]);

    printw1(" FSRSignals : Right Heel     Toe    Left Heel     Toe\n");
    printw5("              %7.2f     %7.2f   %7.2f    %7.2f \n",
	    Fsr_pkt->GetValueRightHeel(), Fsr_pkt->GetValueRightToe(),
	    Fsr_pkt->GetValueLeftHeel(), Fsr_pkt->GetValueLeftToe());
    printw4(" Tracker Pos :   %7.2f X  %7.2f Y  %7.2f Z    \n",
	    Tracker_pkt->GetPosX(), Tracker_pkt->GetPosY(), Tracker_pkt->GetPosZ());
    TRVector ov = Tracker_pkt->GetOrientVector();
    printw5(" Tracker Orient : %7.2f qX  %7.2f qY  %7.2f qZ  %7.2f qW   \n",
	    ov.X(), ov.Y(), ov.Z(), Tracker_pkt->GetOrientRot());


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

#define MAXX +50.0
#define MINX -100.0
#define MAXY +50.0
#define MINY -50.0

/* this defines sqrt(2)/2 */
#define R2O2 0.70710678

double UpdateVW(const UserState_Data &u, VirtualState_Data &v, int slope_ang, Collision_Data &c, Stairs_Data &s, TRTerrain surf_terrain) {
  
  // here we assume the world is a single slope along X axis defined by slope_angle
  // It starts at axis 0,0 (actually up slope is negative X or West)

  double x,y,z;
  TRVector normal,collisionn,collisionp;

  x = u.GetPosEasting();
  y = u.GetPosNorthing();

  double elevation;
  elevation = z = -x*sin(slope_ang/180.0*3.1428);

  normal.X() = sin(slope_ang/180.0*3.1428);
  normal.Y() = 0;
  normal.Z() = cos(slope_ang/180.0*3.1428);

  v.SetSurfaceNormal(normal);
  
  /* worry about collisions with the walls
     This is really ugly code, but it tries
     to do the tom-style basic blending of normals
     at the contact between two surfaces. */
  
  // Put a wall around the whole hill
  collisionn.Z() = 0; // always

  // for code simplicity, assume Collision
  c.SetCollision();

  if (x>MAXX) {
    // banged up against east wall
    collisionp.X() = MAXX;
    collisionp.Z() = -MAXX*sin(slope_ang/180.0*3.1428); 
    if (y>MAXY) {  
      collisionn.X() = -R2O2;
      collisionn.Y() = -R2O2;
      collisionp.Y() = MAXY;
    }
    else if ( y<MINY) {
      collisionn.X() = -R2O2;
      collisionn.Y() = R2O2;
      collisionp.Y() = MINY;
    }
    else {
      collisionn.X() = -1;
      collisionn.Y() = 0;
      collisionp.Y() = y;
    }
  }
  else if (x<MINX) {
    // banged into west wall
    collisionp.X() = MINX;
    collisionp.Z() = -MINX*sin(slope_ang/180.0*3.1428); 
    if (y>MAXY) {  
      collisionn.X() = +R2O2;
      collisionn.Y() = -R2O2;
      collisionp.Y() = MAXY;
    } 
    else if ( y<MINY) {
      collisionn.X() = +R2O2;
      collisionn.Y() = +R2O2;
      collisionp.Y() = MINY;
    }
    else {
      collisionn.X() = 1;
      collisionn.Y() = 0;
      collisionp.Y() = y;
    }
  } 
  else {
    collisionp.Z() = z;
    if (y>MAXY) {  
      // banged into north wall
      collisionn.X() = 0;
      collisionn.Y() = -1;
      collisionp.X() = x;
      collisionp.Y() = MAXY;
    }
    else if (y<MINY) {
      // banged into south wall
      collisionn.X() = 0;
      collisionn.Y() = 1;   
      collisionp.X() = x;
      collisionp.Y() = MINY;
    }
    else {
      c.SetNoCollision();
    }
  }
  if (c.GetCollisionFlag() == COLLISION) {
    // set collision normal and position
    c.SetCollisionNormal(collisionn);
    c.SetCollisionPos(collisionp);
  }
  
  // Now Stairs - putting them between x [-50 -100], letting their
  //  equivalent slope (which currently is ignored by server since it 
  //  implements just one slope of stairs) be the slope entered.
  if (x < -50) {
    // check which direction facing (upstairs is negative x, since that's upslope)
    if (u.GetFacingEasting() > 0) {
      // down
      s.SetDownStairsDir();
    }
    else {
      // up
      s.SetUpStairsDir();
    }
    v.SetTerrain(STAIRS_TERRAIN);
  }
  else {
    // no stairs
    s.SetIgnoreStairsDir();
    v.SetTerrain(surf_terrain);
  }

  return elevation;
}  
