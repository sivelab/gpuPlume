#include "Server_header.h"

/* Variables to use from client_common */
extern char* statstr[64];
extern char* terrstr[4];
extern char* speedkey_string;
extern char* turnkey_string;
extern Server_Controller *server;
extern ServerRequest_Packet * ServerRequest_pkt;
extern int ServerRequest_type;
extern UserState_Packet * UsrState_pkt;
extern VirtualState_Packet * VirtState_pkt;
extern int VirtState_type;
extern Command_Packet * Cmd_pkt;
extern int Cmd_type;
extern NewPos_Packet * NewPos_pkt;
extern NewPosResponse_Packet * NewPosResp_pkt;
extern int NewPos_type;
extern Collision_Packet * Coll_pkt;
extern int Coll_type;

#ifndef RTI_VXWORKS
int main(int argc, char *argv[]) {
#else
int vxmain(int arg1, int arg2) {
#endif
  char              hostname[128];
  short             port;
  double hz_request = 64;

  /* process command line */  
#if defined (RTI_VXWORKS)
  // vxworks allows only 9 int args, so take first arg as port number - I know it's messy
  if (arg1 <= 5000) {
    printf(" \n usage: sp <vxmain address>, port, <hertz>\n\n");
    printf("            portnum must be > 5000\n\n");
    exit(1);
  }

  port = arg1;
  if (arg2 != 0) {
    hz_request = (double)pow(2.0, (int)(log(arg2)/log(2) + 0.5));
  }
#else  
  if (argc < 2
      || (port =(atoi(argv[1]))) <= 5000 ) {

    printf(" \n usage:  ts portnum <hertz>\n\n");
    printf("            portnum must be > 5000\n\n");
    printf(" \n\n and the best way to run this is: \n");
    printf("  (ts portnum <hertz> > /dev/tty) >& /dev/console \n");
    printf("  This way, all stderr messages go to your console window.\n");
    printf("\n\n  In windows, simply run from a console window: \n");
    printf("  ts portnum <hertz>\n");
    exit(1);
  }

  if (argc > 2) {
    // also sent the requested hertz
    hz_request = atof(argv[2]);
    hz_request = (double)pow(2.0, (int)(log(hz_request)/log(2) + 0.5));
  }
#endif
  
  rtcSync my_sync((int)hz_request);

  double timestep = 1 / hz_request;
  int countdown_val = ((int)(5/timestep));
  
  // create server controller
  server = new Server_Controller();

  // Start Server
  //  server->SetVerbosity(1);
  if (server->StartServer(port) != 1) {
    fprintf(stderr," Server didn't start...\n");
    server->PrintLastError();
    exit(1);
  }

  // get hostname
  gethostname(hostname,sizeof(hostname));

  // create all packets
 
  // Register all packets (it creates the packets now)
  RegisterPackets();

  // Set some initial values - use data classes to hold values
  UserState_Data us_data;
  VirtualState_Data vs_data;

  TRVector facing(0,1,0); // face north
  TRVector velocity(0);
  TRVector surface_normal(0,0,1); // flat ground
  double easting = 0;
  double northing = 0;
  TRStatus mystatus = INACTIVE_DISABLED;

  us_data.SetPosEasting(easting);
  us_data.SetPosNorthing(northing);
  
  // not handling any user/eye offsets, except elevation
  us_data.SetUserEasting(0);
  us_data.SetUserNorthing(0);
  us_data.SetUserElevation(1.0);

  us_data.SetEyeEasting(0);
  us_data.SetEyeNorthing(0);
  us_data.SetEyeElevation(1.8);

  us_data.SetFacing(facing);
  us_data.SetVelocity(velocity);
  us_data.SetStatus(mystatus);

  // virtual values initalized until receive first values
  vs_data.SetTerrain(NORMAL_TERRAIN);
  vs_data.SetSurfaceNormal(surface_normal);
  vs_data.SetVirtualForce(TRVector(0));

  // collision packet starts out with no collision

  // Speed value
  double speed = 0.0;

  // counter, and determining if done
  int quit    = 0;
  int counter = 0;

  // variables for handling ramp up and down
  int status_count;
  double status_rampval;

  // dummy for only printing virt state once received
  int received_virtstate    = 0;

  // holds whether collided
  int collided = 0;

  // variables for receiving and handling message
  int sc_result;
  int packet_type;
   
  // Initialize window
  InitWindow();
 
  EraseWindow();

  /* ok, lets do something interesting */
  double newx, newy;
  int ch;
  while(!quit) {
    ch = GetCharacter();
    switch (ch) {
	
    case KEY_UP:
      if (mystatus == ACTIVE_ENABLED) {
	speed += .1;
	if (speed>MAX_SPEED) speed = MAX_SPEED;
      }
      break;
      
    case KEY_DOWN:
      if (mystatus == ACTIVE_ENABLED) {
	speed -= .1;
	if (speed < -MAX_SPEED ) speed = -MAX_SPEED;
      }
      break;
      
    case KEY_LEFT:
      newx = COS5*facing.X() - SIN5*facing.Y();
      newy = +SIN5*facing.X() + COS5*facing.Y();
      facing.X() = newx;
      facing.Y() = newy;
      break;
      
    case KEY_RIGHT:
      newx =  COS5*facing.X() + SIN5*facing.Y();
      newy = -SIN5*facing.X() + COS5*facing.Y();
      facing.X() = newx;
      facing.Y() = newy;
      break;
      
    case 'k': /* KILL */
      if (mystatus & (DISABLED_BIT | RAMPING_DOWN_BIT)) {
	mystatus = INACTIVE_DISABLED;
      }
      else {
	mystatus = INACTIVE_ENABLED;
      }
      status_count = 0;
      speed = 0;
      break;      
      
    case 's': /* START */
      if (mystatus == INACTIVE_ENABLED) {
	status_count = countdown_val;
	mystatus = ACTIVE_ENABLED_RAMPING_UP;
      }
      else if (mystatus == INACTIVE_DISABLED) {
	mystatus = ACTIVE_DISABLED;
      }
      break;
      
    case 'q':
      quit = 1;
    }

    // update changes to us_data
    us_data.SetStatus(mystatus);
    us_data.SetFacing(facing);
  
    /* handle the update in the users virtual position */
    if ((mystatus & ACTIVE_BIT) && (mystatus & ENABLED_BIT)) {
      // if active and enabled (can be ramping as well)
      surface_normal = vs_data.GetSurfaceNormal();
      
      TRVector facing3d = convert2dheadingto3d(surface_normal,facing);

      /* status_rampval is just used to speed
	 us up or slow us down when we
	 are changing modes; */

      velocity = facing3d * speed * status_rampval;
      
      easting  += velocity.X() * timestep;
      northing += velocity.Y() * timestep;
      /* change in elevation depends on the terrain at new spot */

      us_data.SetPosEasting(easting);
      us_data.SetPosNorthing(northing);

      // going to change user/eye easting/northing for fun
      us_data.SetUserEasting(velocity.X() * .1);
      us_data.SetUserNorthing(velocity.Y() * .1);
	
      us_data.SetEyeEasting(velocity.X() * .1);
      us_data.SetEyeNorthing(velocity.Y() * .1);
	
      us_data.SetVelocity(velocity);
    }
    
    // receive message
    sc_result = server->Receive();
    if (sc_result == 1) {
      // received some message, parse it
      packet_type = server->GetPacketType();
      while (packet_type != -1) {
	if (packet_type == ServerRequest_type) {
	  server->GetPacketFromMsg(ServerRequest_pkt);

	  // Parse server request (just query, not actually parsing full bit field
	  //  to see if extra packets requested that I can't return)
	  if (ServerRequest_pkt->IsPacketRequested(UsrState_pkt)) {
	    // requested a user state packet - send it
	    UsrState_pkt->SetData(&us_data);
	    if (server->AddPacket(UsrState_pkt) != 1) {
	      // some error in add packet, print it
	      server->PrintLastError();
	    }
	  }
	}
	else if (packet_type == VirtState_type) {
	  // received a virtual state packet
	  server->GetPacketFromMsg(VirtState_pkt);

	  // update virtual info
	  surface_normal = VirtState_pkt->GetSurfaceNormal();
	  VirtState_pkt->GetCopyOfData(&vs_data);
	  received_virtstate = 1;
	}
	else if (packet_type == Cmd_type) {
	  server->GetPacketFromMsg(Cmd_pkt);
	  switch(Cmd_pkt->GetCmd()) {
	  case ENABLE_CMD:
	    // graphics is enabled
	    if (mystatus == INACTIVE_DISABLED) {
	      mystatus = INACTIVE_ENABLED;
	    }
	    else if (mystatus == ACTIVE_ENABLED_RAMPING_DOWN) {
	      status_count = countdown_val - status_count;
	      mystatus = ACTIVE_ENABLED_RAMPING_UP;
	    }
	    else if (mystatus == ACTIVE_DISABLED && received_virtstate) {
	      status_count = countdown_val;
	      mystatus = ACTIVE_ENABLED_RAMPING_UP;
	    }
	    // otherwise has no affect
	    
	    break;
	  case DISABLE_CMD:
	    if (mystatus == ACTIVE_ENABLED) {
	      status_count = countdown_val;
	      mystatus = ACTIVE_ENABLED_RAMPING_DOWN;
	    }
	    else if (mystatus == ACTIVE_ENABLED_RAMPING_UP) {
	      status_count = countdown_val - status_count;
	      mystatus = ACTIVE_ENABLED_RAMPING_DOWN;
	    }
	    else if (mystatus == INACTIVE_ENABLED) {
	      mystatus = INACTIVE_DISABLED;
	    }
	    break;
	  default:
	    // bad command sent
	    fprintf(stderr, "Received unknown command %i\n", Cmd_pkt->GetCmd());
	    break;
	  }
	  //set new status
	  us_data.SetStatus(mystatus);
	}
	else if (packet_type == NewPos_type) {
	  // received new position
	  server->GetPacketFromMsg(NewPos_pkt);

	  if (mystatus & DISABLED_BIT) {
	    // only can set new position if disabled
	    
	    // get both user and virtual data, though only update
            //  position and facing for user data, other is up to server
            UserState_Data tmp_usd = NewPos_pkt->GetUserState();
	    us_data.SetPos(tmp_usd.GetPos());
            us_data.SetFacing(tmp_usd.GetFacing());
	    vs_data = NewPos_pkt->GetVirtualState();
	    facing = us_data.GetFacing();
	    easting = us_data.GetPosEasting();
	    northing = us_data.GetPosNorthing();
	    surface_normal = vs_data.GetSurfaceNormal();
	    
	    // set valid state for user state packets now
	    us_data.SetValidFlag();

	    // set received some virtual state
	    received_virtstate = 1;
	    
	    // send back a valid response
	    NewPosResp_pkt->SetSuccessFlag();
	    if (server->AddPacket(NewPosResp_pkt) != 1) {
	      // some error in add packet, print it
	      server->PrintLastError();
	    }
	  }
	  else {
	    // not disabled, send back a invalid response
	    NewPosResp_pkt->SetFailFlag();
	    if (server->AddPacket(NewPosResp_pkt) != 1) {
	      // some error in add packet, print it
	      server->PrintLastError();
	    }
	  }
	}
	else if (packet_type == Coll_type) {
	  server->GetPacketFromMsg(Coll_pkt);
	}
	else {
	  // received unknown packet type, disconnect client
	  server->DisconnectClient();
	}
	packet_type = server->GetPacketType();
      }
    }
    else {
      // no message received, or could be error
      switch (sc_result) {
      case REGISTRATION_ERROR:
      case MSG_BUF_BUSY:
      case MSG_ERROR:
	server->PrintLastError();
	break;
      }
      // otherwise no_msg or no_client, no need to print an error
    }

    // send any packets that might be in the queue
    sc_result = server->Send();
    switch (sc_result) {
    case MSG_ERROR:
    case MSG_BUF_BUSY:
    case MSG_BUF_OVERRUN:
    case PACKET_ERROR:
      server->PrintLastError();
      break;
    }

    /* handle state transitions */ 
    
    if (status_count > 0) {
      status_count--;
      if (mystatus & RAMPING_UP_BIT) {
	if (status_count == 0) {
	  mystatus = ACTIVE_ENABLED;
	  status_rampval = 1.0;
	}
	else {
	  status_rampval = (countdown_val-status_count)/((double)countdown_val);
	}
      }
      else {
	// ramping down
	if (status_count == 0) {
	  mystatus = ACTIVE_DISABLED;
	  speed = 0.0;
	  velocity = TRVector(0);
	  status_rampval = 0;
	}
	else {
	  status_rampval = (status_count)/((double)countdown_val);
	}
      }
      // set any new status
      us_data.SetStatus(mystatus);
    }
    

    /* handle a collision, if any */
    if (Coll_pkt->GetCollisionFlag() == COLLISION) {
      collided = 1;      
      TRVector cn = Coll_pkt->GetCollisionNormal();
      TRVector cp = Coll_pkt->GetCollisionPos();

      // just going to not allow 1 m pass wall
      double px = easting - cp.X();
      double py = northing - cp.Y();
      double depth = -(px * cn.X() + py * cn.Y());
      if (depth > 1) {
	// penetrated too far, back them up to be 1 meter in
	easting  += (depth - 1)*cn.X();
	northing += (depth - 1)*cn.Y();
	// update us_data for next get user state
	us_data.SetPosEasting(easting);
	us_data.SetPosNorthing(northing);
      }
      else if (depth < 0) {
	// no longer collided
	collided = 0;
      }
    }
    else {
      collided = 0;
    }
    
    /* display this stuff on the text screen*/
    EraseWindow();

#ifdef RTI_VXWORKS
    // only printing at 1 hz since lack curses, or windows
    if (counter%100 == 0) {
#endif   
    printw2(" Cycle count: %d\n",counter);
    printw3(" Treadport Server Simulator running on %s, port %d \n\n",
		hostname, port);
    printw1(" ------------------------------------- \n");
    printw1(" Keys: \n"); 
    printw2("    %s:  Increase/Decrease Speed \n", speedkey_string);
    printw2("    %s:  Turn \n", turnkey_string);
    printw1("    k            :  kill the treadport.\n");
    printw1("    s            :  start the treadport. \n");
    printw1(" ------------------------------------- \n");    
    printw2(" Status: %s\n",statstr[us_data.GetStatus()]);
    printw2(" Client: %s\n",(server->ClientExists())?("Connected"):("none"));
    printw1(" ------------------------------------- \n");
    printw2(" Treadport Speed:   %7.2lf \n", speed*status_rampval);
    printw3(" Treadport Position:  %7.2lf East %7.2lf North\n",
	   us_data.GetPosEasting(), us_data.GetPosNorthing());
    printw4(" Treadport Facing:    %7.2lf X  %7.2lf Y  %7.2lf Z\n",
	   facing.X(), facing.Y(), facing.Z());
    printw4(" Treadport Velocity:  %7.2lf X  %7.2lf Y  %7.2lf Z\n",
	   velocity.X(),velocity.Y(),velocity.Z());
    printw3(" Treadport Movement & User Scale:   %7.2lf M  %7.2lf U\n", 
	    vs_data.GetMovementScale(), vs_data.GetUserScale());
    printw1(" -------------------------------------\n ");
    
    if (received_virtstate) {
      const char *str;
      printw2(" Terrain: %s\n",
	     terrstr[vs_data.GetTerrain()]);
      str = (collided)?("YEP!"):("");
      printw2(" Collision: %s\n",str);
      TRVector cn = Coll_pkt->GetCollisionNormal();
      TRVector cp = Coll_pkt->GetCollisionPos();

      printw4(" Collision Normal:    %7.2lf X  %7.2lf Y  %7.2lf Z\n", cn.X(),cn.Y(),cn.Z());
      printw4(" Collision Point:     %7.2lf X  %7.2lf Y  %7.2lf Z\n", cp.X(),cp.Y(),cp.Z());

      printw4(" Surface Normal:      %7.2lf X  %7.2lf Y  %7.2lf Z\n",
	     surface_normal.X(),surface_normal.Y(),surface_normal.Z());
      TRVector vf = vs_data.GetVirtualForce();
      printw4(" Virtual Force:       %7.2lf X  %7.2lf Y  %7.2lf Z\n",
	     vf.X(),vf.Y(),vf.Z());
    }
    else {
      printw1(" No Virtual World Info Available. \n");
    }
#ifdef RTI_VXWORKS
    }
#endif

    //#if   defined(sparc_solaris)
    // sleep a bit
    //    struct timespec r;
    //    r.tv_sec = 0;
    //    r.tv_nsec = CYCLE_TIME_NSEC;
    
    //    usleep(r.tv_nsec/1000);
    //#else
    my_sync.block();
    //#endif

    counter++;
  }/* end while */

  StopWindow();

  // Shutdown server
  server->StopServer();
  
  return 1;  
}


