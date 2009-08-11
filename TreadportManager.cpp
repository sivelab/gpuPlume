/**
 *  TreadportManager.cpp
 *  Author: Joshua Clark
 *  Created on: 06/15/2009
 *
 *  Adaptation of treadpotManipulator.cpp (written by Pete Willemsen).
 */

#include "TreadportManager.h"

const std::string TreadportManager::HOST = "roboworks";

TreadportManager::TreadportManager() : 
  _hostname(HOST), 
  _port(5053), 
  _rate(1.0),
  _elevation(0.0),
  m_use_fixed_eyeheight(false),
  _eyeheight(1.8) {

  // Initialize some variables.
  eye_n = 0.0;
  eye_e = 0.0;
  eye_pos.resize(3);
  eye_gaze.resize(3);
  eye_offset.resize(3);

  // Create client controller and connect to server
  client = new Client_Controller();

  switch (client->ConnectToServer(const_cast<char*>(_hostname.c_str()),_port)) {
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
    fprintf(stderr, "Connection successful to %s:%d\n", _hostname.c_str(), _port);
    break;
  }
}

TreadportManager::~TreadportManager() {
  // Nothing to do.
}

std::vector<float> TreadportManager::getEyePosition() {
  return eye_pos;
}

std::vector<float> TreadportManager::getEyeGaze() {
  return eye_gaze;
}

std::vector<float> TreadportManager::getEyeOffset() {
  return eye_offset;
}

void TreadportManager::init(std::vector<float> initial_position, std::vector<float> initial_gaze, std::vector<float> initial_up) {

  //
  // Create all the packets to communicate with the treadport system.
  // Note that their is some limited documentation for these packets
  // within the header files for TreadNetLib.
  //
  RequestUsrState_pkt = new ServerRequest_Packet();
  UsrState_pkt = new UserState_Packet();
  VirtState_pkt = new VirtualState_Packet();
  Cmd_pkt = new Command_Packet();
  NewPos_pkt = new NewPos_Packet();
  NewPosResp_pkt = new NewPosResponse_Packet();
  Coll_pkt = new Collision_Packet();

  //
  // Register all packets (that we are going to send/accpet) 
  // with the treadport server.
  //

  // Request and Response for User State
  UsrState_type = client->Register(UsrState_pkt);
  if(UsrState_type < 1) {
    // Error in Register.
    std::cerr << "Error in register." << std::endl;
    client->PrintLastError();
    exit(1);
  }
  std::cout << "UserState_Packet type is " << UsrState_type << std::endl;

  // Note, we are not saving this type, because this is something we send
  // and thus we don't care about the type
  int result = client->Register(RequestUsrState_pkt);
  if(result < 1) {
    // Error in Register.
    client->PrintLastError();
    exit(1);
  }

  if(RequestUsrState_pkt->AddPacketRequest(UsrState_pkt) == -1) {
    std::cerr << "Some Error Message: Couldn't add request." <<  std::endl;
  }

  // Virtual State
  result = client->Register(VirtState_pkt);
  if(result < 1) {
    client->PrintLastError();
    exit(1);
  }

  // Command Packet
  result = client->Register(Cmd_pkt);
  if(result < 1) {
    client->PrintLastError();
    exit(1);
  }

  // New Position Request and Response
  result = client->Register(NewPos_pkt);
  if(result < 1) {
    client->PrintLastError();
    exit(1);
  }
  NewPosResp_type = client->Register(NewPosResp_pkt);
  if(NewPosResp_type < 1) {
    client->PrintLastError();
    exit(1);
  }
  std::cout << "New Position Response Packet type is " << NewPosResp_type << std::endl;

  // Collision
  result = client->Register(Coll_pkt);
  if(result < 1) {
    client->PrintLastError();
    exit(1);
  }

  std::cout << "Registered all packets." << std::endl;
  
  //
  // Done registering all packets.
  //


  //
  // Set initial values for the treadport. Note that you may have to do some conversion
  // from your program's coordinate system into the treadport coordinate system.
  //
  Vector4 te_start = Vector4(initial_position[0], initial_position[1], initial_position[2], 1.0);
  Vector4 te_facing = Vector4(initial_gaze[0], initial_gaze[1], initial_gaze[2], 1.0);
  Vector4 te_normal = Vector4(initial_up[0], initial_up[1], initial_up[2], 1.0);

  // Note that currently no conversion matrix is being applied, if you need to apply a 
  // rotation matrix to put things into the right coordinate system, this is where it
  // would be done.
  Vector4 tr_start = te_start;
  Vector4 tr_facing = te_facing;
  Vector4 tr_normal = te_normal;
 
  // Place the values into UserState_Data so that they can be passed to the treadport system.
  us_data.SetPosEasting(tr_start[0]);
  us_data.SetPosNorthing(tr_start[1]);  
  us_data.SetFacingEasting(tr_facing[0]);
  us_data.SetFacingNorthing(tr_facing[1]);
 
  // Set some default values, these values were taken from the terrain_render application, 
  // and may need to be changed (if needed).
  TRVector surf_normal(tr_normal[0], tr_normal[1], tr_normal[2]);
  TRTerrain terrain = NORMAL_TERRAIN;
  TRVector virtual_force(0);
  double user_scale = 1;
  _rate = 1.0;

  // Place values into VirtualState_Data so that they can be passed to the treadport system.
  vs_data.SetSurfaceNormal(surf_normal);
  vs_data.SetTerrain(terrain);
  vs_data.SetVirtualForce(virtual_force);
  vs_data.SetMovementScale(_rate);
  vs_data.SetUserScale(user_scale);


  //
  // Done setting initial values.
  //

  //
  // Begin passing the initial values to the treadport system.
  //

  // Disable the Treadport; this must be done in order for the new 
  // position to work.
  TRDisable();
  sleep(1);

  // Because this must be done, we have to check and wait for it to
  // be disabled (i.e. make sure).
  TRGetUserState();
  while(UsrState_pkt->GetStatus() & ENABLED_BIT) {
    std::cout << "Awaiting disabled mode on the treadport..." << std::endl;
    sleep(2);
    TRGetUserState();
  }

  // Start the user in their initial position (i.e. send out
  // intial values to the server).
  TRSetPosition(us_data, vs_data);
  
  // Renable the treadport.
  TREnable();

  //
  // Done passing initial values to the treadport
  //
  std::cout << "Finished sending initial values to the treadport." << std::endl;

}

bool TreadportManager::sync() {
  // Request the latest position and status from the treadport.
  TRGetUserState();

  // If the treadport enters a disabled state whether in active or
  // inactive mode, re-enable the graphics system.
  if(UsrState_pkt->GetStatus() & DISABLED_BIT) {
    TREnable();
  }

  // Check if valid position returned
  if(!UsrState_pkt->PosValid()) {
    // Invalid position returned, which means we lost the server and
    // when it rebooted it no longer has a valid position, set it.
    TRDisable();
    TRSetPosition(us_data, vs_data);
  } else {
    // New data received, get copy of it
    UsrState_pkt->GetCopyOfData(&us_data);
  }

  //
  // Begin setting values that were recieved from the treadport.
  //

  TRVector facing = us_data.GetFacing();

  // Note, currently we are not using velocity.
  // TRVector velocity = us_data.GetVelocity();

  // Store the eye height for frustum calculations.
  // Note using a dynamic eye height (for the frustum calculation
  // may produce too much movement due to head bobbing.
  // 
  // Also, for the time being we are assuming a flat surface that is
  // at a zero height.
  _eyeheight = us_data.GetEyeElevation();
  
  Vector4 tr_current = Vector4(us_data.GetPosEasting() + us_data.GetEyeEasting(),
			       us_data.GetPosNorthing() + us_data.GetEyeNorthing(),
			       0.0, // Note, we are factoring the height value in the eye_offset.
			       1.0);

  // Note, this is where some conversion would happen if we were applying
  // a rotation matrix to undo the initial rotation when we set our initial
  // values (within the init function).
  Vector4 te_current = tr_current;

  // Anything done with the height field would be done here, see treadportManipulator.cpp...
  // For the time being the gpuPlume system is running on a flat surface which means that
  // the surface normal is static (always point up).

  // IMPORTANT: Make sure you have your coordinate system correct. The treadport system
  // does not know how to handle a surface normal that is perpendicular to the up direction
  // and you will get a goofy NaN (not a number) error.
  
  /*
  TRVector surface_normal(0.0, 0.0, 1.0);
  TRTerrain terrain = NORMAL_TERRAIN;
  TRVector virtual_force(0);
  double user_scale = 1;
  _rate = 1.0;

  vs_data.SetSurfaceNormal(surface_normal);
  vs_data.SetTerrain(terrain);
  vs_data.SetVirtualForce(virtual_force);
  vs_data.SetMovementScale(_rate);
  vs_data.SetUserScale(user_scale);
  */

  // Sync the virtual state, including the new normal (but since our normal is currently
  // fixed, there is no need to recalulate it.
  TRSetVirtualState(vs_data);

  // Get the look at direction.
  std::vector<float> tr_facing;
  tr_facing.resize(3);
  tr_facing[0] = facing.X();
  tr_facing[1] = facing.Y();
  tr_facing[2] = facing.Z();
  
  // Again, we aren't applying an conversion matrix to get things back into
  // our coordinate system.
  std::vector<float> te_facing = tr_facing;

  // Set the sync'd values... (this was done in the computePosition method).
  eye_pos[0] = te_current[0];
  eye_pos[1] = te_current[1];
  eye_pos[2] = te_current[2];

  // Note that eye_gaze may need to be normalized...
  eye_gaze[0] = te_facing[0];
  eye_gaze[1] = te_facing[1];
  eye_gaze[2] = te_facing[2];
  
  // set the eye offset in the viewer... so the master can pass the
  // eye offset to the slaves and thus compute the correct frustum
  double eyeN = us_data.GetEyeNorthing();
  double eyeE = us_data.GetEyeEasting();
  
  // ////////////////////////////////////////////
  // uncomment these to test eye offset input on frustum code
  // u, j, k, and m control movement
  // eyeN =  eye_n;
  // eyeE = eye_e;
  //
  // eyeN = sin( val );
  // eyeE = sin( val );
  // double eye_Z = sin( val ); val += 0.01;
  // ////////////////////////////////////////////

  // The dynamic frustum moves too much when using the user eye
  // point offset so, this gives us the chance to restrict the
  // movement.
  //
  // Set the user Northing (direction forward on belt) to 0.0 or,
  // ideally, the offset forward in which the person walks, maybe
  // 0.10m or so??
  // eyeN = 0.0;

  // Set the user Easting position (direction to the right and left
  // of where the user starts on the belt).  Set this to 0.0;
  // eyeE = 0.0;

  double facingN = facing.Y();
  double facingE = facing.X();

  double local_X = eyeN * facingN + eyeE * facingE;  // ie. Easting
  double local_Y = eyeN * facingE - eyeE * facingN;  // ie. Northing

  // set the offset in the viewer --- add in the frustum code....
  eye_offset[0] = -local_Y;
  eye_offset[1] = local_X;
  eye_offset[2] = _eyeheight;

  return true;
}

//
// Begin Treadport Utility Functions.
//

int TreadportManager::TRDisable() {
  Cmd_pkt->SetDisableCmd();
  if (client->SendSinglePacket(Cmd_pkt)) {
    return 1;
  }
  // else error
  client->PrintLastError();
  return 0;
}

int TreadportManager::TREnable() {
  Cmd_pkt->SetEnableCmd();
  if (client->SendSinglePacket(Cmd_pkt)) {
    return 1;
  }
  // else error
  client->PrintLastError();
  return 0;
}

int TreadportManager::TRGetUserState() {
  if (client->SendSinglePacket(RequestUsrState_pkt) &&
      client->ReceiveSinglePacket(UsrState_pkt)) {
    // both worked correctly
    return 1;
  }
  // otherwise some error, print
  client->PrintLastError();
  return 0;
}

int TreadportManager::TRSetPosition(const UserState_Data& usd, const VirtualState_Data& vsd) {
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

int TreadportManager::TRSetVirtualState(VirtualState_Data& vsd) {
  VirtState_pkt->SetData(&vsd);
  if (client->SendSinglePacket(VirtState_pkt)) {
    return 1;
  }
  // else error
  client->PrintLastError();
  return 0;  
}

//
// End Treadport Utility Functions.
//
