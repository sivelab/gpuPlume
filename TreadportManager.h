#ifndef TREADPORT_MANAGER_H
#define TREADPORT_MANAGER_H 1

/**
 *  TreadportManager.h
 *  Author: Joshua Clark
 *  Created on: 06/15/2009
 *
 *  Adaptation of treadpotManipulator.h (written by Pete Willemsen).
 */

#include <string>
#include <iostream>
#include <vector>
#include "graphicsUtil.h"

//
// Treadport network headers
//
#include "Client_Controller.h"
#include "Command_Packet.h"
#include "ServerRequest_Packet.h"
#include "NewPos_Packet.h"
#include "NewPosResponse_Packet.h"
#include "UserState_Packet.h"
#include "VirtualState_Packet.h"
#include "Collision_Packet.h"

#if defined(_MSC_VER)
  #pragma warning( disable : 4786 )
#endif

/**
 * TreadportManager is the class that is responcible for syncing information with
 * the treadport system.
 */
class TreadportManager {

 public:

  static const std::string HOST;  

  float eye_n;
  float eye_e;

  TreadportManager();

  virtual ~TreadportManager();

  /**
   * Retrieves the current position of the user.
   */
  std::vector<float> getEyePosition();

  /**
   * Retrieves the gaze vector which is a vector pointing
   * at where the user is currently facing / looking.
   */
  std::vector<float> getEyeGaze();

  /**
   * Retrieves the local eye offset needed to calculate
   * the view frustum.
   */
  std::vector<float> getEyeOffset();

  /**
   *  init initializes the Treadport system and sets
   *  the initial position, gazex, (which are defined by
   *  northing and easting values) and surface normal
   *  (in our case up) vectors.
   * 
   *  Note that the vectors taken in here are assuming
   *  the standard OpenGL coordinate system where X is
   *  easting (or side to side), Z is northing (forward
   *  and backward), and Y is up.
   */
  void init(std::vector<float> initial_position, std::vector<float> initial_gaze, std::vector<float> initial_up);

  /**
   *  sync is a modified version of the calcMovement and calcPosition functions from
   *  TreadportManipulator. This function will sync get the current
   *  data from the treadport.
   */
  bool sync();

 protected:

  std::vector<float> eye_pos;

  std::vector<float> eye_gaze;

  std::vector<float> eye_offset;

  // Begin TreadPort Helper Functions.

  int TRDisable();

  int TREnable();

  int TRGetUserState();

  int TRSetPosition(const UserState_Data& usd, const VirtualState_Data& vsd);

  int TRSetVirtualState(VirtualState_Data& vsd);

  // Begin Treadport Data Memebers.
  // Note that these are all ported from treadportManager
  // and some may not be used. TODO: Cleanup is needed.

  float _distance;

  std::string _hostname;
 
  short _port;

  float _rate;

  double _elevation;

  bool m_use_fixed_eyeheight;

  float _eyeheight;

  float m_initial_northing;

  float m_initial_easting;

  Client_Controller* client;

  // User state packs.
  ServerRequest_Packet* RequestUsrState_pkt;
  UserState_Packet* UsrState_pkt;
  int UsrState_type;

  UserState_Data us_data;
  VirtualState_Data vs_data;

  // Virtual state packet.
  VirtualState_Packet* VirtState_pkt;

  // Command packet
  Command_Packet * Cmd_pkt; 

  // New position packets.
  NewPos_Packet * NewPos_pkt;
  NewPosResponse_Packet* NewPosResp_pkt;
  int NewPosResp_type;

  // Collision packet.
  Collision_Packet* Coll_pkt;
  
};

#endif // TREADPORT_MANAGER_H
