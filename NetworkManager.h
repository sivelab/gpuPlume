#ifndef NETWORK_MANAGER_H
#define NETWORK_MANAGER_H

/**
 * NetworkManager.h
 *
 * Author: Joshua Clark
 * Created On: Tuesday June 2nd, 2009
 */

#include <arpa/inet.h>
#include <net/if.h>
#include <sys/socket.h>
#include <sys/ioctl.h>

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <unistd.h>

#include "broadcaster.h"
#include "receiver.h"

/**
 * NetworkManager is responcible for managing the network synchronization code
 * within gpuPlume.
 */
class NetworkManager {

public:

  /**
   * Mode is an enum describing the state of the network manager, and thus the
   * network state of gpuPlume. Each NetworkManager object has it's own mode.
   */
  enum Mode {
    DISABLED = -1,
    BROADCAST = 0,
    RECEIVE = 1
  };

  /**
   * NetworkSyncData is a struct representing the data that will be sync'd
   * across the network.
   */
  struct NetworkSyncData {
    long frame_id;
    float eye_pos[3];
    float eye_gaze[3];
    float eye_offset[3];
    float eye_up[3];
    float eye_right[3];
    int gesture_id;
    float wim_scale;
    bool emit_particles;
    float part_emit_pos[3];
    float hand_pos[3];
    float hand_up[3];
    double hand_matrix[16];
    bool draw_layers;
    bool inPauseMode;
  };

  /**
   * Default constructor.
   */
  NetworkManager();

  /**
   * Destructor.
   */
  ~NetworkManager();

  /**
   * init will initialize the NetworkManager. Note, the proper mode must be set
   * before init is called, otherwise the default mode is assumed.
   */
  void init();

  /**
   * sync will synchronize the given NetworkSyncData across the network. Note,
   * initialize must be called before this. This method will either receive or
   * broadcast depending upon the mode which was previously set; if DISABLED is
   * set then this method does nothing.
   *
   * Note, that this method will sleep for 1000 miliseconds before actually
   * doing the synch over the network.
   */
  void sync(NetworkSyncData &data);

  /**
   * getMode is the accessor for mode.
   */
  NetworkManager::Mode getMode();

  /**
   * setMode is the mutator for mode.
   */
  void setMode(NetworkManager::Mode newMode);

  /**
   * getPortNumber is the accessor for portNumber.
   */
  int getPortNumber();

  /**
   * setPortNumber is the mutator for portNumber, note that the default value
   * is 8101.
   */
  void setPortNumber(int newPortNumber);

  /**
   * getBroadcastAddr is the accessor for broadcastAddr.
   */
  std::string getBroadcastAddr();

  /**
   * setBroadcastAddr is the mutator for broadcastAddr.
   */
  void setBroadcastAddr(std::string newAddr);

  /**
   * getDeviceList is the accessor for deviceList.
   */
  std::vector<std::string> getDeviceList();

  /**
   * setDeviceList is the mutator for deviceList.
   */
  void setDeviceList(std::vector<std::string> newDeviceList);

protected:

  /**
   * mode is an enum describing the state of the network manager, and thus the
   * network state of gpuPlume. The default state should be DISABLED or what has
   * been set within the input.txt config file.
   */
  Mode mode;

  /**
   * portNumber is the network port on which the broadcaster and receiver will
   * operate. The default value is 8101.
   */
  int portNumber;

  /**
   * broadcastAddr is the broadcast address. The default value is
   * 192.168.100.255.
   */
  std::string broadcastAddr;

  /**
   * deviceList is a listing of all of the devices to check to receive or
   * broadcast on.
   *
   * We are given a broadcast address, we need to find a device
   * that broadcasts to the given address (This is a roundabout
   * way of doing things, but SUSE doesn't seem to assign
   * network devices to the same cards at boot.)
   *
   * Broadcasting on the loopback device is allowed for testing.
   */
  std::vector<std::string> deviceList;

  /**
   * frameID is the frame ID of the frame that was most recently sent/received.
   */
  int frameID;

  /**
   * bc is the interface to directly sending the data over the network.
   */
  Broadcaster bc;

  /**
   * rc is the interface to directly receiving data over the network.
   */
  Receiver rc;

private:

};

#endif // NETWORK_MANAGER_H
