/**
 * NetworkManager.cpp
 *
 * Author: Joshua Clark
 * Created On: Tuesday June 2nd, 2009
 */

#include "NetworkManager.h"

NetworkManager::NetworkManager() {

  // Set the default mode.
  mode = NetworkManager::DISABLED;

  // Set the default port number.
  portNumber = 8101;

  // Set the default broadcast address.
  // TODO: Read from config file.
  // broadcastAddr = "192.168.100.255"; // Local Bcast Address for VR Lab (UMD)
  broadcastAddr = "172.19.1.255"; // Local Bcast Address for TPAWT (U of U)

  // Set the default device list.
  deviceList.resize(5);
  deviceList[0] = "eth0";
  deviceList[1] = "eth1";
  deviceList[2] = "eth2";
  deviceList[3] = "eth3";
  deviceList[4] = "lo";

  // Initialize frameID.
  frameID = 0;
}

NetworkManager::~NetworkManager() {
  // Nothing to do.
}

void NetworkManager::init() {

  if(mode == BROADCAST) {

    //
    // Find an adapter to broadcast on.
    //

    for(int i = 0; i < (int)deviceList.size(); i++) {
      std::cout << "NetworkManager: Checking if " << deviceList[i]
                << " has broadcast adddress " << broadcastAddr
                << " ..." << std::endl;

      // Lookup broadcast address for the network devices in deviceList
      struct ifreq ifr;
      int sock;
      in_addr_t targetAddr = inet_addr( broadcastAddr.c_str() );
      strncpy( ifr.ifr_name, deviceList[i].c_str(), sizeof(ifr.ifr_name) );

      if( (sock = socket( AF_INET, SOCK_DGRAM, 0 )) < 0 ) {
        perror( "NetworkManager: Socket" );
        exit(1);
      }

      if( (ioctl( sock, SIOCGIFBRDADDR, &ifr)) < 0 ) {
        perror("NetworkManager: Cannot get broadcast device");
        continue;
      }

      // we now have the broadcast address for the device
      in_addr_t deviceAddr = (((sockaddr_in *)&ifr.ifr_broadaddr)->sin_addr.s_addr);

      // check if it is the broadcast address we were looking for
      if(memcmp(&targetAddr, &deviceAddr, sizeof(in_addr_t)) == 0) {

        unsigned char *ptr = (unsigned char *)&deviceAddr;

        std::cout << "NetworkManager: SUCCESS! Broadcasting to "
                  << (unsigned int)ptr[0] << "."
                  << (unsigned int)ptr[1] << "."
                  << (unsigned int)ptr[2] << "."
                  << (unsigned int)ptr[3]
                  << " (with " << deviceList[i] << ")." << std::endl;

        /* FINALLY!  Set up the broadcaster! */
        bc.setEthernetDevice( deviceList[i] );
        bc.setPort( portNumber );
        break;

      } else {
        // if it wasn't the address we were looking for, continue trying...
        unsigned char *ptr = (unsigned char *)&deviceAddr;

        std::cout << "NetworkManager: " << deviceList[i]
                  << " has broadcast address "
                  << (unsigned int)ptr[0] << "."
                  << (unsigned int)ptr[1] << "."
                  << (unsigned int)ptr[2] << "."
                  << (unsigned int)ptr[3]
                  << std::endl;
      } // end if
    } // end for

    // fail if we can not find the broadcast address we were looking for
    // osg::notify(osg::FATAL) << "Settings: Can't find network device with broadcast address " << broadcastAddr << std::endl;
    // exit(1);

  } else if(mode == RECEIVE) {

    //
    // Set the listing network address.
    //
    std::cout << "Settings: Listening to network address " << broadcastAddr << std::endl;
    rc.setHost(broadcastAddr);
    rc.setPort(portNumber);
  }

}

void NetworkManager::sync(NetworkSyncData &data) {

  if(mode == BROADCAST) {
    data.frame_id = ++frameID;
    bc.setBuffer(&data, sizeof(data));
    usleep(1000);
    bc.sync();
  } else if(mode == RECEIVE) {
    NetworkSyncData newData;
    rc.setBuffer(&newData, sizeof(newData));
    rc.sync();
    if(newData.frame_id != (frameID + 1)) {
      std::cerr << "NetworkManager: Data Frame Skip!\n\tReceived Frame ID = "
                << newData.frame_id << ", Last Received Frame ID = "
                << frameID << std::endl;
    }
    frameID = newData.frame_id;
    data = newData;
  }

  // Code to grab data, this should be in display control.
#if 0
  if (gbl_net_mode == 0)
    {
      memcpy(curr_sharedData.eye_pos, eye_pos, sizeof(float) * 3);
      memcpy(curr_sharedData.eye_gaze, eye_gaze, sizeof(float) * 3);
      memcpy(curr_sharedData.eye_up, eye_up, sizeof(float) * 3);
      memcpy(&curr_sharedData.wim_scale, &wim_scale_factor, sizeof(float));
      memcpy(&curr_sharedData.gesture_id, &gesture_id, sizeof(int));
      memcpy(curr_sharedData.hand_pos, current_hand_pos, sizeof(float) * 3);
      memcpy(curr_sharedData.hand_up, current_hand_up, sizeof(float) * 3);
    }
  else
    {
      memcpy(&prev_sharedData, &curr_sharedData, sizeof(NetworkSyncData));;
      memcpy(eye_pos, curr_sharedData.eye_pos, sizeof(float) * 3);
      memcpy(eye_gaze, curr_sharedData.eye_gaze, sizeof(float) * 3);
      memcpy(eye_up, curr_sharedData.eye_up, sizeof(float) * 3);
      memcpy(&wim_scale_factor, &curr_sharedData.wim_scale, sizeof(float));
      memcpy(&gesture_id, &curr_sharedData.gesture_id, sizeof(int));
      memcpy(current_hand_pos, curr_sharedData.hand_pos, sizeof(float) * 3);
      memcpy(current_hand_up, curr_sharedData.hand_up, sizeof(float) * 3);
      // memcpy(curr_sharedData.hand_matrix, hand_matrix, sizeof(float) * 16);
    }
#endif

}

NetworkManager::Mode NetworkManager::getMode() {
  return mode;
}

void NetworkManager::setMode(NetworkManager::Mode newMode) {
  mode = newMode;
}

int NetworkManager::getPortNumber() {
  return portNumber;
}

void NetworkManager::setPortNumber(int newPortNumber) {
  portNumber = newPortNumber;
}

std::string NetworkManager::getBroadcastAddr() {
  return broadcastAddr;
}

void NetworkManager::setBroadcastAddr(std::string newAddr) {
  broadcastAddr = newAddr;
}

std::vector<std::string> NetworkManager::getDeviceList() {
  return deviceList;
}

void NetworkManager::setDeviceList(std::vector<std::string> newDeviceList) {
  if(newDeviceList.size() == 0) {
    std::cerr << "NetworkManager.setDeviceList: WARNING: "
              << "newDeviceList has a size of 0." << std::endl;
  }

  // Resize to fit the new data.
  deviceList.resize(newDeviceList.size());

  // Copy all of the new data.
  for(int i = 0; i < (int)newDeviceList.size(); i++) {
    deviceList[i] = newDeviceList[i];
  }
}
