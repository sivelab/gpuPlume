#ifndef CLIENT_CONTROLLER_H
#define CLIENT_CONTROLLER_H

#include "Net_Controller.h"
#include "Registration_Packet.h"

// To avoid having to use list structures, I will define
//  the max number of packets that can be sent by the 
//  client controller.  Setting it to 5, which will cover
//  3 data packets, a command, and a request
#define MAX_PACKETS_PER_MSG 5

// Defining a few error values
#define SOCKET_ERR -5
#define PACKET_DEF_CLASH -6
#define OUT_OF_PACKET_SLOTS -7
#define NO_SERVER -8
#define RECEIVE_TIMEOUT -9
#define REGISTRATION_ERROR -10
#define SERVER_REG_UNKNOWN -11
#define INCORRECT_RESPONSE_TAG -12
#define INCORRECT_PACKET_RECEIVED -13

/****************************************************************/
/*  Client_Controller Class (derived from Net_Controller)       */
/*                                                              */
/*  This class handles the client side of the networking.  It   */
/*   connects with the server and sends/recieves packets.       */
/*  I have made the class with virtual methods in case later    */
/*   on a new type of client needs to be created using this one */
/*   as a mold (which is possible for the sound client, or a    */
/*   windows version of the client)                             */
/****************************************************************/

// Definition for only polling the client fds, i.e. not
//  waiting (blocking) for any reply
#define POLL 0

class Client_Controller : public Net_Controller {
 public:
  // Destructor
  virtual ~Client_Controller();

  // Default Constructor
  Client_Controller();

  // Register packets (returns type)
  virtual int Register(Packet * new_packet);

  // Start Client, w/ port and host name passed in
  virtual int ConnectToServer(char * hostname, short port);

  // Stop Client - disconnect from server
  virtual void Disconnect();

  // Add packet for sending
  virtual int AddPacket(Packet * new_packet);

  // Send packets out
  virtual int Send();

  // Receive message - set if wants to wait (based on params)
  //  below, or just a quick check
  virtual int Receive(BOOL wait);

  // Set Receive wait/try params
  inline void SetReceiveParams(int num_tries, int wait_time) {
    m_num_receive_tries = num_tries;
    m_max_wait_time = wait_time;
  }

  // just returns whether a server exists
  inline BOOL ServerExists() {return (BOOL)(m_server_fd>=0);}

  // These are simplified functions for sending/receiving a single packet
  // Result 1 - success, 0 - fail on both
  virtual int SendSinglePacket(Packet *);
  virtual int ReceiveSinglePacket(Packet *);

  /*************************************************************/
  /*   Inherited public functions                              */
  /*                                                           */
  /* int GetPacketType()       - get next packet type from msg */
  /* int ParsePacket(Packet *) - get packet from message       */
  /* int PrintLastError()      - Print/return last error       */
  /* int ReturnLastErrorVal()  - return last error             */
  /* int SetVerbosity()        - Set message print verbosity   */
  /*************************************************************/

 protected:
  // File descriptor for the server connection
  // Only allowing 1 now, multiple connections would be too messy to 
  //  handle and complicate all the routines (like add to which socket)
  // If you want multiple connections, have multiple client_controllers
  int m_server_fd;

  // Array of packets to send.  This is a placeholder for verifying
  //  the responses come back for each packet that requests one
  Packet * m_packets_send[MAX_PACKETS_PER_MSG];
  int m_num_packets_sending;
  int m_num_packs_sent;

  // Saved info for server, so can reconnect if lose
  char * m_server_name;
  short m_server_port;

  // Variables for how many times / how long to wait for
  //  receive command (time in usecs)
  int m_num_receive_tries;
  int m_max_wait_time;

  // Tag for all sent messages to ensure the reply from server
  //  is not from a different message
  ushort m_msg_tag;

  /*************************************************************/
  /*   Variables already accessible in class                   */
  /*                                                           */
  /* int m_tread_verbosity;   -- Verbosity for message prints  */
  /*                                                           */
  /* list of extra error_vals assigned                         */
  /* -301 - socket call failed                                 */
  /* -302 - setsockopt call failed                             */
  /* -303 - connect call failed                                */
  /* -304 - has not specified what server/port to connect to   */
  /* -305 - wrong packet type returned for registration        */
  /* -306 - incorrect registration sent back (class name/size) */
  /* -307 - incorrect registration tag sent back               */
  /*************************************************************/

  // Internal routines

  // Trying to reconnect to server
  virtual int TryReconnectToServer();

  // Internal disconnect - which actually closes socket and cleans up
  virtual void CloseServerConnection();

  // Getting address of computer by its name
  virtual long GetAddrByName(char *name);
};

#endif
