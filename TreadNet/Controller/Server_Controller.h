#ifndef SERVER_CONTROLLER_H
#define SERVER_CONTROLLER_H

#include "Net_Controller.h"
#include "Registration_Packet.h"

// Setting max types so I can avoid using lists.  We only
//  have 10 currently (some not even created yet), so this
//  should be enough
#define MAX_PACKET_TYPES 15

// To avoid having to use list structures, I will define
//  the max number of packets that can be sent by the 
//  server controller.  Server really does not need many since
//  it mostly just sends a single data packet, plus a response or
//  two, so 5 should cover it.
#define MAX_PACKETS_PER_MSG 5

// Defining a few error values
#define SOCKET_ERR -5
#define PACKET_DEF_CLASH -6
#define OUT_OF_PACKET_TYPES -7
#define OUT_OF_PACKET_SLOTS -8
#define NO_CLIENT -9
#define REGISTRATION_ERROR -10

/****************************************************************/
/*  Server_Controller Class (derived from Net_Controller)       */
/*                                                              */
/*  This class handles the server side of the networking.  It   */
/*   creates the server network connect, accepts clients, and   */
/*   receives/sends packets.                                    */
/*  I have made the class with virtual methods in case later    */
/*   on a new type of server needs to be created using this one */
/*   as a mold (which is possible for the sound server, or a    */
/*   windows version of the server)                             */
/****************************************************************/

// Definition for only polling the client fds, i.e. not
//  waiting (blocking) for any reply
#define POLL 0

class Server_Controller : public Net_Controller {
 public:
  // Destructor
  virtual ~Server_Controller();

  // Default Constructor
  Server_Controller();

  // Register packets (returns type)
  virtual int Register(Packet * new_packet);

  // Start Server, w/ port passed in
  virtual int StartServer(short port);

  // Stop Server
  virtual void StopServer();

  // Add packet for sending
  virtual int AddPacket(Packet * new_packet);

  // Send packets out
  virtual int Send();

  // Receive message (takes care of most of TRNet) - never
  //  waits for a message...
  virtual int Receive(BOOL wait = FALSE);

  // Disconnect client
  virtual void DisconnectClient();

  // just resturns whether a client exists
  inline BOOL ClientExists() {return (BOOL)(m_client_fd>=0);}

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
  // File descriptor for the port listening for client connections
  int m_listening_fd;

  // File descriptor for the client connection
  // Only allowing 1 now, multiple connections would be too messy to 
  //  handle and complicate all the routines (like add to which socket)
  // If you want multiple connections, have multiple server_controllers
  int m_client_fd;

  // Array of packets to send.  This is a placeholder for calling
  //  PutPacketInMsg when ready for Send()
  Packet * m_packets_send[MAX_PACKETS_PER_MSG];
  int m_num_packets_sending;

  // Last message tag, so know what to send back w/ response
  ushort m_last_msg_tag;

  // Variables for holding registered packet classes
  char *m_reg_packet_names[MAX_PACKET_TYPES];
  uint m_packet_types[MAX_PACKET_TYPES];
  ushort m_packet_sizes[MAX_PACKET_TYPES];
  int m_num_packets_reg;

  /*************************************************************/
  /*   Variables already accessible in class                   */
  /*                                                           */
  /* int m_tread_verbosity;   -- Verbosity for message prints  */
  /*                                                           */
  /* list of extra error_vals assigned                         */
  /* -201 - socket call failed                                 */
  /* -202 - setsockopt call failed                             */
  /* -203 - bind call failed                                   */
  /* -204 - listen call failed                                 */
  /* -205 - accept call failed                                 */
  /* -206 - received unknown reg packet                        */
  /* -207 - received known, but invalid size reg packet        */
  /*************************************************************/

  // Internal routines
  // Check if valid registration
  virtual int CheckValidRegistration();

  // Check for new client
  virtual int CheckForNewClient();

};

#endif
