#include "Server_Controller.h"

/***********************************************************/
/*  Name : Destructor                                      */
/*                                                         */
/*  Description: Destructor                                */
/*                                                         */
/*  Input : none                                           */
/*  Return Value : None                                    */
/***********************************************************/

Server_Controller::~Server_Controller() {
  // must use free to delete the string (since used strdup)
  int i;
  for (i = 0; i < m_num_packets_reg; i++) {
    free(m_reg_packet_names[i]);
  }
}

/***********************************************************/
/*  Name : Server_Controller                               */
/*                                                         */
/*  Description: Default Constructor                       */
/*                                                         */
/*  Input : none                                           */
/*  Return Value : None                                    */
/***********************************************************/

Server_Controller::Server_Controller() : Net_Controller() {
  // no socket connections yet, so set file descriptors to -1
  m_listening_fd = -1;
  m_client_fd = -1;
  
  m_tread_verbosity = 0;
  m_num_packets_reg = 0;
  m_last_msg_tag = 0;

  int i;
  for (i=0; i < MAX_PACKET_TYPES; i++) {
    m_reg_packet_names[i] = NULL;
    m_packet_types[i] = 0;
    m_packet_sizes[i] = 0;
  }

  for (i=0; i < MAX_PACKETS_PER_MSG; i++) {
    m_packets_send[i] = NULL;
  }
  m_num_packets_sending = 0;
}

/***********************************************************/
/*  Name : Register                                        */
/*                                                         */
/*  Description: Register Packet with Server               */
/*    Essentially this specifies what packets the server   */
/*     will recognize, and what their size is (could add   */
/*     some scheme for checking data sizes also!!!)        */
/*                                                         */
/*  Input : Packet ptr                                     */
/*  Return Value : returns Packet type (if successful)     */
/*                 PACKET_ERROR - Invalid packet/packet    */
/*                                  class name             */
/*                 PACKET_DEF_CLASH - if already registered*/
/*                        this class but it was a different*/
/*                        size or type                     */
/*                 OUT_OF_PACKET_TYPES - ran out of types  */
/*                        (increase MAX_PACKET_TYPES)      */
/***********************************************************/

int Server_Controller::Register(Packet * new_packet) {
  if (new_packet == NULL) {
    SET_RETURN_ERROR(PACKET_ERROR, PACKET_ERROR, "Register",
		     "Packet was Null");
  }
  const char * class_name = new_packet->ClassName();
  if (class_name == NULL) {
    SET_RETURN_ERROR(PACKET_ERROR, PACKET_ERROR, "Register",
		     "Packet class name was NULL");
  }

  // Check if already in registry
  int i;
  for (i=0; i < m_num_packets_reg; i++) {
    if (strcmp(m_reg_packet_names[i],class_name) == 0) {
      // found it
      if (new_packet->IsRegistered()) {
	// if already registered, ensure that it matches
	//  what we found in the registry
	if ((new_packet->GetType() != m_packet_types[i]) ||
	    (new_packet->GetSize() != m_packet_sizes[i])) {
	  SET_RETURN_ERROR(PACKET_DEF_CLASH, PACKET_DEF_CLASH, "Register",
			   "Packet registration size and type does not match input packet");
	}
	else {
	  // matches
	  return m_packet_types[i];
	}
      }
      else {
	// Not already registered, set packet registered
	SetPacketType(new_packet,m_packet_types[i]);
	SetPacketRegistered(new_packet);
	return m_packet_types[i];
      }
    }
  }

  // otherwise this is a new packet, create a type for it
  // Currently I am just incrementing (first packet is type
  //  1, then 2, etc.).  Really do not need m_packet_types,
  //  but will keep for now so easier to adjust this scheme
  //  later on.
  if (m_num_packets_reg >= MAX_PACKET_TYPES) {
    // out of types already
    SET_RETURN_ERROR(OUT_OF_PACKET_TYPES, OUT_OF_PACKET_TYPES, "Register",
		     "Ran out of Packet Types, increase MAX_PACKET_TYPES");
  }
  else {
    // Create and Store new registry
    uint new_type = m_num_packets_reg + 1;
    m_reg_packet_names[m_num_packets_reg] = strdup(class_name);
    m_packet_types[m_num_packets_reg] = new_type;
    m_packet_sizes[m_num_packets_reg] = new_packet->GetSize();
    m_num_packets_reg++;

    // Tell packet it is registered
    SetPacketType(new_packet,new_type);
    SetPacketRegistered(new_packet);
  }
  return m_num_packets_reg;
}

/***********************************************************/
/*  Name : StartServer                                     */
/*                                                         */
/*  Description: Start TCP server on the specified port    */
/*                                                         */
/*  Input : port number                                    */
/*  Return Value : if successful (1)                       */
/*                 if failed (SOCKET_ERR)                  */
/***********************************************************/

int Server_Controller::StartServer(short port) {
  // If server already started, return success
  if (m_listening_fd != -1) return 1;
  
#ifdef WIN32
  WORD wVersionRequested = MAKEWORD( 2, 0 );
  WSADATA wsaData;
  int err = WSAStartup( wVersionRequested, &wsaData ); 
  if (err != 0) {
    // Tell the user that we couldn't find a usable 
    // WinSock DLL.
    SET_RETURN_ERROR(SOCKET_ERR, SOCKET_ERR, "ConnectToServer",
		     "Couldn't hi find a useable WinSock DLL.");
  }
#endif

  // Initialize and set server socket address information
  // Using TCP/IP sock streams
  int sockAddrSize = sizeof(struct sockaddr_in);
  struct sockaddr_in server_addr;

  memset((char *)&server_addr, 0, sockAddrSize);
  server_addr.sin_family = AF_INET;
  server_addr.sin_port = htons((unsigned short) port);
  server_addr.sin_addr.s_addr = htonl(INADDR_ANY);

  // Create socket
  if ((m_listening_fd = socket(AF_INET, SOCK_STREAM, 0)) == -1) {
    SET_RETURN_ERROR(SOCKET_ERR, -201, "StartServer", 
		     "Socket call failed");
  }

  // Set socket options so can reuse the port immediately and
  //  to TCP nodelay (i.e. send packets as soon as possible)
  int on = 1;

  if ((setsockopt(m_listening_fd, SOL_SOCKET, SO_REUSEADDR, (char *)&on, sizeof(int)) == -1) ||
      (setsockopt(m_listening_fd, IPPROTO_TCP, TCP_NODELAY, (char *)&on, sizeof(int)) == -1)) {
    close(m_listening_fd);
    m_listening_fd = -1;
    if (m_tread_verbosity) {
      perror("Server_Controller: StartServer (setsockopt)");
    } 
    SET_RETURN_ERROR(SOCKET_ERR, -202, "StartServer",
		     "SetSockOpt call failed");
  }

  // Bind the socket to the server address
  if (bind(m_listening_fd, (struct sockaddr *)&server_addr, sockAddrSize) == -1) {
    close(m_listening_fd);
    m_listening_fd = -1;
    SET_RETURN_ERROR(SOCKET_ERR, -203, "StartServer",
		     "Bind call failed");
  }

  // only allow one queued connection
  if (listen(m_listening_fd, 1) == -1) {
    close(m_listening_fd);
    m_listening_fd = -1;
    SET_RETURN_ERROR(SOCKET_ERR, -204, "StartServer",
		     "Listen call failed");
  }

  if (m_tread_verbosity) 
    fprintf(stderr, "Server_Controller (StartServer) : Server established on port %d.\n", port);

  Cancel_Message();
  m_num_packets_sending = 0;

  return 1;
}

/***********************************************************/
/*  Name : StopServer                                      */
/*                                                         */
/*  Description: Stop TCP server (close all sockets), plus */
/*                remove any pending message               */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

void Server_Controller::StopServer() {
  if (m_client_fd >= 0) {
    // if client socket exists
    close(m_client_fd);
    m_client_fd = -1;
  }
  
  if (m_listening_fd >= 0) {
    // if listening socket exists
    close(m_listening_fd);
    m_listening_fd = -1;
  }

  if (m_tread_verbosity) 
    fprintf(stderr, "Server_Controller (StopServer) : Server stopped\n");

  Cancel_Message();
  m_num_packets_sending = 0;

#ifdef WIN32
  WSACleanup();
#endif
}

/***********************************************************/
/*  Name : AddPacket                                       */
/*                                                         */
/*  Description: Add packet routine.  This implements the  */
/*    server protocol - packets are placemarked for adding */
/*    to message, but not actually added till Send.        */
/*                                                         */
/*  Notes: 1) Consequence of this is the packet''s data can*/
/*     can be altered in between AddPacket and Send.       */
/*         2) This was the easiest way to avoid clashing   */
/*     of message buffers (w/ server AddPacket is a likely */
/*     call while Parsing the received message, thus both  */
/*     needing the message buffer).                        */
/*                                                         */
/*  Input : Packet Ptr                                     */
/*  Return Value : 1 - if successful                       */
/*                 PACKET_ERROR - if packet is invalid     */
/*                       null or unregistered)             */
/*                 OUT_OF_PACKET_SLOTS - if too many       */
/*                       packets in message                */
/***********************************************************/

int Server_Controller::AddPacket(Packet * new_packet) {
  if (m_num_packets_sending >= MAX_PACKETS_PER_MSG) {
    SET_RETURN_ERROR(OUT_OF_PACKET_SLOTS, OUT_OF_PACKET_SLOTS, "AddPacket",
		     "Too many packets being sent - increase MAX_PACKETS_PER_MSG");
  }

  if (new_packet != NULL && new_packet->IsRegistered()) {
    m_packets_send[m_num_packets_sending] = new_packet;
    m_num_packets_sending++;
  }
  else {
    SET_RETURN_ERROR(PACKET_ERROR, -103, "AddPacket",
		     "Unregistered or Null packet");
  }
  return 1;
}

/***********************************************************/
/*  Name : Receive                                         */
/*                                                         */
/*  Description: Main Server routine.  It checks for       */
/*    incoming messages from the client and implements     */
/*    the server protocol (like can not receive until old  */
/*    message has been fully parsed)                       */
/*       It also watches for registration packets, which it*/
/*    handles here (one packet per msg, and nothing else is*/
/*    the protocol for registration packets).              */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : 1 if received some message              */
/*                 NO_MSG if no message                    */
/*                 MSG_BUF_BUSY if have not finished       */
/*                    parsing last msg                     */
/*                 MSG_ERROR if error in read              */
/*                 NO_CLIENT if no client exists           */
/*                 REGISTRATION_ERROR if was a registration*/
/*                    packet, but was invalid              */
/***********************************************************/

int Server_Controller::Receive(BOOL wait) {
  if (MessageInUse()) {
    SET_RETURN_ERROR(MSG_BUF_BUSY, MSG_BUF_BUSY, "Receive",
		     "Message Buffer Busy");
  }

  if (m_client_fd >= 0) {
    // receive message since client exists
    int rec_result = Receive_Message(m_client_fd, POLL);
    switch(rec_result) {
    case NO_MSG:
      // No message was received, nothing to do
      return NO_MSG;
      break;
    case MSG_ERROR:
      // An error occurred when reading the client
      // disconnect from client
      DisconnectClient();
      SET_RETURN_ERROR(MSG_ERROR, MSG_ERROR, "Receive",
		       "Client has Disconnected");
      break;
    default:
      // a message was received
      // save message tag
      m_last_msg_tag = GetMessageTag();

      // Check if registration packet
      int packet_type = GetPacketType();
      if (packet_type == 0) {
	return CheckValidRegistration();
      }
      else {
	return 1;
      }
      break;
    }
  }
  else {
    // no client connected, check if one is ready
    // Note: not checking if we just disconnected it, since the client 
    //  will take awhile to realize and try to reconnect anyway.
    CheckForNewClient();
  }
  
  if (m_client_fd >= 0) {
    // have a client, so just return no message
    return NO_MSG;
  }
  else {
    // no client connected
    SET_RETURN_ERROR2(NO_CLIENT, NO_CLIENT, "Receive",
		     "No client to read from");
  }
}

/***********************************************************/
/*  Name : Send                                            */
/*                                                         */
/*  Description: Send message routine.  This implements the*/
/*    send server protocol - take placemarked packets, add */
/*    to message, and send (unless not completely parsed   */
/*    last message) - After called there will be no packets*/
/*    in queue (will have to add again if there was some   */
/*    error)                                               */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : 1 if successful                         */
/*                 NO_MSG if nothing to send               */
/*                 MSG_ERROR error in send                 */
/*                 NO_CLIENT if no client                  */
/*                 MSG_BUF_BUSY if still using old message */
/*                 MSG_BUF_OVERRUN Msg buffer is too full  */
/*                 PACKET_ERROR - some packet error        */
/***********************************************************/

int Server_Controller::Send() {
  // go through each packet_sending, PutPacketInMsg, then
  //  SendPackets
  if (MessageInUse()) {
    m_num_packets_sending = 0;
    SET_RETURN_ERROR(MSG_BUF_BUSY, MSG_BUF_BUSY, "Send",
		     "Message buffer is busy");
  }

  int i, result = 1;
  for (i=0; i < m_num_packets_sending && result == 1; i++) {
    result = PutPacketInMsg(m_packets_send[i]);
  }

  // reset packets_sending
  m_num_packets_sending = 0;

  if (result == 1) {
    if (m_client_fd >= 0) {
      result = Send_Packets(m_client_fd, m_last_msg_tag);
    }
    else {
      SET_RETURN_ERROR2(NO_CLIENT, NO_CLIENT, "Send",
		       "No client to send to");
    }
  }

  if (result > 0) {
    // some message was sent
    return 1;
  }

  // else some error occurred or there was no message to send.
  //  This is reflected in result, just return it, except if MSG_ERROR
  //  I will also disconnect client
  if (result == MSG_ERROR) {
    DisconnectClient();
    SET_RETURN_ERROR(MSG_ERROR, MSG_ERROR, "Receive",
		     "Client has Disconnected");
  }

  return result;
}


/***********************************************************/
/*  Name : CheckValidRegistration                          */
/*                                                         */
/*  Description: Checks for valid registration, and replies*/
/*    to client w/ the outcome (valid type for packet, or  */
/*    sets tag to error value - 0 for unknown name, and 1  */
/*    for known name, but incorrect size)                  */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : 0 if valid (since want return NO_MSG)   */
/*                 REGISTRATION_ERROR if invalid           */
/***********************************************************/

int Server_Controller::CheckValidRegistration() {
  // registration packet
  Registration_Packet * reg = new Registration_Packet();
  int result = GetPacketFromMsg(reg);
  if (result == 1) {
    // parse valid, get registration info
    const char * reg_class_name = reg->GetRegClassName();
    ushort reg_class_size = reg->GetRegClassSize();
    
    //    printf("wanting to register %s of size %i\n", reg_class_name, reg_class_size);

    // Check if registration matches known packets
    // Set initially to not found
    reg->SetUnknownClassRegTag();
    BOOL valid_reg = FALSE;
    int temp_error_num = -206;

    int i;
    for (i=0; i < m_num_packets_reg; i++) {
      //      printf("Comparing known %s to requested %s\n", m_reg_packet_names[i], reg_class_name);

      if (strcmp(m_reg_packet_names[i], reg_class_name) == 0) {
	// Found the registration
	// Check size
	if (reg_class_size == m_packet_sizes[i]) {
	  // Yes! valid registration
	  reg->SetValidRegTag();
	  reg->SetRegClassType(m_packet_types[i]);
	  valid_reg = TRUE;
	  //	  printf("Valid Registration\n");
	}
	else {
	  // Wrong size
	  //printf("Wrong size server - %i, client %i\n", m_packet_sizes[i], reg_class_size); 
	
	  reg->SetIncorrectSizeRegTag();
	  temp_error_num = -207;
	}
	break;
      }
    }
    // Send registration packet back

    // First must cancel message in case there was anything else
    //  This server can only handle one registration at a time
    Cancel_Message();

    // Add packet and send it
    result = AddPacket(reg);
    if (result != 1) {
      // serious error since all possible errors should be irrelevent
      //  after Cancel_Message, just return registration_error, but keep
      //  the specific error_parameter
      delete reg;
      SET_RETURN_ERROR(REGISTRATION_ERROR, result, 
		       "CheckValidRegistration",
		       "Error in Add_Packet");
    }
      
    result = Send();

    // past this point, no longer need reg
    delete reg;

    if (result != 1) {
      // Again some serious send error, do same
      SET_RETURN_ERROR(REGISTRATION_ERROR, result, 
		       "CheckValidRegistration",
		       "Error in Send");
    }

    if (valid_reg) {
      // was a valid registration, and was sent correctly
      if (m_tread_verbosity > 1) {
	fprintf(stderr, "Server_Controller (CheckValidRegistration) : Registered %s packet\n", reg_class_name);
      }
      return 0;
    }
    else {
      SET_RETURN_ERROR(REGISTRATION_ERROR, temp_error_num, 
		       "CheckValidRegistration",
		       "Received an unknown or wrong size registration packet");
    }
  }
  else {
    // some error in getting registration packet
    delete reg;
    SET_RETURN_ERROR(REGISTRATION_ERROR, result, 
		     "CheckValidRegistration",
		     "Error in GetPacketFromMsg");
  }
}

/***********************************************************/
/*  Name : CheckForNewClient                               */
/*                                                         */
/*  Description: Internal routine for checking for new     */
/*    client.                                              */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : Whether successful or not               */
/*                   1 - client connected                  */
/*                   0 - no client                         */
/*                  SOCKET_ERR - accept error              */
/***********************************************************/

int Server_Controller::CheckForNewClient() {
  // check client non_existent
  if (m_client_fd != -1) {
    if (m_tread_verbosity > 1) 
      fprintf(stderr, "Server_Controller (CheckForNewClient) : Client already exists\n");
    return 1;
  }

  // Otherwise does not exist
  // Set up which socket to watch
  fd_set                readfds;
  FD_ZERO(&readfds);
  FD_SET((uint) m_listening_fd, &readfds);
  
  // Set up poll timeval
  struct timeval        poll;
  poll.tv_usec = 0;
  poll.tv_sec = 0;
  
  // Is there something to read on the socket?
  select(m_listening_fd + 1, &readfds, (fd_set *) 0, 
	 (fd_set *) 0, &poll);
  
  if (FD_ISSET(m_listening_fd, &readfds)) {
    // Reset message buffer information, since any partial message is
    //  now irrelevant
    Cancel_Message();

    // something is coming in on the listening socket
    struct sockaddr_in    clientAddr;
    int                   sockAddrSize = sizeof(struct sockaddr_in);
 
    if ((m_client_fd = accept(m_listening_fd, (struct sockaddr *)&clientAddr, (socklen_t*)&sockAddrSize)) == -1) {
      m_client_fd = -1;
      perror(" accept ");
      SET_RETURN_ERROR(SOCKET_ERR, -205, "CheckForNewClient",
		       "Error in accept call");
    }
    else if (m_tread_verbosity) {
      fprintf(stderr, "Server_Controller (CheckForNewClient) : Accepted client.\n");
      #ifndef RTI_VXWORKS
        struct hostent *ptr;
        long client_addr = ntohl(clientAddr.sin_addr.s_addr);
	ptr = gethostbyaddr((char *)&(client_addr), sizeof(client_addr), AF_INET);

	if (ptr)
	  fprintf(stderr, " from %s on port %hu\n", ptr->h_name, ntohs(clientAddr.sin_port));
      #else
	char host_name_buf[MAXHOSTNAMELEN+1];
	hostGetByAddr(ntohl(clientAddr.sin_addr.s_addr),host_name_buf);
	fprintf(stderr, " from %s\n", host_name_buf);
      #endif
    }
  }
  return m_client_fd >= 0;
}

/***********************************************************/
/*  Name : DisconnectClient                                */
/*                                                         */
/*  Description: Routine for disconnecting from a bad      */
/*    client.                                              */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

void Server_Controller::DisconnectClient() {
  if (m_tread_verbosity) {
    fprintf(stderr, "Server_Controller (DisconnectClient) : Disconnected from client\n");
  
    struct sockaddr_in clientAddr;
    int addrSize = sizeof(clientAddr);
    
    // Get info on the socket (particuarly its address)
    getsockname (m_client_fd, (struct sockaddr *)&clientAddr, (socklen_t*)&addrSize);
    
    // Get name of host
    #ifndef RTI_VXWORKS
      struct hostent *ptr;
      long client_addr = clientAddr.sin_addr.s_addr;
      ptr = gethostbyaddr((char *)&(client_addr), sizeof(client_addr), AF_INET);
    
      if (ptr)
        fprintf(stderr, " which resided on the machine %s\n", ptr->h_name);
    #else
      char host_name_buf[MAXHOSTNAMELEN+1];
      hostGetByAddr(ntohl(clientAddr.sin_addr.s_addr),host_name_buf);
      fprintf(stderr, " which resided on the machine %s\n", host_name_buf);
    #endif
  }

  // Close fd and set back to -1
  close(m_client_fd);
  m_client_fd = -1;
 
  // Reset message buffer information, since any partial message is
  //  now irrelevant
  Cancel_Message();
}
