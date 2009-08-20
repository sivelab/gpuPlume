#include "Client_Controller.h"

/***********************************************************/
/*  Name : Destructor                                      */
/*                                                         */
/*  Description: Destructor                                */
/*                                                         */
/*  Input : none                                           */
/*  Return Value : None                                    */
/***********************************************************/

Client_Controller::~Client_Controller() {
  if (m_server_name != NULL) {
    free(m_server_name);
  }
}

/***********************************************************/
/*  Name : Client_Controller                               */
/*                                                         */
/*  Description: Default Constructor                       */
/*                                                         */
/*  Input : none                                           */
/*  Return Value : None                                    */
/***********************************************************/

Client_Controller::Client_Controller() : Net_Controller() {
  // no socket connections yet, so set file descriptor to -1
  m_server_fd = -1;
  
  m_tread_verbosity = 0;
  
  int i;
  for (i=0; i < MAX_PACKETS_PER_MSG; i++) {
    m_packets_send[i] = NULL;
  }
  m_num_packets_sending = 0;
  m_num_packs_sent = 0;

  m_server_name = NULL;
  m_server_port = -1;
  m_msg_tag = 1;

  // some default receive params - tries alot
  m_num_receive_tries = 3;
  m_max_wait_time = 2000000; // in usecs
}

/***********************************************************/
/*  Name : Register                                        */
/*                                                         */
/*  Description: Register Packet with Server and check if  */
/*    the registration succeeds.                           */
/*                                                         */
/*  Input : Packet ptr                                     */
/*  Return Value : returns Packet type (if successful)     */
/*                 PACKET_ERROR - Invalid packet/packet    */
/*                                  class name             */
/*                 PACKET_DEF_CLASH - if server registry   */
/*                        has different size for this class*/
/*                 REGISTRATION_ERROR - error in           */
/*                   registration process                  */
/*                 SERVER_REG_UNKNOWN - if server does not */
/*                   know this packet class                */
/*                 MSG_BUF_BUSY - buffer is busy           */
/***********************************************************/

int Client_Controller::Register(Packet * new_packet) {
  if (new_packet == NULL) {
    SET_RETURN_ERROR(PACKET_ERROR, PACKET_ERROR, (char *)"Register",
		     (char *)"Packet was Null");
  }
  const char * class_name = new_packet->ClassName();
  //  printf("trying to register packet class %s\n", class_name);
  if (class_name == NULL) {
    SET_RETURN_ERROR(PACKET_ERROR, PACKET_ERROR, (char *)"Register",
		     (char *)"Packet class name was NULL");
  }

  if (MessageInUse() || m_num_packets_sending != 0) {
    SET_RETURN_ERROR(MSG_BUF_BUSY, MSG_BUF_BUSY, (char *)"Register",
		     (char *)"Message buffer is busy");
  }

  Registration_Packet * reg_pack = new Registration_Packet();
  reg_pack->SetRegClassName(class_name);
  reg_pack->SetRegClassSize(new_packet->GetSize());

  int result = AddPacket(reg_pack);
  if (result != 1) {
    // should not occur, but in case (can always read exact error
    //  by printing last error message)
    delete reg_pack;
    return REGISTRATION_ERROR;
  }

  result = Send();

  // no longer need reg_pack
  delete reg_pack;

  if (result != 1) {
    return REGISTRATION_ERROR;
  }

  result = Receive(TRUE);
  if (result != 1) {
    return REGISTRATION_ERROR;
  }

  int check_packet_type = GetPacketType();
  if (check_packet_type != 0) {
    // major error, did not send back the registration packet
    Cancel_Message();
    SET_RETURN_ERROR(REGISTRATION_ERROR, -305, (char *)"Register",
		     (char *)"Returned packet from Registration was wrong type");
  }

  Registration_Packet * returned_reg = new Registration_Packet();
  result = GetPacketFromMsg(returned_reg);
  
  // We are now done with the message, make sure that is all that
  //  came with it
  int msg_end = GetPacketType();
  if (msg_end != -1) {
    // something else? Not good, cancel
    Cancel_Message();
    delete returned_reg;
    SET_RETURN_ERROR(REGISTRATION_ERROR, -305, (char *)"Register",
		     (char *)"More than one packet was sent back w/ Registration");
  }

  if (result == 1) {
    // parse valid, get registration info
    const char * ret_class_name = returned_reg->GetRegClassName();
    ushort ret_class_size = returned_reg->GetRegClassSize();
    
    // ensure sent back the correct registration packet
    if (!(strcmp(ret_class_name, class_name) == 0) || 
	!(ret_class_size == new_packet->GetSize())) {
      // wrong registration packet sent back
      delete returned_reg;
      SET_RETURN_ERROR(REGISTRATION_ERROR, -306, (char *)"Register",
		       (char *)"Wrong registration class/size sent back");
    }
    
    // returned valid registration packet, check if server new packet
    TRRegTag reg_result = returned_reg->GetRegTag();
    uint class_type = returned_reg->GetRegClassType();

    // no longer need returned_reg
    delete returned_reg;

    switch (reg_result) {
    case VALID_REG:
      if (m_tread_verbosity > 1) {
	fprintf(stderr, "Client_Controller (Register) : Registered %s packet\n", class_name);
      }
      // Tell packet it is registered
      SetPacketType(new_packet,class_type);
      SetPacketRegistered(new_packet);
      return class_type;
      break;
    case INVALID_REG_SIZE:
      SET_RETURN_ERROR(PACKET_DEF_CLASH, PACKET_DEF_CLASH, (char *)"Register",
		       (char *)"Server has a different size for this packet class");
      break;
    case UNKNOWN_CLASS:
      SET_RETURN_ERROR(SERVER_REG_UNKNOWN, SERVER_REG_UNKNOWN, (char *)"Register",
		       (char *)"Server does not recognize this packet class name");
      break;
    default:
      // sent back some unknown tag
      SET_RETURN_ERROR(REGISTRATION_ERROR, -307, (char *)"Register",
		       (char *)"Server sent back an unknown Registration Tag Value");
    } 
  }
  else {
    // else just return a registration_error, the call would have set an
    //  appropriate error message that can be read
    return REGISTRATION_ERROR;
  }
}

/***********************************************************/
/*  Name : ConnectToServer                                 */
/*                                                         */
/*  Description: Connect to TCP server on the specified    */
/*                machine at given port                    */
/*                                                         */
/*  Input : server hostname, port number                   */
/*  Return Value : if successful (1)                       */
/*                 if socket failure (SOCKET_ERR)          */
/*                 if no server to connect to (NO_SERVER)  */
/***********************************************************/

int Client_Controller::ConnectToServer(char * hostname, short port) {
  if (m_server_fd != -1) 
    return 1;

#ifdef WIN32
  WORD wVersionRequested = MAKEWORD( 2, 0 );
  WSADATA wsaData;
  int err = WSAStartup( wVersionRequested, &wsaData );
  if (err != 0) {
    // Tell the user that we couldn't find a usable 
    // WinSock DLL.                                  
    SET_RETURN_ERROR(SOCKET_ERR, SOCKET_ERR, (char *)"ConnectToServer",
		     (char *)"Couldn't find a useable WinSock DLL.");
  }
#endif

  // Create socket
  if ((m_server_fd = socket(AF_INET, SOCK_STREAM, 0)) == -1) {
    if (m_tread_verbosity) {
      perror("Client_Controller: ConnectToServer (socket)");
    } 
    SET_RETURN_ERROR(SOCKET_ERR, -301, (char *)"ConnectToServer", 
		     (char *)"Socket call failed");
  }

  // Initialize and set server socket address information
  // Using TCP/IP sock streams
  struct sockaddr_in server_addr;

  // Initialize it to all zeros
  memset((char *)&server_addr, 0, sizeof(server_addr));

  server_addr.sin_family = AF_INET;
  server_addr.sin_port = htons((unsigned short) port);
  server_addr.sin_addr.s_addr = GetAddrByName(hostname);

  // Set socket options to TCP nodelay (i.e. send packets as soon as possible)
  int on = 1;
  
  if (setsockopt(m_server_fd, IPPROTO_TCP, TCP_NODELAY, (char *)&on, sizeof(int)) == -1) {
    close(m_server_fd);
    m_server_fd = -1;
    if (m_tread_verbosity) {
      perror("Client_Controller: ConnectToServer (setsockopt)");
    }
    SET_RETURN_ERROR(SOCKET_ERR, -302, (char *)"ConnectToServer",
		     (char *)"SetSockOpt call failed");
  }
  
  // Connect the socket to the server address
  if (connect(m_server_fd, (struct sockaddr *)&server_addr, sizeof(server_addr)) == -1) {
    close(m_server_fd);
    m_server_fd = -1;
    if (m_tread_verbosity) {
      perror("Client_Controller: ConnectToServer (connect)");
    } 
    SET_RETURN_ERROR(NO_SERVER, -303, (char *)"ConnectToServer",
		     (char *)"Connect call failed - nothing to connect to");
  }
  
  if (m_tread_verbosity) 
    fprintf(stderr, "Client_Controller (ConnectToServer) : Connected to server on %s, port %d. \n", 
	    hostname, port);
  
  if (m_server_name != NULL)
    free(m_server_name);
  m_server_name = strdup(hostname);
  m_server_port = port;

  return 1;
}

/***********************************************************/
/*  Name : Disconnect                                      */
/*                                                         */
/*  Description: Disconnect client (close socket and remove*/
/*                message)                                 */
/*                                                         */
/*  Note - This also removes knowledge of the server, i.e. */
/*    the user has asked to disconnect from server so it   */
/*    will not try to reconnect.                           */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

void Client_Controller::Disconnect() {
  // First call private function for closing the server
  //  connection and cleaning up message
  CloseServerConnection();

  // Also remove knowledge of server
  if (m_server_name != NULL) 
    free(m_server_name);
  m_server_name = NULL;
  m_server_port = -1;
}

/***********************************************************/
/*  Name : AddPacket                                       */
/*                                                         */
/*  Description: Add packet routine.  This implements the  */
/*    client protocol - packets are added right now, thus  */
/*    the message cannot be busy (i.e. must have been      */
/*    totally read/parsed)                                 */
/*                                                         */
/*  Notes: 1) Consequence of this is the packet''s data can*/
/*     NOT be altered in between AddPacket and Send.       */
/*                                                         */
/*  Input : Packet Ptr                                     */
/*  Return Value : 1 - if successful                       */
/*                 PACKET_ERROR - if packet is invalid,    */
/*                       null, unregistered, fills wrong)  */
/*                 OUT_OF_PACKET_SLOTS - if too many       */
/*                       packets in message                */
/*                 MSG_BUF_BUSY if have not finished       */
/*                    parsing last msg                     */
/*                 Message buffer overrun (MSG_BUF_OVERRUN)*/
/***********************************************************/

int Client_Controller::AddPacket(Packet * new_packet) {
  if (MessageInUse()) {
    SET_RETURN_ERROR(MSG_BUF_BUSY, MSG_BUF_BUSY, (char *)"AddPacket",
		     (char *)"Message buffer is busy");
  }

  if (m_num_packets_sending >= MAX_PACKETS_PER_MSG) {
    SET_RETURN_ERROR(OUT_OF_PACKET_SLOTS, OUT_OF_PACKET_SLOTS, (char *)"AddPacket",
		     (char *)"Too many packets being sent - increase MAX_PACKETS_PER_MSG");
  }

  if (new_packet != NULL && new_packet->IsRegistered()) {
    m_packets_send[m_num_packets_sending] = new_packet;
    int result = PutPacketInMsg(m_packets_send[m_num_packets_sending]);
    m_num_packets_sending++;
    if (result < 0) {
      // since all the errors from this set the error values
      //  correctly, I will just return the error
      return result;
    } 
  }
  else {
    SET_RETURN_ERROR(PACKET_ERROR, -103, (char *)"AddPacket",
		     (char *)"Unregistered or Null packet");
  }

  return 1;
}


/* RH - Changes!!!! - Add protocol for checking that all responses come???? */

/***********************************************************/
/*  Name : Receive                                         */
/*                                                         */
/*  Description: This procedure checks for incoming        */
/*    messages from the server and implements the client   */
/*    protocol (like can not receive until old message has */
/*    been fully parsed), and tries multiple times to get  */
/*    a message (if wait is TRUE).                         */
/*                                                         */
/*  Notes(1) : Do not call method with wait = TRUE unless  */
/*          expecting a message back (for example you sent */
/*          a request to server for info).  Otherwise the  */
/*          routine will not return for quite awhile and   */
/*          be pointless since there is nothing to receive */
/*       (2) : The message tag returned must match the last*/
/*          one sent, because the only messages the client */
/*          receives are from requests (need to ensure the */
/*          the response matches the request)              */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : 1 if received some message              */
/*                 NO_MSG if no message and !wait          */
/*                 MSG_BUF_BUSY if have not finished       */
/*                    parsing last msg                     */
/*                 MSG_ERROR if error in read              */
/*                 NO_SERVER if no server exists           */
/*                 RECEIVE_TIMEOUT if no message and wait  */
/*                 INCORRECT_RESPONSE_TAG if incorrect tag */
/***********************************************************/

int Client_Controller::Receive(BOOL wait) {
  if (MessageInUse()) {
    SET_RETURN_ERROR(MSG_BUF_BUSY, MSG_BUF_BUSY, (char *)"Receive",
		     (char *)"Message Buffer Busy");
  }

  if (m_server_fd < 0) {
    // no connection, try to reconnect
    if (!TryReconnectToServer())
      SET_RETURN_ERROR2(NO_SERVER, NO_SERVER, (char *)"Receive",
			(char *)"Not connected to a server");
  }
      
  // receive message since client exists
  int rec_result;
  int wait_time;
  int num_tries = 0;
  if (wait) {
    // want to wait, set params appropriately
    wait_time = m_max_wait_time;
  }
  else {
    // do not want to wait, set on last try and no wait
    num_tries = m_num_receive_tries - 1;
    wait_time = POLL;
  }
  
  while (num_tries < m_num_receive_tries) {
    rec_result = Receive_Message(m_server_fd, wait_time);
    switch(rec_result) {
    case MSG_ERROR:
      // An error occurred when reading the server
      // disconnect from server
      CloseServerConnection();
      SET_RETURN_ERROR(MSG_ERROR, MSG_ERROR, (char *)"Receive",
		       (char *)"Server Connection Lost");
      break;
    case NO_MSG:
      // try again
      break;
    default:
      // a message was received
      // ensure it is the correct response
      ushort response_tag = GetMessageTag();
      if (response_tag != m_msg_tag) {
	// incorrect response, cancel message since it is invalid
	Cancel_Message();
	SET_RETURN_ERROR(INCORRECT_RESPONSE_TAG, INCORRECT_RESPONSE_TAG, (char *)"Receive",
			 (char *)"Tag from response message did not match outgoing tag - incorrect response");
      }
      return 1;
      break;
    }
    num_tries++;
  }
  
  // if got here that means either we timed out (if waiting), or
  //  just no message if not waiting
  if (wait) {
    // have a client, so just return no message
    SET_RETURN_ERROR(RECEIVE_TIMEOUT, RECEIVE_TIMEOUT, (char *)"Receive",
		     (char *)"Did not receive a message in the time allotted");
  }
  else {
    // was not waiting, so just return no message
    return NO_MSG;
  }
}

/***********************************************************/
/*  Name : Send                                            */
/*                                                         */
/*  Description: Send message routine.  This implements the*/
/*    send client protocol - take placemarked packets, add */
/*    to message, and send (unless not completely parsed   */
/*    last message) - After this call the packet queue     */
/*    will be empty, so if an error occurred AddPacket     */
/*    must be called again to send that packet.            */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : 1 if successful                         */
/*                 NO_MSG if nothing to send               */
/*                 MSG_ERROR error in send                 */
/*                 NO_SERVER if no server                  */
/*                 MSG_BUF_BUSY if still using old message */
/***********************************************************/

int Client_Controller::Send() {
  // reset packets sending, and save sent
  int temp_packs_sent = m_num_packets_sending;

  m_num_packs_sent = 0; // assumption none sent creates less code
  m_num_packets_sending = 0;

  //  Check if message in use
  if (MessageInUse()) {
    SET_RETURN_ERROR(MSG_BUF_BUSY, MSG_BUF_BUSY, (char *)"Send",
		     (char *)"Message buffer is busy");
  }

  if (m_server_fd < 0) {
    // no connection, try to reconnect
    if (!TryReconnectToServer()) {
      // need to cancel message so buffer does not fill up w/ unsent packets
      Cancel_Message();
      SET_RETURN_ERROR2(NO_SERVER, NO_SERVER, (char *)"Send",
			(char *)"Not connected to a server");
    }
  }
      
  // increase message tag and Send Packets
  m_msg_tag++;
  int result = Send_Packets(m_server_fd, m_msg_tag);

  if (result > 0) {
    // no problem, send back 1 and set sent
    m_num_packs_sent = temp_packs_sent;
    return 1;
  }

  // else some problem or no message, send back result
  //  plus if MSG_ERROR - close connection
  if (result == MSG_ERROR) {
    // Close connection
    CloseServerConnection();
    SET_RETURN_ERROR(MSG_ERROR, MSG_ERROR, (char *)"Send",
		     (char *)"Server Connection Lost");
  }
  return result;
}

/***********************************************************/
/*  Name : SendSinglePacket                                */
/*                                                         */
/*  Description: Send a single packet to server.  The      */
/*   packet queue must be empty and message buffer free    */
/*   for this to succeed.                                  */
/*                                                         */
/*  Input : ptr to packet to send                          */
/*  Return Value : 1 if successful                         */
/*                 0 if fail (use PrintLastError to get    */
/*                  the reason for failure)                */
/***********************************************************/

int Client_Controller::SendSinglePacket(Packet * send_pkt) {
  if ((m_num_packets_sending == 0) && 
      (AddPacket(send_pkt) == 1) &&
      (Send() == 1)) {
    return 1;
  }
  return 0;
}

/***********************************************************/
/*  Name : ReceiveSinglePacket                             */
/*                                                         */
/*  Description: Receive a single packet from server (used */
/*    in conjuction with SendSinglePacket w/ the packet    */
/*    being a response from server)                        */
/*                                                         */
/*  Notes(1) : This function will expect to receive a      */
/*          packet, so do not call if a request packet was */
/*          not just sent.                                 */
/*       (2) : Any message that is received will be        */
/*          cancelled if it is not correct.                */
/*                                                         */
/*  Input : Ptr to packet which is being received          */
/*  Return Value : 1 if successful                         */
/*                 0 if fail (use PrintLastError to get    */
/*                  the reason for failure)                */
/***********************************************************/

int Client_Controller::ReceiveSinglePacket(Packet * rec_pkt) {
  if (Receive(TRUE) == 1) {
    if (rec_pkt != NULL) {
      // check type
      if (GetPacketType() == (int) rec_pkt->GetType()) {
	// received correct packet, get it
	if (GetPacketFromMsg(rec_pkt) == 1) {
	  // ensure that was the only packet in message
	  if (GetPacketType() != -1) {
	    // since I did get a valid response, I will allow this
	    //  but cancel the rest of the message
	    if (m_tread_verbosity) {
	      fprintf(stderr, "Error (ReceiveSinglePacket) : Received more than one packet\n");
	    }
	    Cancel_Message();
	  }
	  return 1;
	}
      }
      else {
	m_error_val = INCORRECT_PACKET_RECEIVED;
	SetErrorMessage((char *)"ReceiveSinglePacket",(char *)"Received incorrect packet type");
	if (m_tread_verbosity) {
	  fprintf(stderr, "%s\n", m_error_msg);
	}
      }
    }
  }
  // if got here then some error occurred (error val is set correctly
  //  for PrintLastError).  So, cancel message and return 0.
  Cancel_Message();
  return 0;
}

/***********************************************************/
/*  Name : TryReconnectToServer                            */
/*                                                         */
/*  Description: Internal routine for trying to reconnect  */
/*    to server                                            */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : Whether successful or not               */
/*                   1 - server connected                  */
/*                   0 - no server                         */
/***********************************************************/

int Client_Controller::TryReconnectToServer() {
  // check that server is non_existent
  if (m_server_fd != -1) {
    if (m_tread_verbosity > 1) 
      fprintf(stderr, "Error (TryReconnectToServer) : Already connected to Server\n");
    return 1;
  }

  if (m_server_name == NULL || m_server_port < 0) {
    // nothing to connect to - either have not called ConnectToServer or 
    //  called Disconnect since
    SET_RETURN_ERROR(0, -304, (char *)"TryReconnectToServer",
		     (char *)"No server/port specified to try and reconnect to");
  }

  if (m_tread_verbosity > 1)
    fprintf(stderr, "Client_Controller (TryReconnectToServer) : Attempting to reconnect to server %s on %d\n",
	    m_server_name, m_server_port);

  if (ConnectToServer(m_server_name, m_server_port) == 1) {
    // succeeded
    if (m_tread_verbosity)
      fprintf(stderr, "Client_Controller (TryReconnectToServer) : Connection succeeded\n");
    return 1;
  }
  
  // else some error or no server to connect to, just return 0
  return 0;
}

/***********************************************************/
/*  Name : CloserServerConnection                          */
/*                                                         */
/*  Description: Protected routine for closing socket and  */
/*                cleaning up message)                     */ 
/*                                                         */
/*  Note - This will NOT remove knowledge of the server,   */
/*    so that we can try reconnecting to it (i.e. this was */
/*    called because of a read failure, so maybe the server*/
/*    was rebooted - still try to reconnect)               */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

void Client_Controller::CloseServerConnection() {
  if (m_server_fd >= 0) {
    // if server socket exists
    close(m_server_fd);
    m_server_fd = -1;
  }
  
  if (m_tread_verbosity) 
    fprintf(stderr, "Client_Controller (CloseServerConnection) : Client disconnected\n");

  // remove any pending message
  Cancel_Message();
  m_num_packets_sending = 0;
  m_num_packs_sent = 0;

  #ifdef WIN32
  WSACleanup();
  #endif
}

/***********************************************************/
/*  Name : GetAddrByName                                   */
/*                                                         */
/*  Description: Internal routine for trying to reconnect  */
/*    to server                                            */
/*                                                         */
/*  Input : Name of computer                               */
/*  Return Value : Address of computer (if successful)     */
/*                 0 - if can not find computer address    */
/***********************************************************/

long Client_Controller::GetAddrByName(char *name) {
#ifndef RTI_VXWORKS
  struct hostent * ptr;
  ptr = gethostbyname(name);
  if (ptr) {
    // essentially *h_addr_list will look like 192.128.55.2 (char array)
    //  except without the periods.  Treat this array of characters
    //  as a long and return that value.
    return (*(long*)*(ptr->h_addr_list));
  }
  else return 0;
#else
  return hostGetByName(name);
#endif
}
