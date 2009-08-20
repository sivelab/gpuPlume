#include "Net_Controller.h"

/***********************************************************/
/*  Name : Destructor                                      */
/*                                                         */
/*  Description: Destructor                                */
/*                                                         */
/*  Input : none                                           */
/*  Return Value : None                                    */
/***********************************************************/

Net_Controller::~Net_Controller() {
  delete [] m_msg_buffer;
  delete [] m_error_msg;
}

/***********************************************************/
/*  Name : Default Constructor                             */
/*                                                         */
/*  Description: Default Constructor                       */
/*                                                         */
/*  Input : none                                           */
/*  Return Value : None                                    */
/***********************************************************/

Net_Controller::Net_Controller() {
  // Initialize variables
  m_error_val = 0;
  m_error_msg = new char[ERR_MSG_SIZE];

  // Initialize message buffer
  m_msg_buffer = new char[MESSAGE_BUFFER_SIZE];
  m_msg_pointer = 0;
  m_new_message = FALSE;

  // A message size of -1 means we do not yet know the size (just 
  //  started reading the message)
  m_complete_message_size = -1;
}

/***********************************************************/
/*  Name : PrintLastError                                  */
/*                                                         */
/*  Description: Prints last error message and returns last*/
/*                error value                              */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : Last error value                        */
/***********************************************************/

int Net_Controller::PrintLastError() {
  fprintf(stderr, "%s\n", m_error_msg);
  return m_error_val;
}

/***********************************************************/
/*  Name : SetErrorMessage                                 */
/*                                                         */
/*  Description: Sets error message                        */
/*                                                         */
/*  Input : func_name - name of function where error       */
/*                       occurred                          */
/*          err_msg - error message                        */
/*  Return Value : if successful (1)                       */
/*                 if fail (-1) - message too big          */
/***********************************************************/

int Net_Controller::SetErrorMessage(char * func_name, char * err_msg) {
  // error looks like "Error (func_name) : err_msg"
  if (12 + strlen(func_name) + strlen(err_msg) > ERR_MSG_SIZE) {
    return -1;
  }

  strcpy(m_error_msg,"Error (");
  strcat(m_error_msg,func_name);
  strcat(m_error_msg,") : ");
  strcat(m_error_msg,err_msg);
  return 1;
}

/***********************************************************/
/*  Name : Cancel_Message                                  */
/*                                                         */
/*  Description: Cancel message, whether it is something   */
/*   being prepared or received.                           */
/*                                                         */
/*  Note: 1) This is useful if for some reason cannot      */
/*          parse the message correctly, and just want to  */
/*          quite trying (i.e. trash it).                  */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

void Net_Controller::Cancel_Message() {
  m_new_message = FALSE;
  m_msg_pointer = 0;
  m_complete_message_size = -1;
}

/***********************************************************/
/*  Name : GetMessageTag                                   */
/*                                                         */
/*  Description: GetMessageTag from buffer (byte 2&3)      */
/*                                                         */
/*  Note: 1) This just returns the short value in byte 2&3 */
/*          So if there is not a valid message in the buffer*/
/*          then the result will not make sense (all calls */
/*          should happen when valid anyway)               */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

ushort Net_Controller::GetMessageTag() {
  short * val_ptr = (short *)(m_msg_buffer+2);
  return ntohs(*val_ptr);
}

/***********************************************************/
/*  Name : GetPacketType                                   */
/*                                                         */
/*  Description: Get packet type from message              */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : Packet type of next packet in message   */
/*                  or -1 if seen full message             */
/***********************************************************/

int Net_Controller::GetPacketType() {
  if (!m_new_message) {
    return -1;
  }

  if (m_msg_pointer == 0) {
    // first packet attempted to check
    m_msg_pointer = 4;
  }

  if (m_msg_pointer >= m_complete_message_size) {
    // seen full message, set new message to false since completed
    m_new_message = FALSE;
    m_msg_pointer = 0;
    m_complete_message_size = -1;
    return -1;
  }
  
  // else get type, which is just a character
  int packet_type = (int)m_msg_buffer[m_msg_pointer];
  return packet_type;
}


/************************************************************/
/*  Name : GetPacketFromMsg                                 */
/*                                                          */
/*  Description: Parse one packet of the message.           */
/*                                                          */
/*  Notes: The supplied Packet is checked to verify that    */
/*   it is of the correct type (should have called          */
/*   GetPacketType)                                         */
/*                                                          */
/*  Input : None                                            */
/*  Return Value : if successful (1)                        */
/*                 if no message to parse (0)               */
/*                 if packet error (PACKET_ERROR)           */
/*                 if out of message space (MSG_BUF_OVERRUN)*/
/************************************************************/

int Net_Controller::GetPacketFromMsg(Packet * return_packet) {
  if (!m_new_message) {
    // we have not received a new message so nothing to parse
    return 0;
  }

  uint packet_type = (int)m_msg_buffer[m_msg_pointer];
  if (packet_type != return_packet->GetType()) {
    SET_RETURN_ERROR(PACKET_ERROR, -101, (char *)"GetPacketFromMsg",
		     (char *)"Input packet/message type mismatch");
  }

  // else correct type, get the packet
  m_msg_pointer++;

  int packet_size = return_packet->GetSize();
  if (packet_type != 0) {
    // Only can check size for non_registration packets
    if (packet_size + m_msg_pointer > m_complete_message_size) {
      // not enough data
      // reset new message to false
      m_new_message = FALSE;
      SET_RETURN_ERROR(MSG_BUF_OVERRUN, MSG_BUF_OVERRUN, (char *)"GetPacketFromMsg",
		       (char *)"Message buffer overrun");
    }
  }

  int result = return_packet->ReadPacket(m_msg_buffer + m_msg_pointer);
  if (result == 1) {
    // need to get packet_size again, in case registration packet
    packet_size = return_packet->GetSize();

    // set that a new packet has been received
    m_msg_pointer += packet_size;
    if (m_msg_pointer > m_complete_message_size) {
      // error, read past end of message - should only happen with bad
      //  registration packets
      SET_RETURN_ERROR(MSG_BUF_OVERRUN, MSG_BUF_OVERRUN, (char *)"GetPacketFromMsg",
		       (char *)"Message buffer overrun");
    }
    return 1;
  }
  
  // else error in read (change later - what can go wrong???)
  SET_RETURN_ERROR(PACKET_ERROR, -102, (char *)"GetPacketFromMsg",
		   (char *)"Packet Read Error");
}

/***********************************************************/
/*  Name : PutPacketInMsg                                  */
/*                                                         */
/*  Description: Adds packet into message                  */
/*                                                         */
/*  Notes: 1) Checks that packet is registered             */
/*         2) Assumes that packet''s fill routine is       */
/*             implemented and checks that the size it     */
/*             fills matches GetSize()                     */
/*         3) This will override the previous message      */
/*             received, which is why it is a protected    */
/*             function (the child controller cannot call  */
/*             this method unless m_new_message = FALSE)   */
/*                                                         */
/*  Input : Packet to add                                  */
/*  Return Value : Success (1)                             */
/*                 Old Message still in use (MSG_BUF_BUSY) */
/*                 Message buffer overrun (MSG_BUF_OVERRUN)*/
/*                 Invalid/Unregistered/Fill Error on      */
/*                  Packet (PACKET_ERROR)                  */
/***********************************************************/

int Net_Controller::PutPacketInMsg(Packet * packet_add) {
  if (m_new_message) {
    SET_RETURN_ERROR(MSG_BUF_BUSY, MSG_BUF_BUSY, (char *)"PutPacketInMsg", 
		     (char *)"Message Buffer Busy");
  }

  // Fix message pointer if this is the first packet put in message
  if (m_msg_pointer < 4) {
    m_msg_pointer = 4;
  }

  if (packet_add != NULL && packet_add->IsRegistered()) {
    // These hold the size of the current packet, and how many bytes
    //  it filled the message buffer with
    int ind_packet_size, fill_size;

    ind_packet_size = packet_add->GetSize();
    if (m_msg_pointer + ind_packet_size + 1 > MESSAGE_BUFFER_SIZE) {
      // Error we will overrun the buffer
      SET_RETURN_ERROR(MSG_BUF_OVERRUN, MSG_BUF_OVERRUN, (char *)"PutPacketInMsg", 
		       (char *)"Not enough space in message buffer");
    }

    // Otherwise add into message
    m_msg_buffer[m_msg_pointer++] = (char)packet_add->GetType();
    fill_size = packet_add->Fill(m_msg_buffer + m_msg_pointer);
	
    if (fill_size != ind_packet_size) {
      // Error in fill procedure.  It did not fill the
      //  message buffer with the correct amount of data
      SET_RETURN_ERROR(PACKET_ERROR, -104, (char *)"PutPacketInMsg", 
		       (char *)"Packet's fill command error - incorrect fill bytes");
    }
	
    // otherwise it worked correctly, continue on
    m_msg_pointer += ind_packet_size;
  }
  else {
    // invalid
    SET_RETURN_ERROR(PACKET_ERROR, -103, (char *)"PutPacketInMsg", 
		     (char *)"Unregistered Packet or Null Packet");
  }
  return 1;
}

/***********************************************************/
/*  Name : Send_Packets                                    */
/*                                                         */
/*  Description: Send all packets that have been added     */
/*                                                         */
/*  Note:  If an error is returned the file descriptor     */
/*      should be closed so that the other side will be    */
/*      notified there is a problem (it will see the fd    */
/*      closed when it tries to read packets).             */
/*                                                         */
/*  Input : file descriptor for socket to use and tag      */
/*           to attach unto message                        */
/*  Return Value : if successful it returns the # of bytes */
/*                      sent                               */
/*                 if nothing to send (NO_MSG)             */
/*                 or if failed (MSG_ERROR)                */
/*                 if message buf busy (MSG_BUF_BUSY)      */
/***********************************************************/

int Net_Controller::Send_Packets(int fd, ushort tag) {
  // When this routine is called all the packets have already
  //  been added to the message.  Only need to fix the size
  //  (client or server classes set the tag).

  // Ensure that the m_new_message is false, otherwise have not
  //  finished with previous message
  if (m_new_message) {
    SET_RETURN_ERROR(MSG_BUF_BUSY, MSG_BUF_BUSY, (char *)"Send_Packets", (char *)"Message buffer busy");
  }

  ushort message_size = (ushort) m_msg_pointer;
  
  // Fix message size and send
  if (m_msg_pointer > 4) {
    // Need to put in network byte order
    ushort nbo_msg_size = htons(message_size);
    char *ptr_val = (char *)&nbo_msg_size;
    m_msg_buffer[0] = ptr_val[0];
    m_msg_buffer[1] = ptr_val[1];

    ushort nbo_tag = htons(tag);
    ptr_val = (char *)&nbo_tag;
    m_msg_buffer[2] = ptr_val[0];
    m_msg_buffer[3] = ptr_val[1];

    // temporary, print out what sending
    //        for (int i=0; i < message_size; i++) {
    //  if (m_msg_buffer[i] > 30) {
    //    printf("%c",m_msg_buffer[i]);
    //  }
    //  else {
    //    printf(" %i ", (int)m_msg_buffer[i]);
    //  }
    //  }
    //  printf("\n");


    // send message
    int send_result = send(fd, m_msg_buffer, message_size, 0);
    
    // reset message pointer
    m_msg_pointer = 0;

    if (send_result != message_size) {
      // error, not all of the message was sent
      // Note: I am not handling sending all of large messages.  Our
      //  largest packet is still quite small so I do not expect this
      //  to be a problem.  If the packets begin to be quite large
      //  then multiple sends will have to be handled.
      SET_RETURN_ERROR(MSG_ERROR, MSG_ERROR, (char *)"Send_Packets", (char *)"Send error - incorrect bytes sent");
    }
    // otherwise just return the message_size
    return message_size;
  }
  
  // else no packets were ready, so sent 0 bytes
  return NO_MSG;
}

/***********************************************************/
/*  Name : Receive_Message                                 */
/*                                                         */
/*  Description: Receive message on the given socket,      */
/*   waiting no longer than waittime_usec microseconds.    */
/*                                                         */
/*  Note: 1) If MSG_ERROR, the file descriptor should      */
/*          be closed since we are incorrectly reading it. */
/*        2) This function will override the message       */
/*          so if have added packets, should call send     */
/*          before receive.                                */
/*                                                         */
/*  Input : int fd - File descriptor for socket to receive */
/*                    packets from                         */
/*          int waittime_usec - number of usecs to wait    */
/*                               for packets               */
/*  Return Value : size of message                         */
/*                 NO_MSG if no message was received       */
/*                 MSG_ERROR if read error                 */
/*                 if still on old message (MSG_BUF_BUSY)  */
/***********************************************************/

int Net_Controller::Receive_Message(int fd, int waittime_usec) {
  // Ensure that the m_new_message is false, otherwise have not
  //  finished with previous message
  if (m_new_message) {
    SET_RETURN_ERROR(MSG_BUF_BUSY, MSG_BUF_BUSY, 
		     (char *)"Receive_Message", (char *)"Message Buffer Busy");
  }

  // Registration Packet will count as No Packet since
  //  for external purposes should not know about registration
  //  packets.  

  int readresult;
  int read_length;
  
  m_msg_pointer = 0;
  m_complete_message_size = -1;

  while (1) {
    // check if there is something to read on file descriptor
    // If not return no packet
    if (!ReadNonBlocking(fd, waittime_usec)) 
      return NO_MSG;
    
    if (m_complete_message_size < 0) {
      // if we have not read what size the message is then we
      //  are only wanting 2 bytes of information (assuming
      //  the size is a short value).

      // Have to use m_msg_pointer for the rare case that we have
      //  already read one byte of the information
      read_length = 2 - m_msg_pointer;
    }
    else {
      // we have the complete size, read the rest of the message
      read_length = m_complete_message_size - m_msg_pointer;
    }

    readresult = recv(fd, m_msg_buffer + m_msg_pointer,
		      read_length, 0);
  
    if (readresult <= 0) {
      // Either read error, or end of file (i.e. the socket has been
      //  closed on the other side)

      // Reset message pointer and message size
      m_complete_message_size = -1;
      m_msg_pointer = 0;
  
      SET_RETURN_ERROR(MSG_ERROR, MSG_ERROR, (char *)"Receive_Message", 
		       (char *)"Read from file descriptor failed");
    }
    else {
      // received some bytes

      // update message_pointer
      m_msg_pointer += readresult;

      if (m_complete_message_size < 0 && m_msg_pointer >= 2) {
	// If have not received the size previously, but just have
	//  -> get the size (again assuming 2 bytes - ushort)
	m_complete_message_size = ntohs(*(ushort*)(m_msg_buffer));
	if (m_complete_message_size > MESSAGE_BUFFER_SIZE) {
	  // The message we are trying to read is too big, report error
	  SET_RETURN_ERROR(MSG_ERROR, MSG_ERROR, (char *)"Receive_Message", 
			   (char *)"Message sent is too big");
	}	  
      }
      if (m_msg_pointer == m_complete_message_size) {
	// We have received the whole message
	m_msg_pointer = 0;
	m_new_message = TRUE;

	// temporary, print out what received
	//	for (int i=0; i < m_complete_message_size; i++) {
	//  if (m_msg_buffer[i] > 30) {
	//    printf("%c",m_msg_buffer[i]);
	//  }
	//  else {
	//    printf(" %i ", (int)m_msg_buffer[i]);
	//  }
	//}
	//printf("\n");

	// Message has the type of packets stored as a uint right
	//  after size, so return this.
	return m_complete_message_size;
      }
    }
  }
}

/***********************************************************/
/*  Name : ReadNonBlocking                                 */
/*                                                         */
/*  Description: Check socket for read status w/o blocking */
/*                                                         */
/*  Input : int fd - File descriptor for socket to read    */
/*          int microsecs - number of usecs to wait        */
/*  Return Value :  Whether fd is ready to read            */
/***********************************************************/

int Net_Controller::ReadNonBlocking(int fd, int microsecs) {
  // if file descriptor is invalid, obviously not ready to read
  if (fd < 0) return 0;

  // Set up a fd_set for our file descriptor
  fd_set readfds;

  FD_ZERO(&readfds);
  FD_SET((uint)fd, &readfds);

  // Set up how long to wait for select to check read status
  //  (i.e. to make the read nonblocking)
  struct timeval timeout;

  // Note I am assuming less than one second here, but the
  //  server never blocks, and the client should be running
  //  at around 30 Hz.....
  timeout.tv_usec = microsecs;
  timeout.tv_sec = 0;

  // Check fd only for its read status
  select(fd+1, &readfds, (fd_set *)0, (fd_set *)0, 
	 (struct timeval *)&timeout);

  // return if it is ready for read
  return FD_ISSET(fd, &readfds);
}

