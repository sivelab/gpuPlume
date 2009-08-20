#ifndef NET_CONTROLLER_H
#define NET_CONTROLLER_H

#include "Registration_Packet.h"
#include <sys/types.h>

#ifndef WIN32
#include <sys/socket.h>
#include <netinet/tcp.h>
#include <netdb.h>
#else
#define close(x) closesocket(x)
#define socklen_t int
#endif

#ifdef RTI_VXWORKS
#include <sockLib.h>
#include <hostLib.h>
#define socklen_t int
#endif

// Also defining a specific message_type NO_MSG to be zero
//  and MSG_ERROR to be -1
#define NO_MSG 0
#define MSG_ERROR -1         // (read/send) errors
#define MSG_BUF_BUSY -2
#define MSG_BUF_OVERRUN -3
#define PACKET_ERROR -4      // any error on the packet

#define RETURN_ERROR(ret_val) { \
  if (m_tread_verbosity) \
    fprintf(stderr, "%s\n", m_error_msg); \
  return ret_val; \
}

#define SET_RETURN_ERROR(ret_val, error_val, func_name, err_msg) { \
  m_error_val = error_val; \
  SetErrorMessage(func_name, err_msg); \
  RETURN_ERROR(ret_val); \
}

// Another definition for higher level verbosity
#define RETURN_ERROR2(ret_val) { \
  if (m_tread_verbosity > 1) \
    fprintf(stderr, "%s\n", m_error_msg); \
  return ret_val; \
}

#define SET_RETURN_ERROR2(ret_val, error_val, func_name, err_msg) { \
  m_error_val = error_val; \
  SetErrorMessage(func_name, err_msg); \
  RETURN_ERROR2(ret_val); \
}

// Also define the message size (will print out error if
//  possible packet size ever becomes bigger than this value). 
#define MESSAGE_BUFFER_SIZE 1024

// defining how much space for error message
#define ERR_MSG_SIZE 128

/*****************************************************************/
/*                       Message format                          */
/*                                                               */
/*  ushort message_size                          (bytes 0 & 1)   */
/*  ushort tag                                   (bytes 1 & 2)   */
/*                                                               */
/*  Next follows the Packets, which contain two items            */
/*     char   packet_type                        (1 byte)        */
/*            Packet_Data                        (any # bytes)   */
/*                                                               */
/* The packet''s data can be anything and it knows how to handle */
/*   its own data.  If it wants to send its size it can, but is  */
/*   not necessary.  When the controller splits the message into */
/*   packets it asks each packet what size it took up (which     */
/*   is known by a packet) to determine when the next packet     */
/*   starts in the message.                                      */
/*****************************************************************/

class Net_Controller {
 public:
  // virtual destructor 
  virtual ~Net_Controller();

  // Default Constructor
  Net_Controller();

  // Register packets (returns error if mismatch, or unknown)
  //  The bit type will be in the packet itself.
  virtual int Register(Packet *) = 0;

  // Set verbosity level
  void SetVerbosity(int level) {m_tread_verbosity = level;}

  // Returns the packet type for the next message, moving the message
  //  pointer one character forward (-1 if seen full message)
  virtual int GetPacketType();

  // Starting at current message pointer it fills the return_packet
  //  data structure (gets size from return_packet->GetSize())
  //  (-1 if runs out of message space)
  virtual int GetPacketFromMsg(Packet * return_packet);

  // Print last error message, and return the error val, or just
  //  return the val w/o printing
  int PrintLastError();
  inline int ReturnLastErrorVal() { return m_error_val; }

  // These functions are implemented by the child classes, because they
  //  have different protocols.
  virtual int Send() = 0;
  virtual int Receive(BOOL wait) = 0;
  virtual int AddPacket(Packet *) = 0;
  // Notes - for both send should ensure that the parse packet parsed
  //  the whole message
  // For client the addpacket should not work till done with parsing
  //  because we need the original packets around to verify any waiting
  //  for responses have been filled, not necessary for server

 protected:
  // The functions below are the low level message buffer functions.
  // They are protected because the child controller classes will
  //  call them under the hood, while the semantics of how the network 
  //  runs will be handled in the above functions.

  // Adds packet unto message buffer, checking 
  //  that it is registered and fills correctly
  virtual int PutPacketInMsg(Packet * packet_add);

  // Send message across network with packets in buffer
  virtual int Send_Packets(int fd, ushort tag);

  // Receive message across network using specified socket
  virtual int Receive_Message(int fd, int waittime_usec);

  // Cancel current message, whether it is something preparing or
  //  something received
  virtual void Cancel_Message();

  // Get the message tag - needed by server and client controllers
  // (server to send back, client to verify)
  virtual ushort GetMessageTag();

  // Essentially this represents whether a message has been
  //  received and has not yet been fully parsed (or cancelled)
  inline BOOL MessageInUse() { return m_new_message;}

  // Calls to packet protected classes (must be done through net_controller,
  //  since that is the friend - friendship does not pass to child classes!!)
  inline void SetPacketType(Packet * pack, uint val) {pack->SetType(val);}
  inline void SetPacketRegistered(Packet * pack) {pack->SetRegistered();}

  // Setting error message
  int SetErrorMessage(char * func_name, char * err_msg);
  inline void SetErrorVal(int val) {m_error_val = val;}

  /**************************************************************/
  /*                       Class Variables                      */
  /**************************************************************/

  // Printing verbosity (for error and success messages)
  int m_tread_verbosity;

  // Error variables, if want more info (stores last error seen)
  // List of all error vals
  // PACKET_TYPE_ERROR -101
  // PACKET_READ_ERROR -102
  // PACKET_UNREGISTERED_ERROR -103
  // PACKET_FILL_ERROR -104
  int m_error_val;
  char * m_error_msg;

 private:
  // Message buffer and pointer (offset) into it
  char *m_msg_buffer;
  int m_msg_pointer;

  // Set if there is a new message_buffer received that
  //  has not been fully parsed.
  BOOL m_new_message;

  // Size of incoming message
  short m_complete_message_size;

  // Function for checking if socket is ready to be read
  // No reason to override this function so private
  int ReadNonBlocking(int fd, int microsecs);
};

#endif
