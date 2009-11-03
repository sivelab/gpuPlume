/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */


#ifndef __RECEIVER_H__
#define __RECEIVER_H__

#include <string>
#include "InetAddress.h"

////////////////////////////////////////////////////////////
// Receiver.h
//
// Class definition for the recipient of a broadcasted message
//

#ifndef WIN32
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#endif

class Receiver 
{
  public :
  
    Receiver();

    // Setup a receiver to connect to the broadcaster host and port.
    Receiver( std::string &hostname, short port );

    ~Receiver();

    // set host for receiver
    void setHost( const std::string &hostname );

    // setBuffer defines the buffer into which the broadcasted
    // message will be received.
    void setBuffer( void *buffer, const unsigned int size );

    // Define what port to listen and bind to
    void setPort( const short port );

    // Sync does a blocking wait to recieve next message
    void sync( void );

  private :
    bool init( void );
  
  private :
#if defined (WIN32)
    SOCKET _so;
#else
    int _so;
#endif
#if defined (WIN32)
    SOCKADDR_IN saddr;
#else
    struct sockaddr_in saddr;
#endif
    bool _initialized;
    void *_buffer;
    unsigned int _buffer_size;

    InetAddress _broadcast_host;
    short _port;
};
#endif 
