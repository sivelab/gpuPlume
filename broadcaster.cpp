#include <iostream>
#include <fcntl.h>
#include <sys/types.h>

#if !defined (WIN32)
#include <sys/ioctl.h>
#include <sys/uio.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <sys/time.h>
#include <net/if.h>
#include <netdb.h>
#endif

#include <string.h>

#if defined(__linux)
    #include <unistd.h>
    #include <linux/sockios.h>
#elif defined(__FreeBSD__)
    #include <unistd.h>
    #include <sys/sockio.h>
#elif defined(__sgi)
    #include <unistd.h>
    #include <net/soioctl.h>
#elif defined(__CYGWIN__) 
    // nothing needed
#elif defined (__DARWIN_OSX__)
    #include <unistd.h>
    #include <sys/sockio.h>
#elif defined (WIN32)
    #include <winsock.h>
    #include <stdio.h>
#elif defined (__hpux__)
    #include <unistd.h>
#else
    #error Teach me how to build on this system
#endif

#include "broadcaster.h"


Broadcaster::Broadcaster( void )
  : _eth_device( "eth0" )
{
  _port = 0;
  _initialized = false;
  _buffer = 0L;
  _address = 0;
}

Broadcaster::~Broadcaster( void )
{
#if defined (WIN32)
  closesocket( _so);
#else
  close( _so );
#endif
}

bool Broadcaster::init( void )
{
#if defined (WIN32)
  WORD version = MAKEWORD(1,1);
  WSADATA wsaData;
  // First, we start up Winsock
  WSAStartup(version, &wsaData);
#endif

  //
  // Verify that a valid port has been specified.
  //
  if( _port == 0 )
    {
      std::cerr << "Broadcaster::init() - port not defined: " << _port << std::endl;
      return false;
    }
  
  // 
  // Set up the UDP socket for broadcasting.
  //
  if( (_so = socket( AF_INET, SOCK_DGRAM, 0 )) < 0 )
    {
      perror( "Socket" );
      return false;
    }

#if defined (WIN32)
  const BOOL on = TRUE;
#else
  int on = 1;
#endif
  
  //
  // Allow us to reuse the socket address.  Makes for re-running the
  // code much more simple.
  //
#if defined (WIN32)
  setsockopt( _so, SOL_SOCKET, SO_REUSEADDR, (const char *) &on, sizeof(int));
#else
  setsockopt( _so, SOL_SOCKET, SO_REUSEADDR, &on, sizeof(on));
#endif

  //
  // Set the socket options so that we create a broadcast socket.
  //
#if defined (WIN32)
  setsockopt( _so, SOL_SOCKET, SO_BROADCAST, (const char *) &on, sizeof(int));
#else
  setsockopt( _so, SOL_SOCKET, SO_BROADCAST, &on, sizeof(on));
#endif

  bzero( &saddr, sizeof(saddr) );
  saddr.sin_family = AF_INET;
  saddr.sin_port   = htons( _port );  // specify the port as a 
                                      // short in host byte order
  saddr.sin_addr.s_addr = INADDR_ANY; // host interface

  if( _address == 0 )
    {
      
#if !defined (WIN32)
      struct ifreq ifr;
#endif

#if defined (__linux)
      strncpy( ifr.ifr_name, _eth_device.c_str(), sizeof(ifr.ifr_name) );
#elif !defined (WIN32)
      strncpy( ifr.ifr_name, "ef0", sizeof(ifr.ifr_name) );
#endif

#if defined (WIN32) // get the server address
      saddr.sin_addr.s_addr = htonl(INADDR_BROADCAST);
    }
#else

  if (ioctl( _so, SIOCGIFADDR, (void *) &ifr) < 0) {
    perror("ioctl SIOCGIFADDR");
    return false;
  }
  

      if( (ioctl( _so, SIOCGIFBRDADDR, &ifr)) < 0 )
       {
          perror( "Broadcaster::init() Cannot get Broadcast Address" );
          // need to throw expceptions and get out of the program...!
          return false;
       }
  
      saddr.sin_addr.s_addr = (((sockaddr_in *)&ifr.ifr_broadaddr)->sin_addr.s_addr);
    }
  else
    {
      // else clause??
      std::cerr << "else clause" << std::endl;
      saddr.sin_addr.s_addr = _address;
    }
#endif

#ifdef _VERBOSE
    unsigned char *ptr = (unsigned char *)&saddr.sin_addr.s_addr;
    printf( "Broadcast address : %u.%u.%u.%u\n", ptr[0], ptr[1], ptr[2], ptr[3] );
#endif

    _initialized = true;
    return _initialized;
}

void Broadcaster::setHost( const char *hostname )
{
    struct hostent *h;
    if( (h = gethostbyname( hostname )) == 0L )
    {
      std::cerr << "Broadcaster::setHost() - Cannot resolv an address for \"" << hostname << "\"." << std::endl;
      _address = 0;
    }
    else
      _address = *(( unsigned long  *)h->h_addr);
}

void Broadcaster::setEthernetDevice( const std::string& eth_device )
{
  _eth_device = eth_device;
}


void Broadcaster::setPort( const short port )
{
  _port = port;
}

void Broadcaster::setBuffer( void *buffer, const unsigned int size )
{
    _buffer = buffer;
    _buffer_size = size;
}

void Broadcaster::sync( void )
{
  if(!_initialized) init();

  if( _buffer == 0L )
    {
      std::cerr << "Broadcaster::sync() - No buffer" << std::endl;
      return;
    }

#if defined (WIN32)
    unsigned int size = sizeof( SOCKADDR_IN );
    sendto( _so, (const char *)_buffer, _buffer_size, 0, (struct sockaddr *)&saddr, size );
    int err = WSAGetLastError ();
    int *dum = (int*) _buffer;
#else
    unsigned int size = sizeof( struct sockaddr_in );
    sendto( _so, (const void *)_buffer, _buffer_size, 0, (struct sockaddr *)&saddr, size );
#endif

}

