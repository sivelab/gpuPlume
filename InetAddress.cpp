/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */

#include <iostream>

#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>

#include "InetAddress.h"

std::istream& operator>>( std::istream &is, InetAddress& )
{
    return is;
}

std::ostream& operator<<( std::ostream &os, const InetAddress &ia )
{
    return os << ia._ip_addr;
}

InetAddress::InetAddress()
    : _MAX_HOSTADDR_LENGTH(1024), _ip_netaddr(4)
{
}

InetAddress::InetAddress( std::string host )
    : _MAX_HOSTADDR_LENGTH(1024), _ip_netaddr(4)
{
    _host_name = host;
    _init_address();
}

InetAddress::InetAddress( InetAddress &inet )
    : _MAX_HOSTADDR_LENGTH(1024), _ip_netaddr(4)
{
    _host_name = inet.getHostName();
    _init_address();
}

InetAddress InetAddress::operator=( InetAddress &inet )
{
    _host_name = inet.getHostName();
    _init_address();
  
    return *this;
}

InetAddress InetAddress::operator=( const std::string& host )
{
//    std::cout << "InetAddress: " << host << std::endl;

    _host_name = host;
    _init_address();
  
    return *this;
}

void InetAddress::_init_address( void ) 
{
    // struct hostent {
    //   char    *h_name;        /* official name of host */
    //   char    **h_aliases;    /* alias list */
    //   int     h_addrtype;     /* host address type */
    //   int     h_length;       /* length of address */
    //   char    **h_addr_list;  /* list of addresses */
    // }
    register struct hostent *hostptr;

    hostptr = gethostbyname( _host_name.c_str() );
    if (hostptr == 0) {
        throw UnknownHostException();
    }
  
#if 0   // debugging info
    std::cout << "\taliases: ";
    int d_i=0;
    while (hostptr->h_aliases[d_i])
        std::cout << " " << hostptr->h_aliases[d_i++];
    std::cout << std::endl;
    std::cout << "\taddr type = " << hostptr->h_addrtype << ", addr length = " << hostptr->h_length << std::endl;
#endif

    //
    // store the ip net address in the class for later reference
    // - get the first ip address off the h_addr_list
    //
    int i;
    for(i=0; i<hostptr->h_length; i++)
        _ip_netaddr[i] = hostptr->h_addr_list[0][i];
  
    // store the network address in string form
    _ip_addr = inet_ntoa( *(struct in_addr*)(hostptr->h_addr_list[0]) );

#if 0   // debugging info
    char **listptr = hostptr->h_addr_list;
    struct in_addr *ptr;

    switch( hostptr->h_addrtype ) {
        case AF_INET:
    
            while ( (ptr = (struct in_addr*) *listptr++) != 0 ) 
                std::cout << "\tInternet address: " << inet_ntoa(*ptr) << '.';
            std::cout << std::endl;
            break;

        default:
            std::cerr << "unknown address type!" << std::endl;
    }
#endif

}

/** Determines all the IP addresses of a host, given the host's name */
// vector<InetAddress> InetAddress::getAllByName( string host )
// {
//   vector<InetAddress> ret_vec(5);
// }


std::vector<char> InetAddress::getAddress() 
{
    return _ip_netaddr;
}

std::string InetAddress::getHostAddress() 
{
    return _ip_addr;
}

std::string InetAddress::getBroadcastAddress() 
{
    std::string bc_addr;
    return bc_addr;
}

std::string InetAddress::getHostName()
{
    return _host_name;
}

// bool InetAddress::isMulticastAddress()
//{
//   return false;
// }
