/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */


#ifndef __INET_ADDRESS_H__
#define __INET_ADDRESS_H__

#include <string>
#include <vector>

/** This class represents an Internet Protocol(IP) address.  The class
    interface has been adopted from the Java jdk1.2 API java.net
    class interfaces. */

class InetAddress {
  public:
    //
    // EXCEPTIONS
    //
    class UnknownHostException {};

    /** Constructs an an InetAddress based on the local host */
    InetAddress();

    /** Constructs an InetAddress based on the host's name */
    InetAddress( std::string host );

    /** Constructs an InetAddress based on an InetAddress (aka. Copy
        Constructor) */
    InetAddress( InetAddress& inet );
  
    /** Assignment operator */
    InetAddress operator=( InetAddress& inet );
    InetAddress operator=( const std::string& host );

    /** Determines all the IP addresses of a host, given the host's name */
    // static vector<InetAddress> getAllByName( string host );

    /** Returns the raw IP address in network byte order of this
        InetAddress object.  The highest order byte of the address is in
        slot [0] of the return vector */
    std::vector<char> getAddress();

    /** Returns the IP address string "%d.%d.%d.%d" */
    std::string getHostAddress();

    /** Returns the IP address string for broadcast "%d.%d.%d.255" */
    std::string getBroadcastAddress();
  
    /** Returns the hostname for this address */
    std::string getHostName();

    /** Utility routine to check if the InetAddress is an IP multicast address */
    // bool isMulticastAddress();

    friend std::istream& operator>>( std::istream&, InetAddress& );
    friend std::ostream& operator<<( std::ostream&, const InetAddress& );

  private:
    const int _MAX_HOSTADDR_LENGTH;

    std::vector<char> _ip_netaddr;
    std::string _ip_addr;
    std::string _host_name;

    void _init_address( void );

};

#endif // __INET_ADDRESS_H__
