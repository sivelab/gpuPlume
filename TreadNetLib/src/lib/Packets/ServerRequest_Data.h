#ifndef SERVERREQUEST_DATA_H
#define SERVERREQUEST_DATA_H

#include "Packet.h"
#include "Packet_Data.h"

// This class is for ServerRequest packet data
//  There is nothing in the data, only the type matters
//  so all the functions (almost) can be null

class ServerRequest_Data : public Packet_Data {
 public:
  friend class ServerRequest_Packet;

  // Destructor
  ~ServerRequest_Data();
  
  // default constructor
  ServerRequest_Data();

  // copy constructor
  ServerRequest_Data(const ServerRequest_Data& right);

  // Assignment operator
  void operator = (const ServerRequest_Data&);

  // Accessor Functions
  // Allow user to remove and add packets to be requested
  int AddPacketRequest(Packet * req_pkt);
  int RemovePacketRequest(Packet * req_pkt);

  // if want to reset field
  inline void ZeroRequestField() { m_request_field = 0;}
  
  // ask if packet is in request field (1 - yes, 0 - no)
  int IsPacketRequested(Packet * req_pkt);

  // Already provides uint GetSize();
 protected:
  // Protected set and get full bit field functions (for ServerRequest_Packet)
  inline void SetRequestField(uint new_val) { m_request_field = new_val; }
  inline uint GetRequestField() const       { return m_request_field;    }

  // Inherits variable ushort m_size

 private:
  // Request bit field (each registered packet has a number assigned to it,
  //  set that bit in this field to request it from server).  Since uint
  //  it can handle up to 32 registered packets - shouldn't ever be an issue.
  uint m_request_field;
};

#endif
