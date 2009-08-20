#include "ServerRequest_Data.h"

/***********************************************************/
/*  Name : ~ServerRequest_Data                             */
/*                                                         */
/*  Description: Destructor (free memory specific to each  */
/*                instance of ServerRequest_Data class)    */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

ServerRequest_Data::~ServerRequest_Data() {
}

/*********************************************************/
/*  Name : ServerRequest_Data                            */
/*                                                       */
/*  Description: Default Constructor                     */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

ServerRequest_Data::ServerRequest_Data() {
  // Set to default values (no requests)
  m_request_field = 0;

  m_size = sizeof(int);
}

/*********************************************************/
/*  Name : ServerRequest_Data                            */
/*                                                       */
/*  Description: Copy Constructor                        */
/*                                                       */
/*  Input : ServerRequest_Data&                          */
/*  Return Value : None                                  */
/*********************************************************/

ServerRequest_Data::ServerRequest_Data(const ServerRequest_Data& right) 
	: Packet_Data(right) {
  m_request_field = right.GetRequestField();

  m_size = sizeof(int);
}

/*********************************************************/
/*  Name : Operator =                                    */
/*                                                       */
/*  Description: Equals operator                         */
/*                                                       */
/*  Input : ServerRequest_Data&                          */
/*  Return Value : None                                  */
/*********************************************************/

void ServerRequest_Data::operator = (const ServerRequest_Data& right) {
  m_request_field = right.GetRequestField();
}

/*********************************************************/
/*  Name : AddPacketRequest                              */
/*                                                       */
/*  Description: Add packet to be requested              */
/*    Can fail if not registered, or type exceeds 32     */
/*                                                       */
/*  Input : Packet*                                      */
/*  Return Value : Success (1) / Fail (-1)               */
/*********************************************************/

int ServerRequest_Data::AddPacketRequest(Packet * req_pkt) {
  if (req_pkt->IsRegistered()) {
    uint pkt_type = req_pkt->GetType();
    if (pkt_type < 33 && pkt_type > 0) {
      // shift bit over till in correct place (note bits are numbered
      //  1 to 32) and or it to place a 1 in the request field.
      m_request_field |= (1 << (pkt_type - 1));
      return 1;
    }
  }
  return -1;
}

/*********************************************************/
/*  Name : RemovePacketRequest                           */
/*                                                       */
/*  Description: Remove packet to be requested           */
/*    Can fail if not registered, or type exceeds 32     */
/*    If not already added still returns 1 since is gone */
/*                                                       */
/*  Input : Packet*                                      */
/*  Return Value : Success (1) / Fail (-1)               */
/*********************************************************/

int ServerRequest_Data::RemovePacketRequest(Packet * req_pkt) {
  if (req_pkt->IsRegistered()) {
    uint pkt_type = req_pkt->GetType();
    if (pkt_type < 33 && pkt_type > 0) {
      // shift bit over till in correct place (note bits are numbered
      //  1 to 32), negate to make it zero (rest 1s), and and it to 
      //  zero out that bit in the field.
      m_request_field &= ~(1 << (pkt_type - 1));
      return 1;
    }
  }
  return -1;  
}

/*********************************************************/
/*  Name : IsPacketRequested                             */
/*                                                       */
/*  Description: Check if packet is being requested (i.e.*/
/*    in bit field.                                      */
/*                                                       */
/*  It will also return no if the packet is unregistered */
/*                                                       */
/*  Input : Packet*                                      */
/*  Return Value : Yes (1) / No (0)                      */
/*********************************************************/

int ServerRequest_Data::IsPacketRequested(Packet * req_pkt) {
  if (req_pkt->IsRegistered()) {
    uint pkt_type = req_pkt->GetType();
    if (pkt_type < 33 && pkt_type > 0) {
      // shift bit over till in correct place (note bits are numbered
      //  1 to 32), and it with the field, and if the result is non zero
      //  that bit was set, otherwise not
      uint test = m_request_field & (1 << (pkt_type - 1));
      if (test) {
	return 1;
      }
    }
  }
  return 0;  
}
