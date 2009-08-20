#include "ServerRequest_Packet.h"

char * ServerRequest_Packet::m_class_name = strdup("ServerRequest_Packet");

/***********************************************************/
/*  Name : ~ServerRequest_Packet                           */
/*                                                         */
/*  Description: Destructor (free memory specific to each  */
/*                instance of ServerRequest_Packet class)  */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

ServerRequest_Packet::~ServerRequest_Packet() {
  // nothing to do here because Packet already deletes the
  //  data, since nothing else associated with this type of packet.
}

/*********************************************************/
/*  Name : ServerRequest_Packet                          */
/*                                                       */
/*  Description: Default Constructor                     */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

ServerRequest_Packet::ServerRequest_Packet() {
  // Some default values
  m_registered = FALSE;
  m_type = 0;
  m_wait_response = TRUE;

  // Create new data packet
  m_data = new ServerRequest_Data();
}

/*********************************************************/
/*  Name : ServerRequest_Packet                          */
/*                                                       */
/*  Description: Constructor with server request data as */
/*                input                                  */
/*                                                       */
/*  Input : Reference to server request data to          */
/*           initialize m_data                           */
/*  Return Value : None                                  */
/*********************************************************/

ServerRequest_Packet::ServerRequest_Packet(const ServerRequest_Data& sr_data) {
  m_data = new ServerRequest_Data(sr_data);

  // default values
  m_registered = FALSE;
  m_type = 0;
  m_wait_response = TRUE;
}

/*********************************************************/
/*  Name : GetCopyOfData                                 */
/*                                                       */
/*  Description: Returns a copy of m_data through CopyDat*/
/*                (CopyDat must already be created and of*/
/*                 correct type)                         */
/*                                                       */
/*  Input : Pointer to server req data to return copy    */
/*           of m_data.                                  */
/*  Return Value : 1 - successful, -1 - failure          */
/*********************************************************/

int ServerRequest_Packet::GetCopyOfData(Packet_Data * CopyDat) {
  ServerRequest_Data* CopySRDat = dynamic_cast<ServerRequest_Data*>(CopyDat);
  if (CopySRDat != NULL) {
    *(CopySRDat) = *(SERVERREQUEST_DATA);
    return 1;
  }
  return -1;
}

/*********************************************************/
/*  Name : SetData                                       */
/*                                                       */
/*  Description: Sets m_data values based on NewDat      */
/*                (m_data is not deleted, just its values*/
/*                changed)                               */
/*                                                       */
/*  Input : Pointer to server req data which contains the*/
/*           values to set m_data (does not set ptr      */
/*           since that would allow the data storage to  */
/*           exist outside of this class and thus        */
/*           vunerable to being deleted)                 */
/*  Return Value : 1 - successful, -1 - failure          */
/*********************************************************/

int ServerRequest_Packet::SetData(Packet_Data * NewDat) {
  ServerRequest_Data* NewSRDat = dynamic_cast<ServerRequest_Data*>(NewDat);
  if (NewSRDat != NULL) {
    *(SERVERREQUEST_DATA) = *(NewSRDat);
    return 1;
  }
  return -1;
}

/*********************************************************/
/*  Name : Fill                                          */
/*                                                       */
/*  Description: Fill message buffer with data           */
/*                                                       */
/*  Input : pointer to beginning of portion of message   */
/*           buffer to fill data                         */
/*  Return Value : Returns how many bytes were filled in */
/*                  to verify correctness (-1 if fail)   */
/*********************************************************/

int ServerRequest_Packet::Fill(char * mb) {
  // Check that m_data is valid first
  if (m_data != NULL) {
    // Order of ServerRequest packet is
    // Request bit field (uint - 4 bytes)

    return FillInt(mb, SERVERREQUEST_DATA->GetRequestField());
  }
  return -1;
}

/*********************************************************/
/*  Name : Read Packet                                   */
/*                                                       */
/*  Description: Fill Data packet with information from  */
/*                message buffer                         */
/*                                                       */
/*  Input : pointer to beginning of portion of message   */
/*           buffer to receive data                      */
/*  Return Value : Returns success (1) or failure (-1)   */
/*********************************************************/

int ServerRequest_Packet::ReadPacket(char * mb) {
  // Update data with new information
  if (m_data != NULL) {
    // Order of command packet is
    // Request bit field (uint - 4 bytes)

    // Set data to correspond to new values
    SERVERREQUEST_DATA->SetRequestField(GrabInt(mb));
    return 1;
  }
  return -1;
}

