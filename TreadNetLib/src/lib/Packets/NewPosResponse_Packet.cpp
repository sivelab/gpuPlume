#include "NewPosResponse_Packet.h"

char * NewPosResponse_Packet::m_class_name = strdup("NewPosResponse_Packet");

/***********************************************************/
/*  Name : ~NewPosResponse_Packet                          */
/*                                                         */
/*  Description: Destructor (free memory specific to each  */
/*                instance of NewPosResponse_Packet class) */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

NewPosResponse_Packet::~NewPosResponse_Packet() {
  // nothing to do here because Packet already deletes the
  //  data, since nothing else associated with this type of packet.
}

/*********************************************************/
/*  Name : NewPosResponse_Packet                         */
/*                                                       */
/*  Description: Default Constructor                     */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

NewPosResponse_Packet::NewPosResponse_Packet() {
  // Some default values
  m_registered = FALSE;
  m_type = 0;
  m_wait_response = FALSE;

  // Create new data packet
  m_data = new NewPosResponse_Data();
}

/*********************************************************/
/*  Name : NewPosResponse_Packet                         */
/*                                                       */
/*  Description: Constructor with NewPosResponse data as */
/*                input                                  */
/*                                                       */
/*  Input : Reference to NewPosResponse data to          */
/*           initialize m_data                           */
/*  Return Value : None                                  */
/*********************************************************/

NewPosResponse_Packet::NewPosResponse_Packet(const NewPosResponse_Data& npr_data) {
  m_data = new NewPosResponse_Data(npr_data);

  // default values
  m_registered = FALSE;
  m_type = 0;
  m_wait_response = FALSE;
}

/*********************************************************/
/*  Name : GetCopyOfData                                 */
/*                                                       */
/*  Description: Returns a copy of m_data through CopyDat*/
/*                (CopyDat must already be created and of*/
/*                 correct type)                         */
/*                                                       */
/*  Input : Pointer to NewPosResponse data to return copy*/
/*           of m_data.                                  */
/*  Return Value : 1 - successful, -1 - failure          */
/*********************************************************/

int NewPosResponse_Packet::GetCopyOfData(Packet_Data * CopyDat) {
  NewPosResponse_Data* CopyNPRDat = dynamic_cast<NewPosResponse_Data*>(CopyDat);
  if (CopyNPRDat != NULL) {
    *(CopyNPRDat) = *(NEWPOSRESP_DATA);
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
/*  Input : Pointer to NewPosResponse data which contains*/
/*           the values to set m_data (does not set ptr  */
/*           since that would allow the data storage to  */
/*           exist outside of this class and thus        */
/*           vunerable to being deleted)                 */
/*  Return Value : 1 - successful, -1 - failure          */
/*********************************************************/

int NewPosResponse_Packet::SetData(Packet_Data * NewDat) {
  NewPosResponse_Data* NewNPRDat = dynamic_cast<NewPosResponse_Data*>(NewDat);
  if (NewNPRDat != NULL) {
    *(NEWPOSRESP_DATA) = *(NewNPRDat);
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

int NewPosResponse_Packet::Fill(char * mb) {
  // Check that m_data is valid first
  if (m_data != NULL) {
    // Order of NewPosResponse packet is
    //  flag (char - 1 byte)

    return FillChar(mb, (char)NEWPOSRESP_DATA->GetFlag());
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

int NewPosResponse_Packet::ReadPacket(char * mb) {
  // Update data with new information
  if (m_data != NULL) {
    // Order of NewPosResponse packet is
    //  flag (char - 1 byte)

    // Grab and set flag value
    NEWPOSRESP_DATA->SetFlag((TRNewPosFlag)*mb);
    return 1;
  }
  return -1;
}

