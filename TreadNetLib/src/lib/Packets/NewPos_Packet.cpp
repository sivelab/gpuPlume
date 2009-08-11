#include "NewPos_Packet.h"

char * NewPos_Packet::m_class_name = strdup("NewPos_Packet");

/***********************************************************/
/*  Name : ~NewPos_Packet                                  */
/*                                                         */
/*  Description: Destructor (free memory specific to each  */
/*                instance of NewPos_Packet class)         */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

NewPos_Packet::~NewPos_Packet() {
  // nothing to do here because Packet already deletes the
  //  data, since nothing else associated with this type of packet.
}

/*********************************************************/
/*  Name : NewPos_Packet                                 */
/*                                                       */
/*  Description: Default Constructor                     */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

NewPos_Packet::NewPos_Packet() {
  // Some default values
  m_registered = FALSE;
  m_type = 0;
  m_wait_response = TRUE;

  // Create new data packet
  m_data = new NewPos_Data();
}

/*********************************************************/
/*  Name : NewPos_Packet                                 */
/*                                                       */
/*  Description: Constructor with NewPos data as         */
/*                input                                  */
/*                                                       */
/*  Input : Reference to NewPos data to initialize       */
/*           m_data                                      */
/*  Return Value : None                                  */
/*********************************************************/

NewPos_Packet::NewPos_Packet(const NewPos_Data& np_data) {
  m_data = new NewPos_Data(np_data);

  // default values
  m_registered = FALSE;
  m_type = 0;
  m_wait_response = TRUE;
}

/*********************************************************/
/*  Name : NewPos_Packet                                 */
/*                                                       */
/*  Description: Constructor with Virtual and User state */
/*                 data as input                         */
/*                                                       */
/*  Input : Reference to Virtual and User state data to  */
/*           initialize m_data                           */
/*  Return Value : None                                  */
/*********************************************************/

NewPos_Packet::NewPos_Packet(const UserState_Data& us_data, const VirtualState_Data& vs_data) {
  m_data = new NewPos_Data(us_data, vs_data);

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
/*  Input : Pointer to NewPos data to return copy        */
/*           of m_data.                                  */
/*  Return Value : 1 - successful, -1 - failure          */
/*********************************************************/

int NewPos_Packet::GetCopyOfData(Packet_Data * CopyDat) {
  NewPos_Data* CopyNPDat = dynamic_cast<NewPos_Data*>(CopyDat);
  if (CopyNPDat != NULL) {
    *(CopyNPDat) = *(NEWPOS_DATA);
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
/*  Input : Pointer to NewPos data which contains        */
/*           the values to set m_data (does not set ptr  */
/*           since that would allow the data storage to  */
/*           exist outside of this class and thus        */
/*           vunerable to being deleted)                 */
/*  Return Value : 1 - successful, -1 - failure          */
/*********************************************************/

int NewPos_Packet::SetData(Packet_Data * NewDat) {
  NewPos_Data* NewNPDat = dynamic_cast<NewPos_Data*>(NewDat);
  if (NewNPDat != NULL) {
    *(NEWPOS_DATA) = *(NewNPDat);
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

int NewPos_Packet::Fill(char * mb) {
  // Check that m_data is valid first
  if (m_data != NULL) {
    // Pointer to current place in buffer
    char * place_ptr = mb;
    int ptr_add = 0;
    
    // Order of NewPos packet is
    // User State
    // Virtual State

    // To make code simpler I will have to create a couple dummy
    //  packets (could not make NewPos consist of these packets, unless
    //  I always wanted duplicate data)
    UserState_Packet usp(GetUserState());
    VirtualState_Packet vsp(GetVirtualState());

    ptr_add  = usp.Fill(place_ptr);
    ptr_add += vsp.Fill(place_ptr + ptr_add);
    return ptr_add;
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

int NewPos_Packet::ReadPacket(char * mb) {
  // Update data with new information
  if (m_data != NULL) {
    // Pointer to current place in buffer
    char * place_ptr = mb;

    // Order of NewPos packet is
    // User State
    // Virtual State

    // To make code simpler I will have to create a couple dummy
    //  packets and data (could not make NewPos consist of these packets, unless
    //  I always wanted duplicate data)
    UserState_Packet usp;
    VirtualState_Packet vsp;
    UserState_Data usd;
    VirtualState_Data vsd;
    
    // Get and Set user state
    usp.ReadPacket(place_ptr);
    usp.GetCopyOfData(&usd);
    NEWPOS_DATA->SetUserState(usd);
    place_ptr += usp.GetSize();

    // Get and Set virtual state
    vsp.ReadPacket(place_ptr);
    vsp.GetCopyOfData(&vsd);
    NEWPOS_DATA->SetVirtualState(vsd);
    return 1;
  }
  return -1;
}

