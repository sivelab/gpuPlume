#include "Stairs_Packet.h"

char * Stairs_Packet::m_class_name = strdup("Stairs_Packet");

/***********************************************************/
/*  Name : ~Stairs_Packet                                  */
/*                                                         */
/*  Description: Destructor (free memory specific to each  */
/*                instance of Stairs_Packet class)         */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

Stairs_Packet::~Stairs_Packet() {
  // nothing to do here because Packet already deletes the
  //  data, since nothing else associated with this type of packet.
}

/*********************************************************/
/*  Name : Stairs_Packet                                 */
/*                                                       */
/*  Description: Default Constructor                     */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

Stairs_Packet::Stairs_Packet() {
  // Some default values
  m_registered = FALSE;
  m_type = 0;
  m_wait_response = FALSE;

  // Create new data packet
  m_data = new Stairs_Data();
}

/*********************************************************/
/*  Name : Stairs_Packet                                 */
/*                                                       */
/*  Description: Constructor with stairs data as input   */
/*                                                       */
/*  Input : Reference to stairs data to initialize m_data*/
/*  Return Value : None                                  */
/*********************************************************/

Stairs_Packet::Stairs_Packet(const Stairs_Data& stair_data) {
  m_data = new Stairs_Data(stair_data);

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
/*  Input : Pointer to stairs data to return copy        */
/*           of m_data.                                  */
/*  Return Value : 1 - successful, -1 - failure          */
/*********************************************************/

int Stairs_Packet::GetCopyOfData(Packet_Data * CopyDat) {
  Stairs_Data* CopyStairDat = dynamic_cast<Stairs_Data*>(CopyDat);
  if (CopyStairDat != NULL) {
    *(CopyStairDat) = *(STAIR_DATA);
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
/*  Input : Pointer to stairs data which contains the    */
/*           values to set m_data (does not set ptr      */
/*           since that would allow the data storage to  */
/*           exist outside of this class and thus        */
/*           vunerable to being deleted)                 */
/*  Return Value : 1 - successful, -1 - failure          */
/*********************************************************/

int Stairs_Packet::SetData(Packet_Data * NewDat) {
  Stairs_Data* NewStairDat = dynamic_cast<Stairs_Data*>(NewDat);
  if (NewStairDat != NULL) {
    *(STAIR_DATA) = *(NewStairDat);
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

int Stairs_Packet::Fill(char * mb) {
  // Check that m_data is valid first
  if (m_data != NULL) {
    // Pointer to current place in buffer
    char * place_ptr = mb;
    int ptr_add = 0;
    
    // Order of stairs packet is
    //  Stair Cmd (char - 1 byte)
    //  Stair Dir (char - 1 byte)

    ptr_add  = FillChar(place_ptr, (char)STAIR_DATA->GetCmd());
    ptr_add += FillChar(place_ptr + ptr_add, (char)STAIR_DATA->GetStairsDir());
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

int Stairs_Packet::ReadPacket(char * mb) {
  // Update data with new information
  if (m_data != NULL) {
    // Pointer to current place in buffer
    char * place_ptr = mb;

    // Order of stairs packet is
    //  Stair Cmd (char - 1 byte)
    //  Stair Dir (char - 1 byte)

    // Grab and set command
    STAIR_DATA->SetCmd((TRStairs_Cmd)*place_ptr);
    place_ptr += sizeof(char);

    // Grab and set direction
    STAIR_DATA->SetStairsDir((TRStairs_Dir)*place_ptr);
    return 1;
  }
  return -1;
}

