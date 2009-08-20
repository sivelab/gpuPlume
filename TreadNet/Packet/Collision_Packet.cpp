#include "Collision_Packet.h"

char * Collision_Packet::m_class_name = strdup("Collision_Packet");

/***********************************************************/
/*  Name : ~Collision_Packet                               */
/*                                                         */
/*  Description: Destructor (free memory specific to each  */
/*                instance of Collision_Packet class)      */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

Collision_Packet::~Collision_Packet() {
  // nothing to do here because Packet already deletes the
  //  data, since nothing else associated with this type of packet.
}

/*********************************************************/
/*  Name : Collision_Packet                              */
/*                                                       */
/*  Description: Default Constructor                     */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

Collision_Packet::Collision_Packet() {
  // Some default values
  m_registered = FALSE;
  m_type = 0;
  m_wait_response = FALSE;

  // Create new data packet
  m_data = new Collision_Data();
}

/*********************************************************/
/*  Name : Collision_Packet                              */
/*                                                       */
/*  Description: Constructor with Collision data as      */
/*                input                                  */
/*                                                       */
/*  Input : Reference to Collision data to initialize    */
/*           m_data                                      */
/*  Return Value : None                                  */
/*********************************************************/

Collision_Packet::Collision_Packet(const Collision_Data& c_data) {
  m_data = new Collision_Data(c_data);

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
/*  Input : Pointer to Collision data to return copy     */
/*           of m_data.                                  */
/*  Return Value : 1 - successful, -1 - failure          */
/*********************************************************/

int Collision_Packet::GetCopyOfData(Packet_Data * CopyDat) {
  Collision_Data* CopyCDat = dynamic_cast<Collision_Data*>(CopyDat);
  if (CopyCDat != NULL) {
    *(CopyCDat) = *(COLL_DATA);
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
/*  Input : Pointer to Collision data which contains     */
/*           the values to set m_data (does not set ptr  */
/*           since that would allow the data storage to  */
/*           exist outside of this class and thus        */
/*           vunerable to being deleted)                 */
/*  Return Value : 1 - successful, -1 - failure          */
/*********************************************************/

int Collision_Packet::SetData(Packet_Data * NewDat) {
  Collision_Data* NewCDat = dynamic_cast<Collision_Data*>(NewDat);
  if (NewCDat != NULL) {
    *(COLL_DATA) = *(NewCDat);
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

int Collision_Packet::Fill(char * mb) {
  // Check that m_data is valid first
  if (m_data != NULL) {
    // Pointer to current place in buffer
    char * place_ptr = mb;
    int ptr_add = 0;
    
    // Order of Collision packet is
    // Collision Flag - char, 1 byte
    // Collision Normal (E/N/Up) - 3 doubles, 12 bytes
    // Collision Pos (E/N/Elev) - 3 doubles, 12 bytes

    ptr_add  = FillChar(place_ptr, (char)COLL_DATA->GetCollisionFlag());
    ptr_add += FillTRVector(place_ptr + ptr_add, COLL_DATA->GetCollisionNormal());
    ptr_add += FillTRPoint(place_ptr + ptr_add, COLL_DATA->GetCollisionPos());
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

int Collision_Packet::ReadPacket(char * mb) {
  // Update data with new information
  if (m_data != NULL) {
    // Pointer to current place in buffer
    char * place_ptr = mb;

    // Order of Collision packet is
    // Collision Flag - char, 1 byte
    // Collision Normal (E/N/Up) - 3 doubles, 12 bytes
    // Collision Pos (E/N/Elev) - 3 doubles, 12 bytes

    // Grab and set Collision Flag
    COLL_DATA->SetCollisionFlag((TRCollFlag)*place_ptr);
    place_ptr += sizeof(char);

    TRVector temp;

    // Grab and set collision normal
    GrabTRVector(place_ptr, temp);
    COLL_DATA->SetCollisionNormal(temp);
    place_ptr += temp.GetSize();

    // Grab and set collision pos
    GrabTRPoint(place_ptr, temp);
    COLL_DATA->SetCollisionPos(temp);
    return 1;
  }
  return -1;
}

