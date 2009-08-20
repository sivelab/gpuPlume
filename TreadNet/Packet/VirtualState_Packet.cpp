#include "VirtualState_Packet.h"

char * VirtualState_Packet::m_class_name = strdup("VirtualState_Packet");

/***********************************************************/
/*  Name : ~VirtualState_Packet                            */
/*                                                         */
/*  Description: Destructor (free memory specific to each  */
/*                instance of VirtualState_Packet class)   */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

VirtualState_Packet::~VirtualState_Packet() {
  // nothing to do here because Packet already deletes the
  //  data, since nothing else associated with this type of packet.
}

/*********************************************************/
/*  Name : VirtualState_Packet                           */
/*                                                       */
/*  Description: Default Constructor                     */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

VirtualState_Packet::VirtualState_Packet() {
  // Some default values
  m_registered = FALSE;
  m_type = 0;
  m_wait_response = FALSE;

  // Create new data packet
  m_data = new VirtualState_Data();
}

/*********************************************************/
/*  Name : VirtualState_Packet                           */
/*                                                       */
/*  Description: Constructor with VirtualState data as   */
/*                input                                  */
/*                                                       */
/*  Input : Reference to VirtualState data to initialize */
/*           m_data                                      */
/*  Return Value : None                                  */
/*********************************************************/

VirtualState_Packet::VirtualState_Packet(const VirtualState_Data& vs_data) {
  m_data = new VirtualState_Data(vs_data);

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
/*  Input : Pointer to VirtualState data to return copy  */
/*           of m_data.                                  */
/*  Return Value : 1 - successful, -1 - failure          */
/*********************************************************/

int VirtualState_Packet::GetCopyOfData(Packet_Data * CopyDat) {
  VirtualState_Data* CopyVSDat = dynamic_cast<VirtualState_Data*>(CopyDat);
  if (CopyVSDat != NULL) {
    *(CopyVSDat) = *(VIRTSTATE_DATA);
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
/*  Input : Pointer to VirtualState data which contains  */
/*           the values to set m_data (does not set ptr  */
/*           since that would allow the data storage to  */
/*           exist outside of this class and thus        */
/*           vunerable to being deleted)                 */
/*  Return Value : 1 - successful, -1 - failure          */
/*********************************************************/

int VirtualState_Packet::SetData(Packet_Data * NewDat) {
  VirtualState_Data* NewVSDat = dynamic_cast<VirtualState_Data*>(NewDat);
  if (NewVSDat != NULL) {
    *(VIRTSTATE_DATA) = *(NewVSDat);
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

int VirtualState_Packet::Fill(char * mb) {
  // Check that m_data is valid first
  if (m_data != NULL) {
    // Pointer to current place in buffer
    char * place_ptr = mb;
    int ptr_add = 0;
    
    // Order of VirtualState packet is
    // Surface Normal (E/N/Up) - 3 doubles, 12 bytes
    // Terrain - char, 1 byte
    // Virtual Force (E/N/Up) - 3 doubles, 12 bytes
    // Movement scale - float, 4 bytes
    // User Scale - float, 4 bytes

    ptr_add  = FillTRVector(place_ptr, VIRTSTATE_DATA->GetSurfaceNormal());
    ptr_add += FillChar(place_ptr + ptr_add, (char)VIRTSTATE_DATA->GetTerrain());
    ptr_add += FillTRVector(place_ptr + ptr_add, VIRTSTATE_DATA->GetVirtualForce());
    ptr_add += FillFloat(place_ptr + ptr_add, VIRTSTATE_DATA->GetMovementScale());
    ptr_add += FillFloat(place_ptr + ptr_add, VIRTSTATE_DATA->GetUserScale());
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

int VirtualState_Packet::ReadPacket(char * mb) {
  // Update data with new information
  if (m_data != NULL) {
    // Pointer to current place in buffer
    char * place_ptr = mb;

    // Order of VirtualState packet is
    // Surface Normal (E/N/Up) - 3 doubles, 12 bytes
    // Terrain - char, 1 byte
    // Virtual Force (E/N/Up) - 3 doubles, 12 bytes
    // Movement scale - float, 4 bytes
    // User Scale - float, 4 bytes

    // Grab and set surface normal
    TRVector temp;

    GrabTRVector(place_ptr, temp);
    VIRTSTATE_DATA->SetSurfaceNormal(temp);
    place_ptr+=temp.GetSize();

    // Grab and set Terrain
    VIRTSTATE_DATA->SetTerrain((TRTerrain)*place_ptr);
    place_ptr += sizeof(char);

    // Grab and set virtual force
    GrabTRVector(place_ptr, temp);
    VIRTSTATE_DATA->SetVirtualForce(temp);
    place_ptr+=temp.GetSize();

    // Grab/Set movement and user scale
    VIRTSTATE_DATA->SetMovementScale(GrabFloat(place_ptr));
    place_ptr += sizeof(float);
    VIRTSTATE_DATA->SetUserScale(GrabFloat(place_ptr));
    return 1;
  }
  return -1;
}

