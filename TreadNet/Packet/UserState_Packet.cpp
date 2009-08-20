#include "UserState_Packet.h"

char * UserState_Packet::m_class_name = strdup("UserState_Packet");

/***********************************************************/
/*  Name : ~UserState_Packet                               */
/*                                                         */
/*  Description: Destructor (free memory specific to each  */
/*                instance of UserState_Packet class)      */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

UserState_Packet::~UserState_Packet() {
  // nothing to do here because Packet already deletes the
  //  data, since nothing else associated with this type of packet.
}

/*********************************************************/
/*  Name : UserState_Packet                                 */
/*                                                       */
/*  Description: Default Constructor                     */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

UserState_Packet::UserState_Packet() {
  // Some default values
  m_registered = FALSE;
  m_type = 0;
  m_wait_response = FALSE;

  // Create new data packet
  m_data = new UserState_Data();
}

/*********************************************************/
/*  Name : UserState_Packet                              */
/*                                                       */
/*  Description: Constructor with UserState data as input*/
/*                                                       */
/*  Input : Reference to UserState data to initialize    */
/*           m_data                                      */
/*  Return Value : None                                  */
/*********************************************************/

UserState_Packet::UserState_Packet(const UserState_Data& us_data) {
  m_data = new UserState_Data(us_data);

  // default values
  m_type = 0;
  m_wait_response = FALSE;
  m_registered = FALSE;
}

/*********************************************************/
/*  Name : GetCopyOfData                                 */
/*                                                       */
/*  Description: Returns a copy of m_data through CopyDat*/
/*                (CopyDat must already be created and of*/
/*                 correct type)                         */
/*                                                       */
/*  Input : Pointer to UserState data to return copy     */
/*           of m_data.                                  */
/*  Return Value : 1 - successful, -1 - failure          */
/*********************************************************/

int UserState_Packet::GetCopyOfData(Packet_Data * CopyDat) {
  UserState_Data* CopyUSDat = dynamic_cast<UserState_Data*>(CopyDat);
  if (CopyUSDat != NULL) {
    *(CopyUSDat) = *(USERSTATE_DATA);
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
/*  Input : Pointer to UserState data which contains the */
/*           values to set m_data (does not set ptr      */
/*           since that would allow the data storage to  */
/*           exist outside of this class and thus        */
/*           vunerable to being deleted)                 */
/*  Return Value : 1 - successful, -1 - failure          */
/*********************************************************/

int UserState_Packet::SetData(Packet_Data * NewDat) {
  UserState_Data* NewUSDat = dynamic_cast<UserState_Data*>(NewDat);
  if (NewUSDat != NULL) {
    *(USERSTATE_DATA) = *(NewUSDat);
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

int UserState_Packet::Fill(char * mb) {
  // Check that m_data is valid first
  if (m_data != NULL) {
    // Pointer to current place in buffer
    char * place_ptr = mb;
    int ptr_add = 0;
    
    // Order of UserState packet is
    // Position (East then North) - 2 doubles, 8 bytes
    // User Offset (E/N/Elev) - 3 doubles, 12 bytes
    // Eye Offset (E/N/Elev) - 3 doubles, 12 bytes
    // Facing (E/N) - 2 doubles, 8 bytes
    // Velocity (E/N/Vert) - 3 doubles, 12 bytes
    // Status - int, 4 bytes
    // Position Valid Flag - char, 1 byte

    ptr_add  = FillDouble(place_ptr, USERSTATE_DATA->GetPosEasting());
    ptr_add += FillDouble(place_ptr + ptr_add, USERSTATE_DATA->GetPosNorthing());

    ptr_add += FillTRPoint(place_ptr + ptr_add, USERSTATE_DATA->GetUserOffset());
    ptr_add += FillTRPoint(place_ptr + ptr_add, USERSTATE_DATA->GetEyeOffset());

    ptr_add += FillDouble(place_ptr + ptr_add, USERSTATE_DATA->GetFacingEasting());
    ptr_add += FillDouble(place_ptr + ptr_add, USERSTATE_DATA->GetFacingNorthing());
    
    ptr_add += FillTRVector(place_ptr + ptr_add, USERSTATE_DATA->GetVelocity());
    ptr_add += FillInt(place_ptr + ptr_add, USERSTATE_DATA->GetStatus());
    ptr_add += FillChar(place_ptr + ptr_add, USERSTATE_DATA->PosValid());
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

int UserState_Packet::ReadPacket(char * mb) {
  // Update data with new information
  if (m_data != NULL) {
    // Pointer to current place in buffer
    char * place_ptr = mb;

    // Order of UserState packet is
    // Position (East then North) - 2 doubles, 8 bytes
    // User Offset (E/N/Elev) - 3 doubles, 12 bytes
    // Eye Offset (E/N/Elev) - 3 doubles, 12 bytes
    // Facing (E/N) - 2 doubles, 8 bytes
    // Velocity (E/N/Vert) - 3 doubles, 12 bytes
    // Status - int, 4 bytes
    // Position Valid Flag - char, 1 byte

    // Grab and set position values
    USERSTATE_DATA->SetPosEasting(GrabDouble(place_ptr));
    place_ptr += sizeof(double);
    USERSTATE_DATA->SetPosNorthing(GrabDouble(place_ptr));
    place_ptr += sizeof(double);

    // Note TRVector and TRPoint are the same
    TRVector temp;

    // Grab and set User Offset
    GrabTRPoint(place_ptr, temp);
    USERSTATE_DATA->SetUserOffset(temp);
    place_ptr += temp.GetSize();

    // Grab and set Eye Offset
    GrabTRPoint(place_ptr, temp);
    USERSTATE_DATA->SetEyeOffset(temp);
    place_ptr += temp.GetSize();

    // Grab and set facing values
    USERSTATE_DATA->SetFacingEasting(GrabDouble(place_ptr));
    place_ptr += sizeof(double);
    USERSTATE_DATA->SetFacingNorthing(GrabDouble(place_ptr));
    place_ptr += sizeof(double);

    // Grab and set Velocity Offset
    GrabTRVector(place_ptr, temp);
    USERSTATE_DATA->SetVelocity(temp);
    place_ptr += temp.GetSize();

    // Grab and set status
    USERSTATE_DATA->SetStatus((TRStatus)GrabInt(place_ptr));
    place_ptr += sizeof(int);

    // Grab and set Position Valid Flag
    USERSTATE_DATA->SetPosFlag(*place_ptr);
    return 1;
  }
  return -1;
}

