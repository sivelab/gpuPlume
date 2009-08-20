#include "Registration_Packet.h"

char * Registration_Packet::m_class_name = strdup("Registration_Packet");

/***********************************************************/
/*  Name : ~Registration_Packet                            */
/*                                                         */
/*  Description: Destructor (free memory specific to each  */
/*                instance of Registration_Packet class)   */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

Registration_Packet::~Registration_Packet() {
  // nothing to do here because Packet already deletes the
  //  data, since nothing else associated with this type of packet.
}

/*********************************************************/
/*  Name : Registration_Packet                           */
/*                                                       */
/*  Description: Default Constructor                     */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

Registration_Packet::Registration_Packet() {
  // This is the actual type for the registration packet
  // It has the type 0
  m_type = 0;

  m_wait_response = TRUE;

  // Create new data packet
  m_data = new Registration_Data();
  
  // Always registered
  m_registered = TRUE;
}

/*********************************************************/
/*  Name : Registration_Packet                           */
/*                                                       */
/*  Description: Constructor with registration data as   */
/*                input                                  */
/*                                                       */
/*  Input : Pointer to registration data to initialize   */
/*           m_data.                                     */
/*  Return Value : None                                  */
/*********************************************************/

Registration_Packet::Registration_Packet(const Registration_Data& reg_data) {
  m_data = new Registration_Data(reg_data);

  // It has the bit type 0
  m_type = 0;
  m_wait_response = TRUE;
  m_registered = TRUE;
}

/*********************************************************/
/*  Name : GetCopyOfData                                 */
/*                                                       */
/*  Description: Returns a copy of m_data through CopyDat*/
/*                (CopyDat must already be created and of*/
/*                 correct type)                         */
/*                                                       */
/*  Input : Pointer to registration data to return copy  */
/*           of m_data.                                  */
/*  Return Value : 1 - successful, -1 - failure          */
/*********************************************************/

int Registration_Packet::GetCopyOfData(Packet_Data * CopyDat) {
  Registration_Data* CopyRegDat = dynamic_cast<Registration_Data*>(CopyDat);
  if (CopyRegDat != NULL) {
    *(CopyRegDat) = *(REG_DATA);
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
/*  Input : Pointer to registration data which contains  */
/*           the values to set m_data (does not set ptr  */
/*           since that would allow the data storage to  */
/*           exist outside of this class and thus        */
/*           vunerable to being deleted)                 */
/*  Return Value : 1 - successful, -1 - failure          */
/*********************************************************/

int Registration_Packet::SetData(Packet_Data * NewDat) {
  Registration_Data* NewRegDat = dynamic_cast<Registration_Data*>(NewDat);
  if (NewRegDat != NULL) {
    *(REG_DATA) = *(NewRegDat);
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

int Registration_Packet::Fill(char * mb) {
  // Check that m_data is valid first
  if (m_data != NULL) {
    // Pointer to current place in buffer
    char * place_ptr = mb;
    int ptr_add = 0;
    
    // Order of registration packet is
    //  Tag returned by server to notify of success (char - 1 byte)
    //  RegClass Class Type (int - 4 bytes)
    //  RegClass Class Size (short - 2 bytes)
    //  RegClass Class Name (can be any size)

    ptr_add  = FillChar(place_ptr, (char)REG_DATA->GetRegTag());
    ptr_add += FillInt(place_ptr + ptr_add, REG_DATA->GetRegClassType());
    ptr_add += FillShort(place_ptr + ptr_add, REG_DATA->GetRegClassSize());
    ptr_add += FillString(place_ptr + ptr_add, REG_DATA->GetRegClassName());
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

int Registration_Packet::ReadPacket(char * mb) {
  // Update data with new information
  if (m_data != NULL) {
    // Pointer to current place in buffer
    char * place_ptr = mb;

    // Order of registration packet is
    //  Success Tag (char - 1 byte)
    //  RegClass Class Type (int - 4 bytes)
    //  RegClass Class Size (short - 2 bytes)
    //  RegClass Class Name (can be any size)

    // Grab and set success tag
    REG_DATA->SetRegTag((TRRegTag)*place_ptr);
    place_ptr += sizeof(char);

    // Grab and set registration class type
    REG_DATA->SetRegClassType(GrabInt(place_ptr));
    place_ptr += sizeof(int);

    // Grab and set registration class size
    REG_DATA->SetRegClassSize(GrabShort(place_ptr));
    place_ptr += sizeof(short);

    // Grab Packet Class Name (which has null terminator)
    char* name;
    name = strdup(place_ptr);

    // Set data to correspond to new values
    if (REG_DATA->SetRegClassName(name)) {
      // clean up memory 
      free(name);
      return 1;
    }
    // clean up memory
    free(name);
  }
  return -1;
}

