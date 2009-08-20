#include "ByteArray_Packet.h"

char * ByteArray_Packet::m_class_name = strdup("ByteArray_Packet");

/***********************************************************/
/*  Name : ~ByteArray_Packet                               */
/*                                                         */
/*  Description: Destructor (free memory specific to each  */
/*                instance of ByteArray_Packet class)      */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

ByteArray_Packet::~ByteArray_Packet() {
  // nothing to do here because Packet already deletes the
  //  data, since nothing else associated with this type of packet.
}

/*********************************************************/
/*  Name : ByteArray_Packet                              */
/*                                                       */
/*  Description: Default Constructor                     */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

ByteArray_Packet::ByteArray_Packet() {
  // some default values
  m_type = 0;
  m_wait_response = FALSE;
  m_registered = FALSE;

  // Create new data packet
  m_data = new ByteArray_Data();
}

/*********************************************************/
/*  Name : ByteArray_Packet                              */
/*                                                       */
/*  Description: Constructor with bytearray data as      */
/*                input                                  */
/*                                                       */
/*  Input : Pointer to bytearray data to initialize      */
/*           m_data.                                     */
/*  Return Value : None                                  */
/*********************************************************/

ByteArray_Packet::ByteArray_Packet(const ByteArray_Data& ba_data) {
  m_data = new ByteArray_Data(ba_data);

  // default
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
/*  Input : Pointer to bytearray data to return copy     */
/*           of m_data.                                  */
/*  Return Value : 1 - successful, -1 - failure          */
/*********************************************************/

int ByteArray_Packet::GetCopyOfData(Packet_Data * CopyDat) {
  ByteArray_Data* CopyBADat = dynamic_cast<ByteArray_Data*>(CopyDat);
  if (CopyBADat != NULL) {
    *(CopyBADat) = *(BA_DATA);
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
/*  Input : Pointer to bytearray data which contains     */
/*           the values to set m_data (does not set ptr  */
/*           since that would allow the data storage to  */
/*           exist outside of this class and thus        */
/*           vunerable to being deleted)                 */
/*  Return Value : 1 - successful, -1 - failure          */
/*********************************************************/

int ByteArray_Packet::SetData(Packet_Data * NewDat) {
  ByteArray_Data* NewBADat = dynamic_cast<ByteArray_Data*>(NewDat);
  if (NewBADat != NULL) {
    *(BA_DATA) = *(NewBADat);
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

int ByteArray_Packet::Fill(char * mb) {
  // Also check that m_data is valid
  if (m_data != NULL) {
    // Pointer to current place in buffer
    char * place_ptr = mb;
    int ptr_add = 0;
    
    // Order of bytearray packet is
    //  array size (short - 2 bytes)
    //  Byte array (can be any size)

    // get array size
    int size = BA_DATA->GetLastElementInd() + 1;
    ptr_add = FillShort(place_ptr, size);
    ptr_add += FillCharArray(place_ptr + ptr_add, BA_DATA->GetArrayPtr(), size);
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

int ByteArray_Packet::ReadPacket(char * mb) {
  // Update data with new information
  if (m_data != NULL) {
    // Pointer to current place in buffer
    char * place_ptr = mb;

    // Order of bytearray packet is
    //  array size (short - 2 bytes)
    //  Byte array (can be any size)

    // Grab and set array size size
    int array_size = GrabShort(place_ptr);
    BA_DATA->SetLastElementInd(array_size - 1);
    place_ptr += sizeof(short);

    // Grab byte array
    GrabCharArray(place_ptr, BA_DATA->GetArrayPtr(), array_size);

    // set m_read_size correctly
    //    m_read_size = 2 + array_size;
    return 1;
  }
  return -1;
}

