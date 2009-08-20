#include "Tracker_Packet.h"

char * Tracker_Packet::m_class_name = strdup("Tracker_Packet");

/***********************************************************/
/*  Name : ~Tracker_Packet                                 */
/*                                                         */
/*  Description: Destructor (free memory specific to each  */
/*                instance of Tracker_Packet class)        */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

Tracker_Packet::~Tracker_Packet() {
  // nothing to do here because Packet already deletes the
  //  data, and since nothing else associated with this type of packet.
}

/*********************************************************/
/*  Name : Tracker_Packet                                */
/*                                                       */
/*  Description: Default Constructor                     */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

Tracker_Packet::Tracker_Packet() {
  // Some default values
  m_registered = FALSE;
  m_type = 0;
  m_wait_response = FALSE;

  // Create new data packet
  m_data = new Tracker_Data();
}

/*********************************************************/
/*  Name : Tracker_Packet                                */
/*                                                       */
/*  Description: Constructor with Tracker data as input  */
/*                                                       */
/*  Input : Reference to Tracker data to initialize      */
/*           m_data                                      */
/*  Return Value : None                                  */
/*********************************************************/

Tracker_Packet::Tracker_Packet(const Tracker_Data& t_data) {
  m_data = new Tracker_Data(t_data);

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
/*  Input : Pointer to Tracker data to return copy       */
/*           of m_data.                                  */
/*  Return Value : 1 - successful, -1 - failure          */
/*********************************************************/

int Tracker_Packet::GetCopyOfData(Packet_Data * CopyDat) {
  Tracker_Data* CopyTDat = dynamic_cast<Tracker_Data*>(CopyDat);
  if (CopyTDat != NULL) {
    *(CopyTDat) = *(TRACKER_DATA);
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
/*  Input : Pointer to Tracker data which contains the   */
/*           values to set m_data (does not set ptr      */
/*           since that would allow the data storage to  */
/*           exist outside of this class and thus        */
/*           vunerable to being deleted)                 */
/*  Return Value : 1 - successful, -1 - failure          */
/*********************************************************/

int Tracker_Packet::SetData(Packet_Data * NewDat) {
  Tracker_Data* NewTDat = dynamic_cast<Tracker_Data*>(NewDat);
  if (NewTDat != NULL) {
    *(TRACKER_DATA) = *(NewTDat);
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

int Tracker_Packet::Fill(char * mb) {
  // Check that m_data is valid first
  if (m_data != NULL) {
    // Pointer to current place in buffer
    char * place_ptr = mb;
    int ptr_add = 0;
    
    // Order of Tracker packet is
    // Position (X, Y, Z) - 3 doubles, 12 bytes
    // Orientation (Vector, double) - 4 doubles, 16 bytes

    ptr_add  = FillTRVector(place_ptr, TRACKER_DATA->GetPos());
    ptr_add += FillTRQuaternion(place_ptr + ptr_add, TRACKER_DATA->GetOrientation());
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

int Tracker_Packet::ReadPacket(char * mb) {
  // Update data with new information
  if (m_data != NULL) {
    // Pointer to current place in buffer
    char * place_ptr = mb;

    // Order of Tracker packet is
    // Position (X, Y, Z) - 3 doubles, 12 bytes
    // Orientation (Vector, double) - 4 doubles, 16 bytes

    // Grab and set position values
    TRVector temp;

    GrabTRVector(place_ptr, temp);
    TRACKER_DATA->SetPos(temp);
    place_ptr += temp.GetSize();

    TRQuaternion qtemp;
    GrabTRQuaternion(place_ptr, qtemp);
    TRACKER_DATA->SetOrientation(qtemp);
    
    return 1;
  }
  return -1;
}

