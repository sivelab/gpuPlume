#include "FSRSignals_Packet.h"

char * FSRSignals_Packet::m_class_name = strdup("FSRSignals_Packet");

/***********************************************************/
/*  Name : ~FSRSignals_Packet                              */
/*                                                         */
/*  Description: Destructor (free memory specific to each  */
/*                instance of FSRSignals_Packet class)     */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

FSRSignals_Packet::~FSRSignals_Packet() {
  // nothing to do here because Packet already deletes the
  //  data, since nothing else associated with this type of packet.
}

/*********************************************************/
/*  Name : FSRSignals_Packet                             */
/*                                                       */
/*  Description: Default Constructor                     */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

FSRSignals_Packet::FSRSignals_Packet() {
  // Some default values
  m_registered = FALSE;
  m_type = 0;
  m_wait_response = FALSE;

  // Create new data packet
  m_data = new FSRSignals_Data();
}

/*********************************************************/
/*  Name : FSRSignals_Packet                             */
/*                                                       */
/*  Description: Constructor with FSRSignal data as input*/
/*                                                       */
/*  Input : Reference to FSR Signal data to initialize   */
/*           m_data                                      */
/*  Return Value : None                                  */
/*********************************************************/

FSRSignals_Packet::FSRSignals_Packet(const FSRSignals_Data& fsr_data) {
  m_data = new FSRSignals_Data(fsr_data);

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
/*  Input : Pointer to FSRSignals data to return copy    */
/*           of m_data.                                  */
/*  Return Value : 1 - successful, -1 - failure          */
/*********************************************************/

int FSRSignals_Packet::GetCopyOfData(Packet_Data * CopyDat) {
  FSRSignals_Data* CopyFSRDat = dynamic_cast<FSRSignals_Data*>(CopyDat);
  if (CopyFSRDat != NULL) {
    *(CopyFSRDat) = *(FSRSIGNALS_DATA);
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
/*  Input : Pointer to FSRSignals data which contains the*/
/*           values to set m_data (does not set ptr      */
/*           since that would allow the data storage to  */
/*           exist outside of this class and thus        */
/*           vunerable to being deleted)                 */
/*  Return Value : 1 - successful, -1 - failure          */
/*********************************************************/

int FSRSignals_Packet::SetData(Packet_Data * NewDat) {
  FSRSignals_Data* NewFSRDat = dynamic_cast<FSRSignals_Data*>(NewDat);
  if (NewFSRDat != NULL) {
    *(FSRSIGNALS_DATA) = *(NewFSRDat);
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

int FSRSignals_Packet::Fill(char * mb) {
  // Check that m_data is valid first
  if (m_data != NULL) {
    // Pointer to current place in buffer
    char * place_ptr = mb;
    int ptr_add = 0;
    
    // Order of FSRSignals packet is
    // RightHeel - one float, 4? bytes
    // RightToe - one float
    // LeftHeel - one float
    // LeftToe - one float

    ptr_add  = FillFloat(place_ptr, FSRSIGNALS_DATA->GetValueRightHeel());
    ptr_add += FillFloat(place_ptr + ptr_add, FSRSIGNALS_DATA->GetValueRightToe());
    ptr_add += FillFloat(place_ptr + ptr_add, FSRSIGNALS_DATA->GetValueLeftHeel());
    ptr_add += FillFloat(place_ptr + ptr_add, FSRSIGNALS_DATA->GetValueLeftToe());
  
    return ptr_add;
  }
  return -1;
}

/*********************************************************/
/*  Name : Read Packet                                   */
/*                                                       */
/*  Description: Read data from packet?                  */
/*                                                       */
/*                                                       */
/*  Input : pointer to beginning of portion of message   */
/*           buffer to receive data                      */
/*  Return Value : Returns success (1) or failure (-1)   */
/*********************************************************/

int FSRSignals_Packet::ReadPacket(char * mb) {
  // Update data with new information
  if (m_data != NULL) {
    // Pointer to current place in buffer
    char * place_ptr = mb;

    // Order of FSRSignals packet is
    // RightHeel - one float, 4? bytes
    // RightToe - one float
    // LeftHeel - one float
    // LeftToe - one float


    // Grab and set position values
    FSRSIGNALS_DATA->SetValueRightHeel(GrabFloat(place_ptr));
    place_ptr += sizeof(float);
    FSRSIGNALS_DATA->SetValueRightToe(GrabFloat(place_ptr));
    place_ptr += sizeof(float);
    FSRSIGNALS_DATA->SetValueLeftHeel(GrabFloat(place_ptr));
    place_ptr += sizeof(float);
    FSRSIGNALS_DATA->SetValueLeftToe(GrabFloat(place_ptr));
    place_ptr += sizeof(float);

    return 1;
  }
  return -1;
}
