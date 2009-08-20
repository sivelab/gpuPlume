#include "Command_Packet.h"

char * Command_Packet::m_class_name = strdup("Command_Packet");

/***********************************************************/
/*  Name : ~Command_Packet                                  */
/*                                                         */
/*  Description: Destructor (free memory specific to each  */
/*                instance of Command_Packet class)         */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

Command_Packet::~Command_Packet() {
  // nothing to do here because Packet already deletes the
  //  data, since nothing else associated with this type of packet.
}

/*********************************************************/
/*  Name : Command_Packet                                 */
/*                                                       */
/*  Description: Default Constructor                     */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

Command_Packet::Command_Packet() {
  // Some default values
  m_registered = FALSE;
  m_type = 0;
  m_wait_response = FALSE;

  // Create new data packet
  m_data = new Command_Data();
}

/*********************************************************/
/*  Name : Command_Packet                                 */
/*                                                       */
/*  Description: Constructor with stairs data as input   */
/*                                                       */
/*  Input : Reference to stairs data to initialize m_data*/
/*  Return Value : None                                  */
/*********************************************************/

Command_Packet::Command_Packet(const Command_Data& cmd_data) {
  m_data = new Command_Data(cmd_data);

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
/*  Input : Pointer to command data to return copy       */
/*           of m_data.                                  */
/*  Return Value : 1 - successful, -1 - failure          */
/*********************************************************/

int Command_Packet::GetCopyOfData(Packet_Data * CopyDat) {
  Command_Data* CopyCmdDat = dynamic_cast<Command_Data*>(CopyDat);
  if (CopyCmdDat != NULL) {
    *(CopyCmdDat) = *(CMD_DATA);
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
/*  Input : Pointer to cmds data which contains the    */
/*           values to set m_data (does not set ptr      */
/*           since that would allow the data storage to  */
/*           exist outside of this class and thus        */
/*           vunerable to being deleted)                 */
/*  Return Value : 1 - successful, -1 - failure          */
/*********************************************************/

int Command_Packet::SetData(Packet_Data * NewDat) {
  Command_Data* NewCmdDat = dynamic_cast<Command_Data*>(NewDat);
  if (NewCmdDat != NULL) {
    *(CMD_DATA) = *(NewCmdDat);
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

int Command_Packet::Fill(char * mb) {
  // Check that m_data is valid first
  if (m_data != NULL) {
    // Order of command packet is
    // Cmd (char - 1 byte)

    return FillChar(mb, (char)CMD_DATA->GetCmd());
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

int Command_Packet::ReadPacket(char * mb) {
  // Update data with new information
  if (m_data != NULL) {
    // Order of command packet is
    //  Cmd (char - 1 byte)

    // Set data to correspond to new values
    CMD_DATA->SetCmd((TRCommand)*mb);
    return 1;
  }
  return -1;
}

