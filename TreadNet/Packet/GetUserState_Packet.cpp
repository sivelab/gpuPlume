#include "GetUserState_Packet.h"

char * GetUserState_Packet::m_class_name = strdup("GetUserState_Packet");

/***********************************************************/
/*  Name : ~GetUserState_Packet                            */
/*                                                         */
/*  Description: Destructor (free memory specific to each  */
/*                instance of GetUserState_Packet class)   */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

GetUserState_Packet::~GetUserState_Packet() {
  // nothing to do here because Packet already deletes the
  //  data, since nothing else associated with this type of packet.
}

/*********************************************************/
/*  Name : GetUserState_Packet                           */
/*                                                       */
/*  Description: Default Constructor                     */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

GetUserState_Packet::GetUserState_Packet() {
  // Some default values
  m_registered = FALSE;
  m_type = 0;
  m_wait_response = TRUE;

  // Create new data packet
  m_data = new GetUserState_Data();
}

/*********************************************************/
/*  Name : GetUserState_Packet                           */
/*                                                       */
/*  Description: Constructor with stairs data as input   */
/*                                                       */
/*  Input : Reference to stairs data to initialize m_data*/
/*  Return Value : None                                  */
/*********************************************************/

GetUserState_Packet::GetUserState_Packet(const GetUserState_Data& gus_data) {
  m_data = new GetUserState_Data();

  // default values
  m_registered = FALSE;
  m_type = 0;
  m_wait_response = TRUE;
}

