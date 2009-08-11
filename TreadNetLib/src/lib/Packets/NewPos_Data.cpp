#include "NewPos_Data.h"

/***********************************************************/
/*  Name : ~NewPos_Data                                    */
/*                                                         */
/*  Description: Destructor (free memory specific to each  */
/*                instance of NewPos_Data class)           */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

NewPos_Data::~NewPos_Data() {
  // nothing to do since no "new" called on construction
}

/*********************************************************/
/*  Name : NewPos_Data                                   */
/*                                                       */
/*  Description: Default Constructor                     */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

NewPos_Data::NewPos_Data() {
  // When the user state and virtual state are created
  //  they will automatically set default values

  // set size based on what the other data give
  m_size = m_user_state.GetSize() + m_virtual_state.GetSize();
}

/*********************************************************/
/*  Name : NewPos_Data                                   */
/*                                                       */
/*  Description: Copy Constructor                        */
/*                                                       */
/*  Input : NewPos_Data&                                 */
/*  Return Value : None                                  */
/*********************************************************/

NewPos_Data::NewPos_Data(const NewPos_Data& right) 
	: Packet_Data(right) {
  m_user_state = right.GetUserState();
  m_virtual_state = right.GetVirtualState();

  // set size based on what the other data give
  m_size = m_user_state.GetSize() + m_virtual_state.GetSize();
}

/*********************************************************/
/*  Name : NewPos_Data                                   */
/*                                                       */
/*  Description: Constructor w/ indiv. data sent in      */
/*                                                       */
/*  Input : UserState_Data& and VirtualState_Data&       */
/*  Return Value : None                                  */
/*********************************************************/

NewPos_Data::NewPos_Data(const UserState_Data& us, const VirtualState_Data& vs) {
  m_user_state = us;
  m_virtual_state = vs;

  // set size based on what the other data give
  m_size = m_user_state.GetSize() + m_virtual_state.GetSize();
}

/*********************************************************/
/*  Name : Operator =                                    */
/*                                                       */
/*  Description: Equals operator                         */
/*                                                       */
/*  Input : NewPos_Data&                                 */
/*  Return Value : None                                  */
/*********************************************************/

void NewPos_Data::operator = (const NewPos_Data& right) {
  m_user_state = right.GetUserState();
  m_virtual_state = right.GetVirtualState();
}
