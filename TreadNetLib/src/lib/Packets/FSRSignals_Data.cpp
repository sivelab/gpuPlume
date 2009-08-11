#include "FSRSignals_Data.h"

/***********************************************************/
/*  Name : ~FSRSignals_Data                                */
/*                                                         */
/*  Description: Destructor (free memory specific to each  */
/*                instance of FSRSignals_Data class)       */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

FSRSignals_Data::~FSRSignals_Data() {
  // nothing to do since no "new" called on construction
}

/*********************************************************/
/*  Name : FSRSignals_Data                               */
/*                                                       */
/*  Description: Default Constructor                     */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

FSRSignals_Data::FSRSignals_Data() {
  // Set to default values
  m_right_heel = 0.0;
  m_right_toe = 0.0;
  m_left_heel = 0.0;
  m_left_toe = 0.0;
 
  m_size = 4*sizeof(float);
}

/*********************************************************/
/*  Name : FSRSignals_Data                               */
/*                                                       */
/*  Description: Copy Constructor                        */
/*                                                       */
/*  Input : FSRSignals_Data&                             */
/*  Return Value : None                                  */
/*********************************************************/

FSRSignals_Data::FSRSignals_Data(const FSRSignals_Data& right) 
	: Packet_Data(right) {
  m_right_heel = right.GetValueRightHeel();
  m_right_toe = right.GetValueRightToe();
  m_left_heel = right.GetValueLeftHeel();
  m_left_toe = right.GetValueLeftToe();

  m_size = 4*sizeof(float);
}

/*********************************************************/
/*  Name : Operator =                                    */
/*                                                       */
/*  Description: Equals operator                         */
/*                                                       */
/*  Input : UserState_Data&                              */
/*  Return Value : None                                  */
/*********************************************************/

void FSRSignals_Data::operator = (const FSRSignals_Data& right) {

  m_right_heel = right.GetValueRightHeel();
  m_right_toe = right.GetValueRightToe();
  m_left_heel = right.GetValueLeftHeel();
  m_left_toe = right.GetValueLeftToe();
  
}
