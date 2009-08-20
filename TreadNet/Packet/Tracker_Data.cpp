#include "Tracker_Data.h"

/***********************************************************/
/*  Name : ~Tracker_Data                                   */
/*                                                         */
/*  Description: Destructor (free memory specific to each  */
/*                instance of Tracker_Data class)          */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

Tracker_Data::~Tracker_Data() {
  // nothing to do since no "new" called on construction
}

/*********************************************************/
/*  Name : Tracker_Data                                  */
/*                                                       */
/*  Description: Default Constructor                     */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

Tracker_Data::Tracker_Data() {
  // Set to default values (origin, no rotation)
  TRVector vzero(0);

  m_pos = vzero;
  m_orient.qV() = vzero;
  m_orient.qW() = 1;

  m_size = m_pos.GetSize() + m_orient.GetSize();
}

/*********************************************************/
/*  Name : Tracker_Data                                  */
/*                                                       */
/*  Description: Copy Constructor                        */
/*                                                       */
/*  Input : Tracker_Data&                                */
/*  Return Value : None                                  */
/*********************************************************/

Tracker_Data::Tracker_Data(const Tracker_Data& right) 
	: Packet_Data(right) {
  m_pos = right.GetPos();
  m_orient = right.GetOrientation();

  m_size = m_pos.GetSize() + m_orient.GetSize();
}

/*********************************************************/
/*  Name : Operator =                                    */
/*                                                       */
/*  Description: Equals operator                         */
/*                                                       */
/*  Input : Tracker_Data&                                */
/*  Return Value : None                                  */
/*********************************************************/

void Tracker_Data::operator = (const Tracker_Data& right) {
  m_pos = right.GetPos();
  m_orient = right.GetOrientation();
}
