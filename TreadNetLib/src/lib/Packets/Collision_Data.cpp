#include "Collision_Data.h"

/***********************************************************/
/*  Name : ~Collision_Data                                 */
/*                                                         */
/*  Description: Destructor (free memory specific to each  */
/*                instance of Collision_Data class)        */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

Collision_Data::~Collision_Data() {
  // nothing to do since no "new" called on construction
}

/*********************************************************/
/*  Name : Collision_Data                                */
/*                                                       */
/*  Description: Default Constructor                     */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

Collision_Data::Collision_Data() {
  // Set to default values
  m_collision_normal.X(0);
  m_collision_normal.Y(0);
  m_collision_normal.Z(0);

  m_collision_pos.X(0);
  m_collision_pos.Y(0);
  m_collision_pos.Z(0);

  m_collision_flag = NO_COLLISION;

  // sending flag as char
  m_size = m_collision_normal.GetSize() + sizeof(char) +
    m_collision_pos.GetSize();
}

/*********************************************************/
/*  Name : Collision_Data                                */
/*                                                       */
/*  Description: Copy Constructor                        */
/*                                                       */
/*  Input : Collision_Data&                              */
/*  Return Value : None                                  */
/*********************************************************/

Collision_Data::Collision_Data(const Collision_Data& right) 
	: Packet_Data(right) {
  m_collision_normal = right.GetCollisionNormal();
  m_collision_pos = right.GetCollisionPos();
  m_collision_flag = right.GetCollisionFlag();

  m_size = m_collision_normal.GetSize() + sizeof(char) +
    m_collision_pos.GetSize();
}

/*********************************************************/
/*  Name : Operator =                                    */
/*                                                       */
/*  Description: Equals operator                         */
/*                                                       */
/*  Input : Collision_Data&                              */
/*  Return Value : None                                  */
/*********************************************************/

void Collision_Data::operator = (const Collision_Data& right) {
  m_collision_normal = right.GetCollisionNormal();
  m_collision_pos = right.GetCollisionPos();
  m_collision_flag = right.GetCollisionFlag();
}
