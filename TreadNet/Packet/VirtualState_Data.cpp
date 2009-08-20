#include "VirtualState_Data.h"

/***********************************************************/
/*  Name : ~VirtualState_Data                              */
/*                                                         */
/*  Description: Destructor (free memory specific to each  */
/*                instance of VirtualState_Data class)     */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

VirtualState_Data::~VirtualState_Data() {
  // nothing to do since no "new" called on construction
}

/*********************************************************/
/*  Name : VirtualState_Data                             */
/*                                                       */
/*  Description: Default Constructor                     */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

VirtualState_Data::VirtualState_Data() {
  // Set to default values
  m_surface_normal.X(0);
  m_surface_normal.Y(0);
  m_surface_normal.Z(1);

  m_terrain = NORMAL_TERRAIN;
  
  m_virtual_force.X(0);
  m_virtual_force.Y(0);
  m_virtual_force.Z(0);

  m_movement_scale = 1;
  m_user_scale = 1;

  // Sending terrain as char to save space
  m_size = m_surface_normal.GetSize() + sizeof(char) +
    m_virtual_force.GetSize() + sizeof(m_movement_scale) +
    sizeof(m_user_scale);
}

/*********************************************************/
/*  Name : VirtualState_Data                             */
/*                                                       */
/*  Description: Copy Constructor                        */
/*                                                       */
/*  Input : VirtualState_Data&                           */
/*  Return Value : None                                  */
/*********************************************************/

VirtualState_Data::VirtualState_Data(const VirtualState_Data& right) 
	: Packet_Data(right) {
  m_surface_normal = right.GetSurfaceNormal();
  m_terrain = right.GetTerrain();
  m_virtual_force = right.GetVirtualForce();
  m_movement_scale = right.GetMovementScale();
  m_user_scale = right.GetUserScale();

  // Sending terrain as char to save space
  m_size = m_surface_normal.GetSize() + sizeof(char) +
    m_virtual_force.GetSize() + sizeof(m_movement_scale) +
    sizeof(m_user_scale);
}

/*********************************************************/
/*  Name : Operator =                                    */
/*                                                       */
/*  Description: Equals operator                         */
/*                                                       */
/*  Input : VirtualState_Data&                           */
/*  Return Value : None                                  */
/*********************************************************/

void VirtualState_Data::operator = (const VirtualState_Data& right) {
  m_surface_normal = right.GetSurfaceNormal();
  m_terrain = right.GetTerrain();
  m_virtual_force = right.GetVirtualForce();
  m_movement_scale = right.GetMovementScale();
  m_user_scale = right.GetUserScale();
}
