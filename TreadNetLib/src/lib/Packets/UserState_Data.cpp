#include "UserState_Data.h"

/***********************************************************/
/*  Name : ~UserState_Data                                 */
/*                                                         */
/*  Description: Destructor (free memory specific to each  */
/*                instance of UserState_Data class)        */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

UserState_Data::~UserState_Data() {
  // nothing to do since no "new" called on construction
}

/*********************************************************/
/*  Name : UserState_Data                                */
/*                                                       */
/*  Description: Default Constructor                     */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

UserState_Data::UserState_Data() {
  // Set to default values
  m_pos.X(0);
  m_pos.Y(0);
  m_pos.Z(0);
  m_user_offset.X(0);
  m_user_offset.Y(0);
  m_user_offset.Z(0);
  m_eye_offset.X(0);
  m_eye_offset.Y(0);
  m_eye_offset.Z(0);
  m_facing.X(0);
  m_facing.Y(1);
  m_facing.Z(0);
  m_velocity.X(0);
  m_velocity.Y(0);
  m_velocity.Z(0);
  m_status = OFF;
  m_valid = 0;

  // Not sending every value, thus m_size is not exactly
  //  what you would expect (not sending pos.Z and facing.Z)
  m_size = 13*sizeof(double) + sizeof(m_status) + sizeof(m_valid);
}

/*********************************************************/
/*  Name : UserState_Data                                */
/*                                                       */
/*  Description: Copy Constructor                        */
/*                                                       */
/*  Input : UserState_Data&                              */
/*  Return Value : None                                  */
/*********************************************************/

UserState_Data::UserState_Data(const UserState_Data& right) 
	: Packet_Data(right) {
  m_pos = right.GetPos();
  m_user_offset = right.GetUserOffset();
  m_eye_offset = right.GetEyeOffset();
  m_facing = right.GetFacing();
  m_velocity = right.GetVelocity();
  m_status = right.GetStatus();
  m_valid = right.PosValid();

  // Not sending every value, thus m_size is not exactly
  //  what you would expect (not sending pos.Z and facing.Z)
  m_size = 13*sizeof(double) + sizeof(m_status) + sizeof(m_valid);
}

/*********************************************************/
/*  Name : Operator =                                    */
/*                                                       */
/*  Description: Equals operator                         */
/*                                                       */
/*  Input : UserState_Data&                              */
/*  Return Value : None                                  */
/*********************************************************/

void UserState_Data::operator = (const UserState_Data& right) {
  m_pos = right.GetPos();
  m_user_offset = right.GetUserOffset();
  m_eye_offset = right.GetEyeOffset();
  m_facing = right.GetFacing();
  m_velocity = right.GetVelocity();
  m_status = right.GetStatus();
  m_valid = right.PosValid();
}


/*********************************************************/
/*  Name : GetLocalUserOffset                            */
/*                                                       */
/*  Description: Gets User Offset in Treadport local     */
/*   coordinate system                                   */
/*                                                       */
/*  Note: it uses m_facing for the conversion, so don't  */
/*   modify it after setting user_offset or this will    */
/*   return the value incorrectly.                       */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : TRPoint - UserOffset                  */
/*********************************************************/

TRPoint UserState_Data::GetLocalUserOffset() const {
  // user local conversion function
  TRPoint luo = m_user_offset;
  ConvertVWToLocal(luo);
  return luo;
}


/*********************************************************/
/*  Name : GetLocalEyeOffset                             */
/*                                                       */
/*  Description: Gets Eye Offset in Treadport local      */
/*   coordinate system                                   */
/*                                                       */
/*  Note: it uses m_facing for the conversion, so don't  */
/*   modify it after setting user_offset or this will    */
/*   return the value incorrectly.                       */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : TRPoint - EyeOffset                   */
/*********************************************************/

TRPoint UserState_Data::GetLocalEyeOffset() const {
  // user local conversion function
  TRPoint leo = m_eye_offset;
  ConvertVWToLocal(leo);
  return leo;
}


/*********************************************************/
/*  Name : SetLocalEyeOffset                             */
/*                                                       */
/*  Description: Sets Eye Offset in Treadport local      */
/*   coordinate system                                   */
/*                                                       */
/*  Note: it uses m_facing for the conversion, so make   */
/*   sure to set it correctly before calling this        */
/*   function or you'll get the wrong conversion to      */
/*   virtual world coordinates.                          */
/*                                                       */
/*  Input : TRPoint& - EyeOffset                         */
/*  Return Value : None                                  */
/*********************************************************/

void UserState_Data::SetLocalUserOffset(const TRPoint& uoff) {
  // Set user offset, and use local conversion function
  m_user_offset = uoff;
  ConvertLocalToVW(m_user_offset);
}

/*********************************************************/
/*  Name : SetLocalUserOffset                            */
/*                                                       */
/*  Description: Sets User Offset in Treadport local     */
/*   coordinate system                                   */
/*                                                       */
/*  Note: it uses m_facing for the conversion, so make   */
/*   sure to set it correctly before calling this        */
/*   function or you'll get the wrong conversion to      */
/*   virtual world coordinates.                          */
/*                                                       */
/*  Input : TRPoint&                                     */
/*  Return Value : None                                  */
/*********************************************************/

void UserState_Data::SetLocalEyeOffset(const TRPoint& eoff) {
  // Set user offset, and use local conversion function
  m_eye_offset = eoff;
  ConvertLocalToVW(m_eye_offset);
}


/*********************************************************/
/*  Name : ConvertLocalToVW                              */
/*                                                       */
/*  Description: Converts a vector/point from local      */
/*   treadport coordinate system to virtual world system */
/*                                                       */
/*  Input : TRVector&                                    */
/*  Return Value : None                                  */
/*********************************************************/

void UserState_Data::ConvertLocalToVW(TRVector& v) const {
  // Explanation of conversion:
  //  Facing represents the local Treadport X axis in the virtual 
  //  world, the Y axis (for left handed coordinate systems)
  //  is just [-X_axis.y X_axix.x]
  //
  // Thus the rotation matrix to convert from treadport frame
  //  to virtual world frame is (g = graphics, t = treadport)
  //
  // g       | g      g    |    |  Facing.x   -Facing.y |
  //  R    = |  X      Y   | =  |                       |
  //   t     |   t      t  |    |  Facing.y    Facing.x |
  //
  // Also note this is a 2D conversion, because the Z axis is the same
  //  in both coordinate frames
  
  double newx, newy;
  newx = v.X() * m_facing.X() - v.Y() * m_facing.Y();
  newy = v.X() * m_facing.Y() + v.Y() * m_facing.X();
  v.X() = newx;
  v.Y() = newy;
}

/*********************************************************/
/*  Name : ConvertVWToLocal                              */
/*                                                       */
/*  Description: Converts a vector/point from virtual    */
/*   world coordinate system to local treadport system   */
/*                                                       */
/*  Input : TRVector&                                    */
/*  Return Value : None                                  */
/*********************************************************/

void UserState_Data::ConvertVWToLocal(TRVector& v) const {
  // Explanation of conversion:
  // Take Rotation matrix from above and calculate the opposite
  //  conversion (which we want) is
  //
  // t       |  Facing.x   Facing.y |
  //  R    = |                      |
  //   g     | -Facing.y   Facing.x |
  //
  
  double newx, newy;
  newx =  v.X() * m_facing.X() + v.Y() * m_facing.Y();
  newy = -v.X() * m_facing.Y() + v.Y() * m_facing.X();
  v.X() = newx;
  v.Y() = newy;
}
