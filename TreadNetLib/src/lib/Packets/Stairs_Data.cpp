#include "Stairs_Data.h"

/***********************************************************/
/*  Name : ~Stairs_Data                                    */
/*                                                         */
/*  Description: Destructor (free memory specific to each  */
/*                instance of Stairs_Data class)           */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

Stairs_Data::~Stairs_Data() {
}

/*********************************************************/
/*  Name : Stairs_Data                                   */
/*                                                       */
/*  Description: Default Constructor                     */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

Stairs_Data::Stairs_Data() {
  // Set to default values
  m_stairs_dir = IGNORE_STAIRS;
  m_stairs_cmd = NO_STAIRS_CMD;

  // sending both as chars
  m_size = sizeof(char) + sizeof(char);
}

/*********************************************************/
/*  Name : Stairs_Data                                   */
/*                                                       */
/*  Description: Copy Constructor                        */
/*                                                       */
/*  Input : Stairs_Data&                                 */
/*  Return Value : None                                  */
/*********************************************************/

Stairs_Data::Stairs_Data(const Stairs_Data& right) 
	: Packet_Data(right) {
  m_stairs_dir = right.GetStairsDir();
  m_stairs_cmd = right.GetCmd();

  m_size = sizeof(char) + sizeof(char);
}

/*********************************************************/
/*  Name : Operator =                                    */
/*                                                       */
/*  Description: Equals operator                         */
/*                                                       */
/*  Input : Stairs_Data&                                 */
/*  Return Value : None                                  */
/*********************************************************/

void Stairs_Data::operator = (const Stairs_Data& right) {
  m_stairs_dir = right.GetStairsDir();
  m_stairs_cmd = right.GetCmd();
}
