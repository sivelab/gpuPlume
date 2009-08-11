#include "Command_Data.h"

/***********************************************************/
/*  Name : ~Command_Data                                    */
/*                                                         */
/*  Description: Destructor (free memory specific to each  */
/*                instance of Command_Data class)           */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

Command_Data::~Command_Data() {
}

/*********************************************************/
/*  Name : Command_Data                                   */
/*                                                       */
/*  Description: Default Constructor                     */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

Command_Data::Command_Data() {
  // Set to default values
  m_command = NO_CMD;

  // sending command as a char
  m_size = sizeof(char);
}

/*********************************************************/
/*  Name : Command_Data                                   */
/*                                                       */
/*  Description: Copy Constructor                        */
/*                                                       */
/*  Input : Command_Data&                                 */
/*  Return Value : None                                  */
/*********************************************************/

Command_Data::Command_Data(const Command_Data& right) 
	: Packet_Data(right) {
  m_command = right.GetCmd();

  m_size = sizeof(char);
}

/*********************************************************/
/*  Name : Operator =                                    */
/*                                                       */
/*  Description: Equals operator                         */
/*                                                       */
/*  Input : Command_Data&                                 */
/*  Return Value : None                                  */
/*********************************************************/

void Command_Data::operator = (const Command_Data& right) {
  m_command = right.GetCmd();
}
