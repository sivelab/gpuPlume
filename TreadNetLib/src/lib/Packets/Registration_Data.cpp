#include "Registration_Data.h"

/***********************************************************/
/*  Name : ~Registration_Data                              */
/*                                                         */
/*  Description: Destructor (free memory specific to each  */
/*                instance of Registration_Data class)     */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

Registration_Data::~Registration_Data() {
  free(m_regclass_name);
}

/*********************************************************/
/*  Name : Registration_Data                             */
/*                                                       */
/*  Description: Default Constructor                     */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

Registration_Data::Registration_Data() {
  // Set to default values
  m_regclass_name = NULL;
  m_regclass_data_size = 0;
  m_regclass_type = 0;
  m_tag = UNKNOWN_CLASS;

  // set size to zero so know not initialized
  m_size = 0;
}

/*********************************************************/
/*  Name : Registration_Data                             */
/*                                                       */
/*  Description: Copy Constructor                        */
/*                                                       */
/*  Input : Registration_Data&                           */
/*  Return Value : None                                  */
/*********************************************************/

Registration_Data::Registration_Data(const Registration_Data& right) 
	: Packet_Data(right) {
  m_regclass_name = strdup(right.GetRegClassName());
  m_regclass_data_size = right.GetRegClassSize();
  m_regclass_type = right.GetRegClassType();
  m_tag = right.GetRegTag();

  CalcDataSize();

}

/*********************************************************/
/*  Name : Operator =                                    */
/*                                                       */
/*  Description: Equals operator                         */
/*                                                       */
/*  Input : Registration_Data&                           */
/*  Return Value : None                                  */
/*********************************************************/

void Registration_Data::operator = (const Registration_Data& right) {
  if (m_regclass_name != NULL) {
    free(m_regclass_name);
  }
  m_regclass_name = strdup(right.GetRegClassName());
  m_regclass_data_size = right.GetRegClassSize();
  m_regclass_type = right.GetRegClassType();
  m_tag = right.GetRegTag();

  CalcDataSize();
}

/*********************************************************/
/*  Name : SetRegClassName                               */
/*                                                       */
/*  Description: Sets packet name to register            */
/*                                                       */
/*  Input : char *                                       */
/*  Return Value : 0 if failed, 1 otherwise              */
/*********************************************************/

int Registration_Data::SetRegClassName(const char *c_name) {
  if (m_regclass_name != NULL) {
    free(m_regclass_name);
  }
  m_regclass_name = strdup(c_name);

  // Need to reset the size for the packet since length
  //  of name could have changed.
  CalcDataSize();

  return m_regclass_name != NULL;
}

/*********************************************************/
/*  Name : CalcDataSize                                  */
/*                                                       */
/*  Description: Recalculates how many bytes this data   */
/*                packet will take up and sets m_size.   */
/*               Internal routine to centralize size     */
/*                calculation to reduce chance of error. */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None (sets m_size)                    */
/*********************************************************/

void Registration_Data::CalcDataSize() {
  // Size equals length of packet_name (plus null character) and
  //  the 3 other packet info variables - with tag being a char
  if (m_regclass_name != NULL) {
    m_size = strlen(m_regclass_name) + sizeof(m_regclass_data_size) +
      sizeof(m_regclass_type) + sizeof(char) + 1;
  }
  else {
    m_size = sizeof(m_regclass_data_size) + sizeof(m_regclass_type) +
      sizeof(char);
  }
}
