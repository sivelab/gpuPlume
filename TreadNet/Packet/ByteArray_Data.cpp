#include "ByteArray_Data.h"

/***********************************************************/
/*  Name : ~ByteArray_Data                                 */
/*                                                         */
/*  Description: Destructor (free memory specific to each  */
/*                instance of ByteArray_Data class)        */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

ByteArray_Data::~ByteArray_Data() {
  // none to free
}

/*********************************************************/
/*  Name : ByteArray_Data                                */
/*                                                       */
/*  Description: Default Constructor                     */
/*                                                       */
/*  Input : None                                         */
/*  Return Value : None                                  */
/*********************************************************/

ByteArray_Data::ByteArray_Data() {
  // Set to default values
  for (int i=0; i < ARRAY_SIZE; i++) {
    m_array[i] = 0;
  }

  // default size is just the size of 0 array + 2 for sending size
  m_size = 2;
  m_LastElementInd = -1;
}

/*********************************************************/
/*  Name : ByteArray_Data                                */
/*                                                       */
/*  Description: Copy Constructor                        */
/*                                                       */
/*  Input : ByteArray_Data&                              */
/*  Return Value : None                                  */
/*********************************************************/

ByteArray_Data::ByteArray_Data(const ByteArray_Data& right) 
	: Packet_Data(right) {
  const char* rptr = right.GetArrayPtr();
  for (int i=0; i < ARRAY_SIZE; i++) {
    m_array[i] = rptr[i];
  }

  m_size = right.GetSize();
  m_LastElementInd = right.GetLastElementInd();
}

/*********************************************************/
/*  Name : Operator =                                    */
/*                                                       */
/*  Description: Equals operator                         */
/*                                                       */
/*  Input : ByteArray_Data&                              */
/*  Return Value : None                                  */
/*********************************************************/

void ByteArray_Data::operator = (const ByteArray_Data& right) {
  const char* rptr = right.GetArrayPtr();
  for (int i=0; i < ARRAY_SIZE; i++) {
    m_array[i] = rptr[i];
  }

  m_size = right.GetSize();
  m_LastElementInd = right.GetLastElementInd();
}

/*********************************************************/
/*  Name : SetByte                                       */
/*                                                       */
/*  Description: Sets a value in the array               */
/*                                                       */
/*  Input : int, char                                    */
/*  Return Value : 0 if failed, 1 otherwise              */
/*********************************************************/

int ByteArray_Data::SetByte(int index, char value) {
  if (index >= 0 && index < ARRAY_SIZE) {
    m_array[index] = value;
    return 1;
  }
  return 0;
}

/*********************************************************/
/*  Name : GetArrayPtr (const and non-const version)     */
/*                                                       */
/*  Description: Returns ptr to array at index           */
/*                                                       */
/*  Input : int index                                    */
/*  Return Value : Ptr (0 if bad index)                  */
/*********************************************************/

char * ByteArray_Data::GetArrayPtr(int index) {
  if (index >= 0 && index < ARRAY_SIZE) {
    return &m_array[index];
  }
  return 0;
}

const char * ByteArray_Data::GetArrayPtr(int index) const {
  if (index >= 0 && index < ARRAY_SIZE) {
    return &m_array[index];
  }
  return 0;
}

/*********************************************************/
/*  Name : SetLastElementInd                             */
/*                                                       */
/*  Description: Sets the index of the last element      */
/*                                                       */
/*  Input : int                                          */
/*  Return Value : 0 if failed, 1 otherwise              */
/*********************************************************/

int ByteArray_Data::SetLastElementInd(int val) {
  if (val < 0) {
    // too small of index
    m_LastElementInd = -1;
    m_size = 0;
    return 0;
  }
  else if (val >= ARRAY_SIZE) {
    // too large, but set m_size to max
    m_LastElementInd = ARRAY_SIZE - 1;
    m_size = ARRAY_SIZE;
    return 0;
  }
  // good
  m_LastElementInd = val;
  // array size is val + 1, plus 2 bytes for sending the array size
  m_size = val + 3;
  return 1;
}
