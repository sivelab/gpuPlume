#include "Packet.h"

/***********************************************************/
/*  Name : Packet                                          */
/*                                                         */
/*  Description: Default Constructor                       */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

Packet::Packet() {
  // just set some initial values
  m_type = 0;
  m_wait_response = FALSE;
  m_registered = FALSE;
  m_data = NULL;
}

/***********************************************************/
/*  Name : ~Packet                                         */
/*                                                         */
/*  Description: Destructor (free memory specific to each  */
/*                instance of packet class)                */
/*                                                         */
/*  Input : None                                           */
/*  Return Value : None                                    */
/***********************************************************/

Packet::~Packet() {
  delete m_data;
}

/***********************************************************/
/*  Name : Fill Short                                      */
/*                                                         */
/*  Description: Fill message buffer with a short value.   */
/*                Message buffer is in network byte order. */
/*                                                         */
/*  Input : mb - pointer to beginning of portion of        */
/*                message buffer to insert value at.       */
/*          val - short value to insert                    */
/*  Return Value : Returns number of bytes put into buffer */
/***********************************************************/

int Packet::FillShort(char * mb, short val) {
  short nbo_val = htons(val);
  char *ptr_val = (char *) &nbo_val;
  mb[0] = ptr_val[0];
  mb[1] = ptr_val[1];
  return 2;
}

/***********************************************************/
/*  Name : Fill Int(eger)                                  */
/*                                                         */
/*  Description: Fill message buffer with a integer value. */
/*                Message buffer is in network byte order. */
/*                                                         */
/*  Input : mb - pointer to beginning of portion of        */
/*                message buffer to insert value at.       */
/*          val - integer value to insert                  */
/*  Return Value : Returns number of bytes put into buffer */
/***********************************************************/

int Packet::FillInt(char * mb, int val) {
  long nbo_val = htonl(val);
  char *ptr_val = (char *) &nbo_val;
  mb[0] = ptr_val[0];
  mb[1] = ptr_val[1];
  mb[2] = ptr_val[2];
  mb[3] = ptr_val[3];
  return 4;
}

/***********************************************************/
/*  Name : Fill Double                                     */
/*                                                         */
/*  Description: Fill message buffer with a double value.  */
/*                Message buffer is in network byte order. */
/*                                                         */
/*  Input : mb - pointer to beginning of portion of        */
/*                message buffer to insert value at.       */
/*          val - double value to insert                   */
/*  Return Value : Returns number of bytes put into buffer */
/***********************************************************/

int Packet::FillDouble(char * mb, double val) {
  double nbo_val = htond(val);
  char *ptr_val = (char *) &nbo_val;
  mb[0] = ptr_val[0];
  mb[1] = ptr_val[1];
  mb[2] = ptr_val[2];
  mb[3] = ptr_val[3];
  mb[4] = ptr_val[4];
  mb[5] = ptr_val[5];
  mb[6] = ptr_val[6];
  mb[7] = ptr_val[7];
  return 8;
}

/***********************************************************/
/*  Name : Fill Float                                      */
/*                                                         */
/*  Description: Fill message buffer with a float value.   */
/*                Message buffer is in network byte order. */
/*                                                         */
/*  Input : mb - pointer to beginning of portion of        */
/*                message buffer to insert value at.       */
/*          val - float value to insert                    */
/*  Return Value : Returns number of bytes put into buffer */
/***********************************************************/

int Packet::FillFloat(char * mb, float val) {
  float nbo_val = htonf(val);
  char *ptr_val = (char *) &nbo_val;
  mb[0] = ptr_val[0];
  mb[1] = ptr_val[1];
  mb[2] = ptr_val[2];
  mb[3] = ptr_val[3];
  return 4;
}

/***********************************************************/
/*  Name : Fill Char                                       */
/*                                                         */
/*  Description: Fill message buffer with a character.     */
/*                                                         */
/*  Input : mb - pointer to beginning of portion of        */
/*                message buffer to insert value at.       */
/*          val - char value to insert                     */
/*  Return Value : Returns number of bytes put into buffer */
/***********************************************************/

int Packet::FillChar(char * mb, char val) {
  *mb = val;
  return 1;
}

/***********************************************************/
/*  Name : Fill String                                     */
/*                                                         */
/*  Description: Fill message buffer with a string value,  */
/*                including the null character.            */
/*               No need to reorder since just a sequence  */
/*                of bytes.                                */
/*                                                         */
/*  Input : mb - pointer to beginning of portion of        */
/*                message buffer to insert value at.       */
/*          val - string value to insert                   */
/*  Return Value : Returns number of bytes put into buffer */
/***********************************************************/

int Packet::FillString(char * mb, const char* val) {
  if (val == NULL) {
    return 0;
  }
  int num_char = strlen(val) + 1;
  int i;
  for (i=0; i < num_char; i++) {
    mb[i] = val[i];
  }
  return num_char;
}

/***********************************************************/
/*  Name : Fill CharArray                                  */
/*                                                         */
/*  Description: Fill message buffer with a character array*/
/*               No need to reorder since just a sequence  */
/*                of bytes.                                */
/*                                                         */
/*  Input : mb - pointer to beginning of portion of        */
/*                message buffer to insert value at.       */
/*          val - string value to insert (can include      */
/*                 nulls, so can't use strlen)             */
/*          size - size of array to insert                 */
/*  Return Value : Returns number of bytes put into buffer */
/***********************************************************/

int Packet::FillCharArray(char * mb, const char* val, int size) {
  if (val == NULL) {
    return 0;
  }
  int i;
  for (i=0; i < size; i++) {
    mb[i] = val[i];
  }
  return size;
}

/***********************************************************/
/*  Name : Fill TRVector                                   */
/*                                                         */
/*  Description: Fill message buffer with a TRVector       */
/*                                                         */
/*  Input : mb - pointer to beginning of portion of        */
/*                message buffer to insert value at.       */
/*          val - TRVector to insert                       */
/*  Return Value : Returns number of bytes put into buffer */
/***********************************************************/

int Packet::FillTRVector(char * mb, const TRVector& val) {
  // order X, Y, Z
  int num_bytes = FillDouble(mb, val.X());
  num_bytes += FillDouble(mb + sizeof(double), val.Y());
  num_bytes += FillDouble(mb + 2*sizeof(double), val.Z());
  return num_bytes;
}

/***********************************************************/
/*  Name : Fill TRPoint                                    */
/*                                                         */
/*  Description: Fill message buffer with a TRPoint        */
/*                                                         */
/*  Input : mb - pointer to beginning of portion of        */
/*                message buffer to insert value at.       */
/*          val - TRPoint to insert                        */
/*  Return Value : Returns number of bytes put into buffer */
/***********************************************************/

int Packet::FillTRPoint(char * mb, const TRPoint& val) {
  // order X, Y, Z
  int num_bytes = FillDouble(mb, val.X());
  num_bytes += FillDouble(mb + sizeof(double), val.Y());
  num_bytes += FillDouble(mb + 2*sizeof(double), val.Z());
  return num_bytes;
}

/***********************************************************/
/*  Name : Fill TRQuaternion                               */
/*                                                         */
/*  Description: Fill message buffer with a TRQuaternion   */
/*                                                         */
/*  Input : mb - pointer to beginning of portion of        */
/*                message buffer to insert value at.       */
/*          val - TRQuaternion to insert                   */
/*  Return Value : Returns number of bytes put into buffer */
/***********************************************************/

int Packet::FillTRQuaternion(char * mb, const TRQuaternion& val) {
  // order Vector (x,y,z) then double
  int num_bytes = FillTRVector(mb, val.qV());
  num_bytes    += FillDouble(mb + num_bytes, val.qW());
  return num_bytes;
}

/***********************************************************/
/*  Name : Grab Short                                      */
/*                                                         */
/*  Description: Grab a short value from the message       */
/*                buffer, which is in network byte order.  */
/*               Must convert back to host byte order.     */
/*                                                         */
/*  Input : pointer to beginning of portion of message     */
/*           buffer to grab value from.                    */
/*  Return Value : Returns short value                     */
/***********************************************************/

short Packet::GrabShort(char * mb) {
  short* val_ptr = (short *)mb;
  return ntohs(*val_ptr);
}

/***********************************************************/
/*  Name : Grab Int(eger)                                  */
/*                                                         */
/*  Description: Grab an integer value from the message    */
/*                buffer, which is in network byte order.  */
/*               Must convert back to host byte order.     */
/*                                                         */
/*  Input : pointer to beginning of portion of message     */
/*           buffer to grab value from.                    */
/*  Return Value : Returns integer value                   */
/***********************************************************/

int Packet::GrabInt(char * mb) {
  int* val_ptr = (int *)mb;
  return ntohl(*val_ptr);
}

/***********************************************************/
/*  Name : Grab Double                                     */
/*                                                         */
/*  Description: Grab a double value from the message      */
/*                buffer, which is in network byte order.  */
/*               Must convert back to host byte order.     */
/*                                                         */
/*  Input : pointer to beginning of portion of message     */
/*           buffer to grab value from.                    */
/*  Return Value : Returns double value                    */
/***********************************************************/

double Packet::GrabDouble(char * mb) {
  double* val_ptr = (double *)mb;
  return ntohd(*val_ptr);
}

/***********************************************************/
/*  Name : Grab Float                                      */
/*                                                         */
/*  Description: Grab a float value from the message       */
/*                buffer, which is in network byte order.  */
/*               Must convert back to host byte order.     */
/*                                                         */
/*  Input : pointer to beginning of portion of message     */
/*           buffer to grab value from.                    */
/*  Return Value : Returns float value                     */
/***********************************************************/

float Packet::GrabFloat(char * mb) {
  float* val_ptr = (float *)mb;
  return ntohf(*val_ptr);
}

/***********************************************************/
/*  Name : Grab CharArray                                  */
/*                                                         */
/*  Description: Grab a character array the message buffer */
/*                                                         */
/*  Input : mb - pointer to beginning of portion of message*/
/*           buffer to grab value from.                    */
/*          dst - ptr to where copy (place) data read      */
/*          size - how much to copy                        */
/*  Return Value : Returns float value                     */
/***********************************************************/

void Packet::GrabCharArray(char * mb, char * dst, int size) {
  if (dst == NULL) {
    return;
  }
  int i;
  for (i=0; i < size; i++) {
    dst[i] = mb[i];
  }  
}

/***********************************************************/
/*  Name : Grab TRVector                                   */
/*                                                         */
/*  Description: Grab a TRVector value from the message    */
/*                buffer, which is in network byte order.  */
/*               Must convert back to host byte order.     */
/*                                                         */
/*  Input : mb - pointer to beginning of portion of message*/
/*           buffer to grab value from.                    */
/*          val - Reference to TRVector (for saving val)   */
/*  Return Value : None                                    */
/***********************************************************/

void Packet::GrabTRVector(char * mb, TRVector & val) {
  // order X, Y, Z
  val.X(GrabDouble(mb));
  int ptr = sizeof(double);
  val.Y(GrabDouble(mb+ptr));
  ptr += sizeof(double);
  val.Z(GrabDouble(mb+ptr));  
}

/***********************************************************/
/*  Name : Grab TRPoint                                    */
/*                                                         */
/*  Description: Grab a TRPoint value from the message     */
/*                buffer, which is in network byte order.  */
/*               Must convert back to host byte order.     */
/*                                                         */
/*  Input : mb - pointer to beginning of portion of message*/
/*           buffer to grab value from.                    */
/*          val - Reference to TRPoint (for saving val)    */
/*  Return Value : None                                    */
/***********************************************************/

void Packet::GrabTRPoint(char * mb, TRPoint & val) {
  // order X, Y, Z
  val.X(GrabDouble(mb));
  int ptr = sizeof(double);
  val.Y(GrabDouble(mb+ptr));
  ptr += sizeof(double);
  val.Z(GrabDouble(mb+ptr));  
}

/***********************************************************/
/*  Name : Grab TRQuaternion                               */
/*                                                         */
/*  Description: Grab a TRQuaternion value from the message*/
/*                buffer, which is in network byte order.  */
/*               Must convert back to host byte order.     */
/*                                                         */
/*  Input : mb - pointer to beginning of portion of message*/
/*           buffer to grab value from.                    */
/*          val - Reference to TRQuaternion (for saving)   */
/*  Return Value : None                                    */
/***********************************************************/

void Packet::GrabTRQuaternion(char * mb, TRQuaternion & val) {
  // order vector (x,y,z) then double
  TRVector temp;
  GrabTRVector(mb, temp);
  val.qV(temp);
  // GrabTRVector(mb, val.qV()) might work, should try
  
  val.qW(GrabDouble(mb+temp.GetSize()));
}

/**********************************************************/
/*              Private Functions                         */
/**********************************************************/

/***********************************************************/
/*  Name : htond                                           */
/*                                                         */
/*  Description: Convert double value from host byte order */
/*                to network byte order.                   */
/*                                                         */
/*  Input : double value                                   */
/*  Return Value : double value in network byte order      */
/***********************************************************/

double Packet::htond(const double& d) {
  if (BYTESWAPPING) {
    double retval;
    char *cptrout = (char*)&retval;
    char *cptrin = (char*)&d;
    cptrout[0] = cptrin[7];
    cptrout[1] = cptrin[6];
    cptrout[2] = cptrin[5];
    cptrout[3] = cptrin[4];
    cptrout[4] = cptrin[3];
    cptrout[5] = cptrin[2];
    cptrout[6] = cptrin[1];
    cptrout[7] = cptrin[0];
    return retval;
  }
  else {
    return d;
  }
}

/***********************************************************/
/*  Name : htonf                                           */
/*                                                         */
/*  Description: Convert float value from host byte order  */
/*                to network byte order.                   */
/*                                                         */
/*  Input : float value                                    */
/*  Return Value : float value in network byte order       */
/***********************************************************/

float Packet::htonf(const float& f) {
  if (BYTESWAPPING) {
    float retval;
    char *cptrout = (char*)&retval;
    char *cptrin = (char*)&f;
    cptrout[0] = cptrin[3];
    cptrout[1] = cptrin[2];
    cptrout[2] = cptrin[1];
    cptrout[3] = cptrin[0];
    return retval;
  }
  else {
    return f;
  }
}

