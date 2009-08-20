#include "Treadport_Types.h"

#ifdef RTI_VXWORKS
// Define strdup for vxworks
char * strdup(const char *str) {
  char * new_string = new char[strlen(str)+1];
  strcpy(new_string,str);
  return new_string;
}
#endif


/*********************************************************/
/*  Name : TRTerrain ++ operator (prefix and postfix)    */
/*                                                       */
/*  Description: Add one to the terrain                  */
/*                                                       */
/*  Input : TRTerrain&, int                              */
/*  Return Value : TRTerrain                             */
/*********************************************************/

TRTerrain operator ++ (TRTerrain& t, int) {
  return t = (TRTerrain)(t + 1);
}

TRTerrain operator ++ (TRTerrain& t) {
  return t = (TRTerrain)(t + 1);
}
