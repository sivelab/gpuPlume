#ifndef PACKET_H
#define PACKET_H

#include <sys/types.h>

#ifndef WIN32
#include <netinet/in.h>
#endif

#include "Packet_Data.h"

// determine if byteswapping is needed for packets to
//  conform to network byte order
#define BYTESWAPPING (htons((short)1) == ((short)256))

class Packet {
  friend class Net_Controller;
  // Make copy routine???
 public:
  // Destructor
  virtual ~Packet();

  // Default constructor
  Packet();

  // Returns Class Name, used for identifying packet
  // The class name should be the same for all packets of a given type,
  //  thus make it a static member in the child classes.
  virtual const char * ClassName() = 0;

  // Access Functions

  // Returns size of packet
  inline ushort GetSize() {return m_data->GetSize();}
  
  // Gets copy of data
  virtual int GetCopyOfData(Packet_Data *) = 0;

  // Return if registered
  inline BOOL IsRegistered() {return m_registered;} 

  // Return packet type
  inline uint GetType() {return m_type;}

  // SetFunction (only public one is for data)
  virtual int SetData(Packet_Data *) = 0;

 protected:
  // Set Type (Just increments based on server registration)
  virtual void SetType(uint val) {m_type = val;}

  // SetRegistered
  inline void SetRegistered() {m_registered = TRUE;}

  // Sending functions

  // Notify network controller that packet is waiting for response
  inline BOOL WaitingForResponse() {return m_wait_response;}

  // Fill the message buffer up with data
  virtual int Fill(char *) = 0;

  // Fill routines for putting different types into buffer (short, 
  //  int, double, float, and string)
  virtual int FillShort(char *, short);
  virtual int FillInt(char *, int);
  virtual int FillDouble(char *, double);
  virtual int FillFloat(char *, float);
  virtual int FillString(char *, const char *);
  virtual int FillChar(char *, char);
  virtual int FillTRVector(char *, const TRVector&);
  virtual int FillTRPoint(char *, const TRPoint&);
  virtual int FillTRQuaternion(char *, const TRQuaternion&);
  virtual int FillCharArray(char * mb, const char * val, int);
  
  // Receiving packet functions

  // Read packet received (gives pointer to message buffer section
  //  containing packet)
  virtual int ReadPacket(char *) = 0;
  
  // Fills the data values from the message buffer for different types
  //  (short, int, double, float)
  virtual short GrabShort(char *);
  virtual int GrabInt(char *);
  virtual double GrabDouble(char *);
  virtual float GrabFloat(char *);
  virtual void GrabCharArray(char * mb, char * dst, int); 
  // these are slightly different to avoid having to create multiple
  //  copies of a class.
  virtual void GrabTRVector(char *, TRVector &);
  virtual void GrabTRPoint(char *, TRPoint &);
  virtual void GrabTRQuaternion(char *, TRQuaternion &);

 protected:
  // Packet data portion
  Packet_Data* m_data;

  // Type, so can identify which packets are part of message
  uint m_type;

  // Set if packet needs a response (unset when receives it)
  BOOL m_wait_response;

  // Set if registered this packet
  BOOL m_registered;

 private:
  // Function to change double to network byte order
  double htond(const double&);

  // Function to change float to network byte order
  float htonf(const float&);

  // Functions to perform the reverse operation
  //  (network to host)
  inline double ntohd(const double& val) {return htond(val);}
  inline float ntohf(const float& val) {return htonf(val);}
};

#endif
