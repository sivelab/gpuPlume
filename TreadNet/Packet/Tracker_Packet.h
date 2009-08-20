#ifndef TRACKER_PACKET_H
#define TRACKER_PACKET_H

#include "Packet.h"
#include "Tracker_Data.h"

// A quick definition to simplify the code below
#define TRACKER_DATA dynamic_cast<Tracker_Data*>(m_data)

class Tracker_Packet : public Packet {
 public:
  // Destructor
  ~Tracker_Packet();

  // Default Constructor
  Tracker_Packet();

  // Constructor w/ data passed in
  // Note: Data will be copied to a new Tracker_data class
  Tracker_Packet(const Tracker_Data&);

  /*************************************************************/
  /* Functions already contained (see Packet.h for details)    */
  /*                                                           */
  /* int GetSize();                                            */
  /* uint GetType();                                           */
  /* BOOL IsRegistered();                                      */
  /*************************************************************/

  // Returns Class Name, for identifying packet
  const char * ClassName() {return Tracker_Packet::m_class_name;}

  // Packet_Data functions - these affect the values of m_data, but
  //  does not delete/create m_data.
  int GetCopyOfData(Packet_Data * CopyDat);
  int SetData(Packet_Data * NewDat);

  // Access Functions - just pass along to the data class
  //  (both whole and particular)
  inline TRPoint GetPos() const              {return TRACKER_DATA->GetPos();}
  inline TRQuaternion GetOrientation() const {return TRACKER_DATA->GetOrientation();}
  
  inline double   GetPosX() const         {return TRACKER_DATA->GetPosX();}
  inline double   GetPosY() const         {return TRACKER_DATA->GetPosY();}
  inline double   GetPosZ() const         {return TRACKER_DATA->GetPosZ();}
  inline double   GetOrientRot() const    {return TRACKER_DATA->GetOrientRot();}
  inline TRVector GetOrientVector() const {return TRACKER_DATA->GetOrientVector();}

  // Set functions - again pass along
  inline void SetPos(const TRPoint& pos)                 {TRACKER_DATA->SetPos(pos);}
  inline void SetOrientation(const TRQuaternion& orient) {TRACKER_DATA->SetOrientation(orient);}

  inline void SetPosX(double x)                  {TRACKER_DATA->SetPosX(x);}
  inline void SetPosY(double y)                  {TRACKER_DATA->SetPosY(y);}
  inline void SetPosZ(double z)                  {TRACKER_DATA->SetPosZ(z);}
  inline void SetOrientRot(double r)             {TRACKER_DATA->SetOrientRot(r);}
  inline void SetOrientVector(const TRVector &v) {TRACKER_DATA->SetOrientVector(v);}
  
 protected:
  // Fill the message buffer up with the stair data
  virtual int Fill(char *);

  // Read a stair packet from message buffer
  virtual int ReadPacket(char *);

  /*************************************************************/
  /* Procedures from Packet                                    */
  /*                                                           */
  /* BOOL WaitingForResponse();                                */
  /* int FillShort(char *, const short&);                      */
  /* ... FillInt, FillDouble, FillFloat, and FillString        */
  /* short GrabShort(char *);                                  */
  /* ..... GrabInt, GrabDouble, GrabFloat                      */
  /*************************************************************/

  // Variable Section

  // Class Name
  static char* m_class_name;

  /*************************************************************/
  /*  Variables already contained in class from Packet         */
  /*                                                           */
  /* Packet_Data* m_data;                                      */
  /* unsigned int m_type;                                      */
  /* BOOL m_wait_response;                                     */
  /* BOOL m_registered;                                        */
  /*************************************************************/

};

#endif
