#ifndef COLLISION_PACKET_H
#define COLLISION_PACKET_H

#include "Packet.h"
#include "Collision_Data.h"

// A quick definition to simplify the code below
#define COLL_DATA dynamic_cast<Collision_Data*>(m_data)

class Collision_Packet : public Packet {
 public:
  // Destructor
  ~Collision_Packet();

  // Default Constructor
  Collision_Packet();

  // Constructor w/ data passed in
  // Note: Data will be copied to a new Collision_data class
  Collision_Packet(const Collision_Data&);

  /*************************************************************/
  /* Functions already contained (see Packet.h for details)    */
  /*                                                           */
  /* int GetSize();                                            */
  /* uint GetType();                                           */
  /* BOOL IsRegistered();                                      */
  /*************************************************************/

  // Returns Class Name, for identifying packet
  const char * ClassName() {return Collision_Packet::m_class_name;}

  // Packet_Data functions - these affect the values of m_data, but
  //  does not delete/create m_data.
  int GetCopyOfData(Packet_Data * CopyDat);
  int SetData(Packet_Data * NewDat);

  // Access Functions - just pass along to the data class
  inline TRCollFlag GetCollisionFlag() const {return COLL_DATA->GetCollisionFlag();}
  inline TRVector GetCollisionNormal() const {return COLL_DATA->GetCollisionNormal();}
  inline TRPoint GetCollisionPos() const {return COLL_DATA->GetCollisionPos();}
  
  // Set functions - again pass along
  inline void SetNoCollision() {COLL_DATA->SetNoCollision();}
  inline void SetCollision() {COLL_DATA->SetCollision();}

  inline void SetCollisionNormal(const TRVector& cn) {COLL_DATA->SetCollisionNormal(cn);}
  inline void SetCollisionPos(const TRVector& cp) {COLL_DATA->SetCollisionPos(cp);}

 protected:
  // Do not want users to set any flag value
  inline void SetCollisionFlag(TRCollFlag f) {COLL_DATA->SetCollisionFlag(f);}

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
