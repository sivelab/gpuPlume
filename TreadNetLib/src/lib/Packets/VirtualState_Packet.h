#ifndef VIRTUALSTATE_PACKET_H
#define VIRTUALSTATE_PACKET_H

#include "Packet.h"
#include "VirtualState_Data.h"

// A quick definition to simplify the code below
#define VIRTSTATE_DATA dynamic_cast<VirtualState_Data*>(m_data)

class VirtualState_Packet : public Packet {
  friend class NewPos_Packet;

 public:
  // Destructor
  ~VirtualState_Packet();

  // Default Constructor
  VirtualState_Packet();

  // Constructor w/ data passed in
  // Note: Data will be copied to a new VirtualState_data class
  VirtualState_Packet(const VirtualState_Data&);

  /*************************************************************/
  /* Functions already contained (see Packet.h for details)    */
  /*                                                           */
  /* int GetSize();                                            */
  /* uint GetType();                                           */
  /* BOOL IsRegistered();                                      */
  /*************************************************************/

  // Returns Class Name, for identifying packet
  const char * ClassName() {return VirtualState_Packet::m_class_name;}

  // Packet_Data functions - these affect the values of m_data, but
  //  does not delete/create m_data.
  int GetCopyOfData(Packet_Data * CopyDat);
  int SetData(Packet_Data * NewDat);

  // Access Functions - just pass along to the data class
  inline TRVector GetSurfaceNormal() const {return VIRTSTATE_DATA->GetSurfaceNormal();}
  inline TRTerrain GetTerrain() const {return VIRTSTATE_DATA->GetTerrain();}
  inline TRVector GetVirtualForce() const {return VIRTSTATE_DATA->GetVirtualForce();}
  inline float GetMovementScale() const {return VIRTSTATE_DATA->GetMovementScale();}
  inline float GetUserScale() const {return VIRTSTATE_DATA->GetUserScale();}
  
  // Set functions - again pass along
  inline void SetSurfaceNormal(const TRVector& sn) {VIRTSTATE_DATA->SetSurfaceNormal(sn);}
  inline void SetTerrain(TRTerrain ter) {VIRTSTATE_DATA->SetTerrain(ter);}
  inline void SetVirtualForce(const TRVector& vf) {VIRTSTATE_DATA->SetVirtualForce(vf);}
  inline void SetMovementScale(float ms) {VIRTSTATE_DATA->SetMovementScale(ms);}
  inline void SetUserScale(float us) {VIRTSTATE_DATA->SetUserScale(us);}

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
