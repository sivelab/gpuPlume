#ifndef STAIRS_DATA_H
#define STAIRS_DATA_H

#include "Packet_Data.h"

enum TRStairs_Dir {
  IGNORE_STAIRS,
  UP_STAIRS,
  DOWN_STAIRS
};

enum TRStairs_Cmd {
  NO_STAIRS_CMD,
  SWITCH_STAIRSLOPE_CMD,
  INC_CONSTFACT_CMD,
  DEC_CONSTFACT_CMD
};

// This class is for stairs packet data

class Stairs_Data : public Packet_Data {
 public:
  friend class Stairs_Packet;

  // Destructor
  ~Stairs_Data();
  
  // default constructor
  Stairs_Data();

  // copy constructor
  Stairs_Data(const Stairs_Data&);

  // Assignment operator
  void operator = (const Stairs_Data&);

  // Accessor functions
  inline TRStairs_Cmd GetCmd() const {return m_stairs_cmd;}
  inline TRStairs_Dir GetStairsDir() const {return m_stairs_dir;}

  // Set functions
  
  // Set Commands
  inline void SetSwitchStairsCmd() {m_stairs_cmd = SWITCH_STAIRSLOPE_CMD;}
  inline void SetIncConstFactCmd() {m_stairs_cmd = INC_CONSTFACT_CMD;    }
  inline void SetDecConstFactCmd() {m_stairs_cmd = DEC_CONSTFACT_CMD;    }
  inline void SetNoStairsCmd()     {m_stairs_cmd = NO_STAIRS_CMD;        }

  // Set Stairs Direction
  inline void SetUpStairsDir()     {m_stairs_dir = UP_STAIRS;    }
  inline void SetDownStairsDir()   {m_stairs_dir = DOWN_STAIRS;  }
  inline void SetIgnoreStairsDir() {m_stairs_dir = IGNORE_STAIRS;}

  // Already provides uint GetSize();

 protected:
  // Do not want user to be able to set any command or dir values
  inline void SetCmd(TRStairs_Cmd cmd)       {m_stairs_cmd = cmd;}
  inline void SetStairsDir(TRStairs_Dir dir) {m_stairs_dir = dir;}


  // Inherits variable ushort m_size

 private:
  // Stairs command
  TRStairs_Cmd m_stairs_cmd;

  // Stair direction
  TRStairs_Dir m_stairs_dir;
};

#endif
