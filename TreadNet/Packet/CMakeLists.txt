cmake_minimum_required (VERSION 2.6)

# Create a library called "Packet"
add_library (Packet Packet.cpp ByteArray_Data.cpp Command_Data.cpp GetUserState_Packet.cpp ServerRequest_Packet.cpp Tracker_Packet.cpp UserState_Data.cpp ByteArray_Packet.cpp Command_Packet.cpp NewPos_Data.cpp Registration_Data.cpp Stairs_Data.cpp Treadport_Types.cpp UserState_Packet.cpp Collision_Data.cpp FSRSignals_Data.cpp NewPos_Packet.cpp Registration_Packet.cpp Stairs_Packet.cpp TRQuaternion.cpp VirtualState_Data.cpp Collision_Packet.cpp FSRSignals_Packet.cpp NewPosResponse_Packet.cpp ServerRequest_Data.cpp Tracker_Data.cpp TRVector.cpp VirtualState_Packet.cpp)

