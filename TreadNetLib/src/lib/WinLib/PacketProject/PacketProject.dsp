# Microsoft Developer Studio Project File - Name="PacketProject" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=PacketProject - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "PacketProject.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "PacketProject.mak" CFG="PacketProject - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "PacketProject - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "PacketProject - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "PacketProject - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "PacketProject - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GR /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"Debug\TRNetPacketLibrary.lib"

!ENDIF 

# Begin Target

# Name "PacketProject - Win32 Release"
# Name "PacketProject - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\Packets\ByteArray_Data.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\ByteArray_Packet.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\Collision_Data.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\Collision_Packet.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\Command_Data.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\Command_Packet.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\FSRSignals_Data.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\FSRSignals_Packet.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\GetUserState_Packet.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\NewPos_Data.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\NewPos_Packet.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\NewPosResponse_Packet.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\Packet.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\Registration_Data.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\Registration_Packet.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\ServerRequest_Data.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\ServerRequest_Packet.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\Stairs_Data.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\Stairs_Packet.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\Tracker_Data.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\Tracker_Packet.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\Treadport_Types.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\TRQuaternion.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\TRVector.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\UserState_Data.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\UserState_Packet.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\VirtualState_Data.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Packets\VirtualState_Packet.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\Packets\ByteArray_Data.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\ByteArray_Packet.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\Collision_Data.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\Collision_Packet.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\Command_Data.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\Command_Packet.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\FSRSignals_Data.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\FSRSignals_Packet.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\GetUserState_Data.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\GetUserState_Packet.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\NewPos_Data.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\NewPos_Packet.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\NewPosResponse_Data.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\NewPosResponse_Packet.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\Packet.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\Packet_Data.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\Registration_Data.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\Registration_Packet.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\ServerRequest_Data.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\ServerRequest_Packet.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\Stairs_Data.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\Stairs_Packet.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\Tracker_Data.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\Tracker_Packet.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\Treadport_Types.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\TRQuaternion.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\TRVector.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\UserState_Data.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\UserState_Packet.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\VirtualState_Data.h
# End Source File
# Begin Source File

SOURCE=..\..\Packets\VirtualState_Packet.h
# End Source File
# End Group
# End Target
# End Project
