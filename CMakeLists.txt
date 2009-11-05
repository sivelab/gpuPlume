# The name of our project is "GPUPLUME". CMakeLists files in this project can 
# refer to the root source directory of the project as ${GPUPLUME_SOURCE_DIR} and 
# to the root binary directory of the project as ${GPUPLUME_BINARY_DIR}. 
cmake_minimum_required (VERSION 2.6) 
project (GPUPLUME) 

set(GPUPLUME_MAJOR_VERSION 0)
set(GPUPLUME_MINOR_VERSION 9)
set(GPUPLUME_PATCH_VERSION 0)
set(GPUPLUME_SOVERSION 55)

set(GPUPLUME_VERSION ${GPUPLUME_MAJOR_VERSION}.${GPUPLUME_MINOR_VERSION}.${GPUPLUME_PATCH_VERSION})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Add defines for Windows.
if(WIN32)
  add_definitions(-DWIN32)
endif(WIN32)

# Find OpenGL
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

# If we are in Windows, set the GLEW path manually.
if(WIN32)
  set(GLEW_INCLUDE_PATH "${GPUPLUME_SOURCE_DIR}/glew/include")
  set(GLEW_LIBRARY "${GPUPLUME_SOURCE_DIR}/glew/lib")
else(WIN32)
  find_package(GLEW REQUIRED)
endif(WIN32)

include_directories(${OPENGL_INCLUDE_DIR})
include_directories(${GLUT_INCLUDE_DIRS})
include_directories(${GLEW_INCLUDE_PATH})

set(LIBS ${LIBS} ${OPENGL_LIBRARIES})
set(LIBS ${LIBS} ${GLUT_LIBRARIES})
set(LIBS ${LIBS} ${GLEW_LIBRARY})

# set(CMAKE_CXX_FLAGS "-g3 -Wall")

if(APPLE)
    find_library(CARBON_LIBRARY Carbon)
    find_library(COCOA_LIBRARY Cocoa)
endif(APPLE)

# Recurse into the "fileParser" subdirectory to build the fileParser
# library. This does not actually cause another cmake executable to
# run. The same process will walk through the project's entire
# directory structure.  
add_subdirectory (fileParser)

# build a library for handling all of the QUIC files for reading and
# writing.
add_subdirectory (quicLoader)

# Recuse into the TreadNet subdirectory to build the libraries to
# connect to the TPAWT.
add_subdirectory (TreadNet/Packet)
add_subdirectory (TreadNet/Controller)

# Make sure the compiler can find include files from our Hello library. 
include_directories (${GPUPLUME_SOURCE_DIR}/fileParser/standard) 

# Make sure the linker can find the Hello library once it is built. 
link_directories (${GPUPLUME_SOURCE_DIR}/fileParser)

# Make sure the compiler can find the include files for the TreadNet library.
include_directories (${GPUPLUME_SOURCE_DIR}/TreadNet/Packet)
include_directories (${GPUPLUME_SOURCE_DIR}/TreadNet/Controller)

# Make sure the linker can find the TreadNetLib.
link_directories (${GPUPLUME_SOURCE_DIR}/TreadNet/Packet)
link_directories (${GPUPLUME_SOURCE_DIR}/TreadNet/Controller)

# Add executable called "gpuPlume" that is built from the source files 
# The source files are placed in the list after the executable name. 
add_executable (gpuPlume 
	       gpuPlume.cpp 
	       ArgumentParsing.h 
	       ArgumentParsing.cpp
	       framebufferObject.cpp 
	       renderbuffer.cpp 
	       GLSL.cpp 
	       glErrorUtil.cpp 
	       particleControl.cpp 
	       pointEmitter.cpp 
	       sphereEmitter.cpp 
	       particleEmitter.cpp 
	       displayControl.cpp 
	       plumeControl.cpp 
	       Timer.cpp 
	       collectionBox.cpp 
	       util.cpp 
	       simulation.cpp 
	       streamLine.cpp 
	       Random.cpp 
	       nonGaussianModel.cpp 
	       GaussianModel.cpp 
	       Gaussian_2shaders_Model.cpp 
	       ReflectionModel.cpp 
	       MultipleBuildingsModel.cpp 
	       GeomTest.cpp 
	       PathLine.cpp 
	       Contour.cpp 
	       VisualPlane.cpp 
	       lineEmitter.cpp 
	       planeEmitter.cpp 
	       IsoSurface.cpp 
	       NetworkManager.cpp 
	       broadcaster.cpp 
	       receiver.cpp 
	       InetAddress.cpp 
	       TreadportManager.cpp 
	       graphicsUtil.cpp 
	       rgb.cpp 
	       CmdOptionInterpreter.cpp
) 

# Link the executable to the TreadNet library. 
target_link_libraries (gpuPlume NetController)
target_link_libraries (gpuPlume ClientController)
target_link_libraries (gpuPlume Packet)
target_link_libraries (gpuPlume ServerController)

# Link the executable to the file parser library. 
target_link_libraries (gpuPlume StandardFileParser) 
target_link_libraries (gpuPlume ${LIBS})

