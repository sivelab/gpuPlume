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

SET(LIBSIVELAB_PATH "${CMAKE_SOURCE_DIR}/../libsivelab" CACHE STRING "Path to the directory where the libsivelab repository is located.")
INCLUDE_DIRECTORIES (${LIBSIVELAB_PATH}) 
LINK_DIRECTORIES (${LIBSIVELAB_PATH}/build/quicutil) 
LINK_DIRECTORIES (${LIBSIVELAB_PATH}/build/util) 

# From http://stackoverflow.com/questions/2368811/how-to-set-warning-level-in-cmake
if(MSVC)
  # Force to always compile with W4
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  endif()
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Update if necessary
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic")
endif()

# Find OpenGL
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

# If we are in Windows, set the GLEW path manually.
if(WIN32)
  add_definitions(-DWIN32)
  set(GLEW_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/extLibs/glew/include")
  set(GLEW_LIBRARY "${CMAKE_SOURCE_DIR}/extLibs/glew/lib/glew32s.lib")

  set(GETOPT_FILES "getopt_win32.c")
else(WIN32)
  if(APPLE)
    set(GLEW_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/extLibs/glew/include")
    set(GLEW_LIBRARY "${CMAKE_SOURCE_DIR}/extLibs/glew/lib/libGLEW.a")
  endif(APPLE)
endif(WIN32)
FIND_PACKAGE(GLEW REQUIRED)

IF (APPLE)
FIND_PACKAGE(PNG REQUIRED)
IF(PNG_FOUND)
   INCLUDE_DIRECTORIES(${PNG_INCLUDE_DIR})
   # TARGET_LINK_LIBRARIES(imageviewer ${PNG_LIBRARY})
ENDIF(PNG_FOUND)
ENDIF(APPLE)

include_directories(${OPENGL_INCLUDE_DIR})
include_directories(${GLUT_INCLUDE_DIRS})
include_directories(${GLEW_INCLUDE_PATH})

set(LIBS ${LIBS} ${OPENGL_LIBRARIES})
set(LIBS ${LIBS} ${GLUT_LIBRARIES})
set(LIBS ${LIBS} ${GLEW_LIBRARY})

# set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-DGPUPLUME_VERSION=\"${GPUPLUME_VERSION}\"")

if(APPLE)
    find_library(CARBON_LIBRARY Carbon)
    find_library(COCOA_LIBRARY Cocoa)
endif(APPLE)

# Recuse into the TreadNet subdirectory to build the libraries to
# connect to the TPAWT.
add_subdirectory (TreadNet/Packet)
add_subdirectory (TreadNet/Controller)

# Make sure the compiler can find the include files for the TreadNet library.
include_directories (${GPUPLUME_SOURCE_DIR}/TreadNet/Packet)
include_directories (${GPUPLUME_SOURCE_DIR}/TreadNet/Controller)

# Make sure the linker can find the TreadNetLib.
link_directories (${GPUPLUME_SOURCE_DIR}/TreadNet/Packet)
link_directories (${GPUPLUME_SOURCE_DIR}/TreadNet/Controller)

# Add executable called "gpuPlume" that is built from the source files 
# The source files are placed in the list after the executable name. 
add_executable (gpuPlume 
               ${GETOPT_FILES}
	       gpuPlume.cpp 
	       ArgumentParsing.cpp ArgumentParsing.h 
	       framebufferObject.cpp framebufferObject.h
	       renderbuffer.cpp renderbuffer.h
	       GLSL.cpp GLSL.h
	       glErrorUtil.cpp glErrorUtil.h
	       particleControl.cpp particleControl.h
	       pointEmitter.cpp pointEmitter.h
	       sphereEmitter.cpp sphereEmitter.h
	       particleEmitter.cpp particleEmitter.h
	       displayControl.cpp displayControl.h
	       plumeControl.cpp plumeControl.h
	       Timer.cpp Timer.h
	       collectionBox.cpp collectionBox.h
	       util.cpp util.h
	       simulation.cpp 
	       streamLine.cpp 
	       Random.cpp Random.h
#	       nonGaussianModel.cpp 
#	       GaussianModel.cpp 
#	       Gaussian_2shaders_Model.cpp 
#	       ReflectionModel.cpp 
	       MultipleBuildingsModel.cpp MultipleBuildingsModel.h
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
	       CmdOptionInterpreter.cpp CmdOptionInterpreter.h
	       writepng.cpp writepng.h
	       PNGImage.cpp PNGImage.h
) 

# Link the executable to the TreadNet library. 
target_link_libraries (gpuPlume NetController)
target_link_libraries (gpuPlume ClientController)
target_link_libraries (gpuPlume Packet)
target_link_libraries (gpuPlume ServerController)

# Link the executable to the file parser library. 
# target_link_libraries (gpuPlume quicLoader) 
# target_link_libraries (gpuPlume StandardFileParser) 
target_link_libraries (gpuPlume sive-quicutil) 
target_link_libraries (gpuPlume ${LIBS})

target_link_libraries (gpuPlume ${PNG_LIBRARY})
target_link_libraries (gpuPlume z)

add_subdirectory (populationSampleOpt)

# add_executable (cudaPlume
#  cudaPlume.cpp
#   util.cpp util.h
# )
# target_link_libraries (cudaPlume sive-quicutil) 
