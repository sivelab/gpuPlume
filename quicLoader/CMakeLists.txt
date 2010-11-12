cmake_minimum_required(VERSION 2.6)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3 -Wall")

include_directories (${GPUPLUME_SOURCE_DIR}/fileParser/standard) 

add_library (quicLoader
	    QUICDataFile.h
	    QUBuildings.h
	    QUBuildings.cpp
	    QUMetParams.h
	    QUMetParams.cpp
	    QUSensor.h
	    QUSensor.cpp
	    QUSimparams.h
	    QUSimparams.cpp
	    QPParams.h
	    QPParams.cpp
	    QPSource.h
	    QPSource.cpp
)
