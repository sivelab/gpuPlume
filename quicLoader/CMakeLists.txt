cmake_minimum_required(VERSION 2.6)

include_directories (${GPUPLUME_SOURCE_DIR}/fileParser/standard) 

add_library (quicLoader
	    QUICDataFile.h
	    QUBuildings.h
	    QUBuildings.cpp
	    QUMetParams.h
	    QUMetParams.cpp
	    QUSimparams.h
	    QUSimparams.cpp
	    QPParams.h
	    QPParams.cpp
	    QPSource.h
	    QPSource.cpp
)
