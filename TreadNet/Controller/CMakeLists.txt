cmake_minimum_required (VERSION 2.6)

# Make sure the compiler can find include files from the Packet library.
include_directories (../Packet)

# Make sure the linker can find the Packet library once it is built. 
link_directories (../Packet)

# Create the ClientController library
add_library (ClientController Client_Controller.cpp)

# Create the ServerController library
add_library (ServerController Server_Controller.cpp)

# Create the NetController library
add_library (NetController Net_Controller.cpp)
