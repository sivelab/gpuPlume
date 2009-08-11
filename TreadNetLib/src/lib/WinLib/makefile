#This is our specific makefile for automatically installing the libraries
# if specifying target install (not for compiling - that must be done
# through the appropriate program)

ROOT = ../../../
include $(ROOT)/Make/MakeDefs

ARCH = intel_win

CLIB = TRNetClientLibrary.lib
SLIB = TRNetServerLibrary.lib
PLIB = TRNetPacketLibrary.lib

#supports only install and install_dist

default clean veryclean :

install: libdir
	@echo "  Installing Windows Libraries"
	-cp ClientProject/Debug/$(CLIB) $(ROOT)/libs/$(ARCH)
	-@chmod g+rw $(ROOT)/libs/$(ARCH)/$(CLIB)
	-cp ServerProject/Debug/$(SLIB) $(ROOT)/libs/$(ARCH)
	-@chmod g+rw $(ROOT)/libs/$(ARCH)/$(SLIB)
	-cp PacketProject/Debug/$(PLIB) $(ROOT)/libs/$(ARCH)
	-@chmod g+rw $(ROOT)/libs/$(ARCH)/$(PLIB)

install_dist: distdirs
	cvs tag -c $(DIST_TAG)
	-cp ClientProject/Debug/$(CLIB) $(FULL_DIST_DIR)/libs/$(ARCH)
	-cp ServerProject/Debug/$(SLIB) $(FULL_DIST_DIR)/libs/$(ARCH)
	-cp PacketProject/Debug/$(PLIB) $(FULL_DIST_DIR)/libs/$(ARCH)

libdir:
	-@mkdir -p $(ROOT)/libs/$(ARCH)
	-@chmod g+rwxs $(ROOT)/libs/$(ARCH)

distdirs : 
	-@mkdir -p $(FULL_DIST_DIR)
	-@mkdir -p $(FULL_DIST_DIR)/libs/$(ARCH)
