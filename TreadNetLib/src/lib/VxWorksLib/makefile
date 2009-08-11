#This is our specific makefile for automatically installing the libraries
# if specifying target install (not for compiling - that must be done
# through the appropriate program)

ROOT = ../../../
include $(ROOT)/Make/MakeDefs

ARCH = intel_vx5.5

CLIB = TRNetClientControllerLib.a
SLIB = TRNetServerControllerLib.a
PLIB = TRNetPacketLibrary.a

#supports only install and install_dist

default clean veryclean :

install: libdir
	@echo "  Installing VxWorks Libraries"
	-cp ClientProject/PENTIUM4gnu/$(CLIB) $(ROOT)/libs/$(ARCH)
	-@chmod g+rw $(ROOT)/libs/$(ARCH)/$(CLIB)
	-cp ServerProject/PENTIUM4gnu/$(SLIB) $(ROOT)/libs/$(ARCH)
	-@chmod g+rw $(ROOT)/libs/$(ARCH)/$(SLIB)
	-cp PacketProject/PENTIUM4gnu/$(PLIB) $(ROOT)/libs/$(ARCH)
	-@chmod g+rw $(ROOT)/libs/$(ARCH)/$(PLIB)

install_dist: distdirs
	cvs tag -c $(DIST_TAG)
	-cp ClientProject/PENTIUM4gnu/$(CLIB) $(FULL_DIST_DIR)/libs/$(ARCH)
	-cp ServerProject/PENTIUM4gnu/$(SLIB) $(FULL_DIST_DIR)/libs/$(ARCH)
	-cp PacketProject/PENTIUM4gnu/$(PLIB) $(FULL_DIST_DIR)/libs/$(ARCH)

libdir:
	-@mkdir -p $(ROOT)/libs/$(ARCH)
	-@chmod g+rwxs $(ROOT)/libs/$(ARCH)

distdirs : 
	-@mkdir -p $(FULL_DIST_DIR)
	-@mkdir -p $(FULL_DIST_DIR)/libs/$(ARCH)
