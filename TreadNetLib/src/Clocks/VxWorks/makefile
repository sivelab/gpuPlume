#This is our specific makefile for automatically installing the libraries
# if specifying target install (not for compiling - that must be done
# through the appropriate program)

ROOT = ../../../
include $(ROOT)/Make/MakeDefs

ARCH = intel_vx5.5

LIB = VxSyncLib.a

#supports only install and install_dist

default clean veryclean :

install: libdir
	@echo "  Installing VxWorks Libraries/Headers"
	-cp SyncProject/PENTIUM4gnu/$(LIB) $(ROOT)/libs/$(ARCH)
	-cp vxSync.h $(ROOT)/include/vxSync.h
	-@chmod g+rw $(ROOT)/libs/$(ARCH)/$(LIB)
	-@chmod g+rw $(ROOT)/include/vxSync.h

install_dist: distdirs
	cvs tag -c $(DIST_TAG)
	-cp SyncProject/PENTIUM4gnu/$(LIB) $(FULL_DIST_DIR)/libs/$(ARCH)
	-cp vxSync.h $(FULL_DIST_DIR)/include/vxSync.h

libdir:
	-@mkdir -p $(ROOT)/libs/$(ARCH)
	-@chmod g+rwxs $(ROOT)/libs/$(ARCH)

distdirs : 
	-@mkdir -p $(FULL_DIST_DIR)
	-@mkdir -p $(FULL_DIST_DIR)/libs/$(ARCH)
