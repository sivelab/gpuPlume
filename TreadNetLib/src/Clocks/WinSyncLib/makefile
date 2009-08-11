#This is our specific makefile for automatically installing the libraries
# if specifying target install (not for compiling - that must be done
# through the appropriate program)

ROOT = ../../../
include $(ROOT)/Make/MakeDefs

ARCH = intel_win

LIB = WinSyncLib.lib

#supports only install and install_dist

default clean veryclean :

install: libdir
	@echo "  Installing Windows Libraries/Headers"
	-cp Debug/$(LIB) $(ROOT)/libs/$(ARCH)
	-cp winSync.h $(ROOT)/include/winSync.h
	-@chmod g+rw $(ROOT)/libs/$(ARCH)/$(LIB)
	-@chmod g+rw $(ROOT)/include/winSync.h

install_dist: distdirs
	cvs tag -c $(DIST_TAG)
	-cp Debug/$(LIB) $(FULL_DIST_DIR)/libs/$(ARCH)
	-cp winSync.h $(FULL_DIST_DIR)/include/winSync.h

libdir:
	-@mkdir -p $(ROOT)/libs/$(ARCH)
	-@chmod g+rwxs $(ROOT)/libs/$(ARCH)

distdirs : 
	-@mkdir -p $(FULL_DIST_DIR)
	-@mkdir -p $(FULL_DIST_DIR)/libs/$(ARCH)
