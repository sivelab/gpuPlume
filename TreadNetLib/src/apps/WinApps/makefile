#This is our specific makefile for automatically installing the apps
# if specifying target install (not for compiling - that must be done
# through the appropriate program)

ROOT = ../../../
include $(ROOT)/Make/MakeDefs

ARCH = intel_win

#variables for each app's directory, add .exe to get executable
TCAPP     = tc
SCAPP     = sc
TCSSAPP   = tc_single_slope
FTCAPP    = full_packet_tc
TSAPP     = ts
TSNPAPP   = ts_noprint
TSTSAPP   = two_servers
TSTSNPAPP = two_servers_np
FTSAPP    = full_packet_ts

CLIENT_APPS = $(TCAPP) $(SCAPP) $(TCSSAPP) $(FTCAPP)
SERVER_APPS = $(TSAPP) $(TSNPAPP) $(TSTSAPP) $(TSTSNPAPP) $(FTSAPP)

#supports only install and install_dist

default clean veryclean :

install : appdir $(CLIENT_APPS) $(SERVER_APPS)
	@echo " Finished installing Windows Test Apps "

$(CLIENT_APPS) : 
	-cp ClientProjects/$@/Debug/$@.exe $(ROOT)/bin/$(ARCH)
	-@chmod g+rwx $(ROOT)/bin/$(ARCH)/$@.exe

$(SERVER_APPS) : 
	-cp ServerProjects/$@/Debug/$@.exe $(ROOT)/bin/$(ARCH)
	-@chmod g+rwx $(ROOT)/bin/$(ARCH)/$@.exe

install_dist: distdirs
	cvs tag -c $(DIST_TAG)
	-cp ClientProjects/$(TCAPP)/Debug/$(TCAPP).exe $(FULL_DIST_DIR)/bin/$(ARCH)
	-cp ClientProjects/$(SCAPP)/Debug/$(SCAPP).exe $(FULL_DIST_DIR)/bin/$(ARCH)
	-cp ClientProjects/$(TCSSAPP)/Debug/$(TCSSAPP).exe $(FULL_DIST_DIR)/bin/$(ARCH)
	-cp ClientProjects/$(FTCAPP)/Debug/$(FTCAPP).exe $(FULL_DIST_DIR)/bin/$(ARCH)
	-cp ServerProjects/$(TSAPP)/Debug/$(TSAPP).exe $(FULL_DIST_DIR)/bin/$(ARCH)
	-cp ServerProjects/$(TSNPAPP)/Debug/$(TSNPAPP).exe $(FULL_DIST_DIR)/bin/$(ARCH)
	-cp ServerProjects/$(TSTSAPP)/Debug/$(TSTSAPP).exe $(FULL_DIST_DIR)/bin/$(ARCH)
	-cp ServerProjects/$(TSTSNPAPP)/Debug/$(TSTSNPAPP).exe $(FULL_DIST_DIR)/bin/$(ARCH)
	-cp ServerProjects/$(FTSAPP)/Debug/$(FTSAPP).exe $(FULL_DIST_DIR)/bin/$(ARCH)
	-@chmod g+rx $(FULL_DIST_DIR)/bin/$(ARCH)/$(TCAPP).exe
	-@chmod g+rx $(FULL_DIST_DIR)/bin/$(ARCH)/$(SCAPP).exe
	-@chmod g+rx $(FULL_DIST_DIR)/bin/$(ARCH)/$(TCSSAPP).exe
	-@chmod g+rx $(FULL_DIST_DIR)/bin/$(ARCH)/$(FTCAPP).exe
	-@chmod g+rx $(FULL_DIST_DIR)/bin/$(ARCH)/$(TSAPP).exe
	-@chmod g+rx $(FULL_DIST_DIR)/bin/$(ARCH)/$(TSNPAPP).exe
	-@chmod g+rx $(FULL_DIST_DIR)/bin/$(ARCH)/$(TSTSAPP).exe
	-@chmod g+rx $(FULL_DIST_DIR)/bin/$(ARCH)/$(TSTSNPAPP).exe
	-@chmod g+rx $(FULL_DIST_DIR)/bin/$(ARCH)/$(FTSAPP).exe

appdir:
	-@mkdir -p $(ROOT)/bin/$(ARCH)
	-@chmod g+rwxs $(ROOT)/bin/$(ARCH)

distdirs : 
	-@mkdir -p $(FULL_DIST_DIR)
	-@mkdir -p $(FULL_DIST_DIR)/bin/$(ARCH)
	-@chmod g+rwxs $(FULL_DIST_DIR)/bin/$(ARCH)

