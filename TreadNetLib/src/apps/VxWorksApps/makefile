#This is our specific makefile for automatically installing the apps
# if specifying target install (not for compiling - that must be done
# through the appropriate program)

ROOT = ../../../
include $(ROOT)/Make/MakeDefs

ARCH = intel_vx5.5

#variables for each app's directory, add .out to get executable
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

install: appdir $(CLIENT_APPS) $(SERVER_APPS)
	@echo " Finished installing VxWorks Test Apps "

$(CLIENT_APPS) : 
	-cp ClientProjects/$@/PENTIUM4gnu/$@.out $(ROOT)/bin/$(ARCH)
	-@chmod g+rwx $(ROOT)/bin/$(ARCH)/$@.out

$(SERVER_APPS) : 
	-cp ServerProjects/$@/PENTIUM4gnu/$@.out $(ROOT)/bin/$(ARCH)
	-@chmod g+rwx $(ROOT)/bin/$(ARCH)/$@.out

install_dist: distdirs
	cvs tag -c $(DIST_TAG)
	-cp ClientProjects/$(TCAPP)/PENTIUM4gnu/$(TCAPP).out $(FULL_DIST_DIR)/bin/$(ARCH)
	-cp ClientProjects/$(SCAPP)/PENTIUM4gnu/$(SCAPP).out $(FULL_DIST_DIR)/bin/$(ARCH)
	-cp ClientProjects/$(TCSSAPP)/PENTIUM4gnu/$(TCSSAPP).out $(FULL_DIST_DIR)/bin/$(ARCH)
	-cp ClientProjects/$(FTCAPP)/PENTIUM4gnu/$(FTCAPP).out $(FULL_DIST_DIR)/bin/$(ARCH)
	-cp ServerProjects/$(TSAPP)/PENTIUM4gnu/$(TSAPP).out $(FULL_DIST_DIR)/bin/$(ARCH)
	-cp ServerProjects/$(TSNPAPP)/PENTIUM4gnu/$(TSNPAPP).out $(FULL_DIST_DIR)/bin/$(ARCH)
	-cp ServerProjects/$(TSTSAPP)/PENTIUM4gnu/$(TSTSAPP).out $(FULL_DIST_DIR)/bin/$(ARCH)
	-cp ServerProjects/$(TSTSNPAPP)/PENTIUM4gnu/$(TSTSNPAPP).out $(FULL_DIST_DIR)/bin/$(ARCH)
	-cp ServerProjects/$(FTSAPP)/PENTIUM4gnu/$(FTSAPP).out $(FULL_DIST_DIR)/bin/$(ARCH)
	-@chmod g+rx $(FULL_DIST_DIR)/bin/$(ARCH)/$(TCAPP).out
	-@chmod g+rx $(FULL_DIST_DIR)/bin/$(ARCH)/$(SCAPP).out
	-@chmod g+rx $(FULL_DIST_DIR)/bin/$(ARCH)/$(TCSSAPP).out
	-@chmod g+rx $(FULL_DIST_DIR)/bin/$(ARCH)/$(FTCAPP).out
	-@chmod g+rx $(FULL_DIST_DIR)/bin/$(ARCH)/$(TSAPP).out
	-@chmod g+rx $(FULL_DIST_DIR)/bin/$(ARCH)/$(TSNPAPP).out
	-@chmod g+rx $(FULL_DIST_DIR)/bin/$(ARCH)/$(TSTSAPP).out
	-@chmod g+rx $(FULL_DIST_DIR)/bin/$(ARCH)/$(TSTSNPAPP).out
	-@chmod g+rx $(FULL_DIST_DIR)/bin/$(ARCH)/$(FTSAPP).out

appdir:
	-@mkdir -p $(ROOT)/bin/$(ARCH)
	-@chmod g+rwxs $(ROOT)/bin/$(ARCH)

distdirs : 
	-@mkdir -p $(FULL_DIST_DIR)
	-@mkdir -p $(FULL_DIST_DIR)/bin/$(ARCH)
	-@chmod g+rwxs $(FULL_DIST_DIR)/bin/$(ARCH)

