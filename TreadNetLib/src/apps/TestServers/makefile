ALLTARGETS = default clean install veryclean install_dist

$(ALLTARGETS):
	gmake -f prog_makefile $@ CFILES='ts.cpp Server_Common.cpp'
	gmake -f prog_makefile $@ CFILES='ts_noprint.cpp Server_Common.cpp'
	gmake -f prog_makefile $@ CFILES='two_servers.cpp Server_Common.cpp'
	gmake -f prog_makefile $@ CFILES='two_servers_np.cpp Server_Common.cpp'
	gmake -f prog_makefile $@ CFILES='full_packet_ts.cpp Server_Common.cpp'

