ALLTARGETS = default clean install veryclean install_dist

$(ALLTARGETS):
	make -f prog_makefile $@ CFILES='sc.cpp Client_Common.cpp'
	make -f prog_makefile $@ CFILES='tc.cpp Client_Common.cpp'
	make -f prog_makefile $@ CFILES='tc_single_slope.cpp Client_Common.cpp'
	make -f prog_makefile $@ CFILES='full_packet_tc.cpp Client_Common.cpp'

