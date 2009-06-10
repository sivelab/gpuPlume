#PROG =	quicplume_LINUX_64.exe
PROG =	quicplume_LINUX_32.exe
#PROG =	quicplume_MACI.exe

SRCS =	quicplume.f95 evaporation.f95 
OBJS =	quicplume.o evaporation.o  

LIBS =	

CC = cc
CFLAGS = -O
# compile with absoft

#FC = f95
#FFLAGS = -lib -f -s -N15 -8 -N27 -O2
#F95 = f95
#F95FLAGS = -s -N27 -O3 -OPT:Olimit=0 -m64 -mcmodel=medium #
#LDFLAGS = -m64

# Compile with PGI

#F95 = pgf95
#F95FLAGS = -fastsse

# compile with g95

FC = g95
FFLAGS = -O2 #-fdefault-double-8 -fdefault-real-8 -fdefault-integer-8 -ffree-form -g
F95 = g95
F95FLAGS = -O2 -ffree-form -g

# compile with intel fortran

#FC = ifort
#FFLAGS = -O2 -r8 -i8 -autodouble -xW -free -convert little_endian -i-static
#F95 = ifort
#F95FLAGS = -O2 -r8 -xW -free -convert little_endian -i-static

all: $(PROG)

$(PROG): $(OBJS)
	$(F95) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

DATE := /bin/date
TODAY := `$(DATE) +%d-%h-%y`
TARFILE := quic-urb.$(TODAY).tar
COMPRESS := gzip

tar:
	tar -cf $(TARFILE) Makefile $(SRCS)
	$(COMPRESS) -v $(TARFILE)
clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f95

.f95.o:
	$(F95) $(F95FLAGS) -c $<
