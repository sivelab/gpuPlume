EXEC	= gpuPlume

CC	= g++ -g -O2 -Wall -DNDEBUG
CSRC 	= gpuPlume.cpp framebufferObject.cpp renderbuffer.cpp GLSL.cpp glErrorUtil.cpp particleControl.cpp particleEmitter.cpp displayControl.cpp plumeControl.cpp Timer.cpp
COBJS   = $(CSRC:.cpp=.o)

PLUME_DIR = Modular_QUICPLUME
FSRCS =	$(PLUME_DIR)/DataModule.f90 $(PLUME_DIR)/ReadFiles.f90
FOBJS =  $(FSRCS:.f90=.o)
FOBJSS = DataModule.o ReadFiles.o

# LIB=-lGLEW -framework GLUT -framework OpenGL -framework Foundation
# LIB_PATH=
# INCLUDE_PATH=

LIB=-Wl,-rpath=/home/cs/vr/software/glew-1.3.4.dist/lib -lGLEW -lglut -lGL 
LIB_PATH=-L/home/cs/vr/software/glew-1.3.4.dist/lib
INCLUDE_PATH=-I/home/cs/vr/software/glew-1.3.4.dist/include

# compile with gfortran
FCLOC	= /usr
FC	= $(FCLOC)/bin/gfortran
FCFLAGS	= -O3 -fdefault-double-8 -fdefault-real-8 -fdefault-integer-8 -ffree-form -ftree-vectorize

FLAGS = $(INCLUDE_PATH) $(LIB_PATH) $(LIB)

$(EXEC): $(COBJS) $(FOBJS)
	$(CC) -o $(EXEC) $(COBJS) $(FOBJSS) $(FLAGS) -L$(FLOC)/lib -lgfortran

%.o : %.cpp
	$(CC) -c $(INCLUDE_PATH) $<

%.o : %.f90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f $(EXEC) *.o
