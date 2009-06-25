EXEC	= gpuPlume

CC	= g++ -g -Wall

CSRC 	= gpuPlume.cpp framebufferObject.cpp renderbuffer.cpp GLSL.cpp glErrorUtil.cpp particleControl.cpp pointEmitter.cpp sphereEmitter.cpp particleEmitter.cpp displayControl.cpp plumeControl.cpp Timer.cpp collectionBox.cpp util.cpp simulation.cpp streamLine.cpp Random.cpp nonGaussianModel.cpp GaussianModel.cpp Gaussian_2shaders_Model.cpp ReflectionModel.cpp MultipleBuildingsModel.cpp  GeomTest.cpp PathLine.cpp Contour.cpp VisualPlane.cpp lineEmitter.cpp planeEmitter.cpp IsoSurface.cpp

COBJS   = $(CSRC:.cpp=.o)

# LIB=-lGLEW -framework GLUT -framework OpenGL -framework Foundation
# LIB_PATH=
# INCLUDE_PATH=

LIB=-lGLEW -lglut -lGL -L/home/cs/vr/users/willemsn/quicurbCUDA/fileParser/standard -lStandardFileParser
LIB_PATH=
INCLUDE_PATH=


FLAGS = $(INCLUDE_PATH) $(LIB_PATH) $(LIB)

$(EXEC): $(COBJS)
	$(CC) -o $(EXEC) $(COBJS) $(FLAGS) 

%.o : %.cpp
	$(CC) -c $(INCLUDE_PATH) $<

clean:
	rm -f $(EXEC) *.o
