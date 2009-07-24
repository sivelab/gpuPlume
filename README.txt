GPUPLUME

GETTING STARTED
---------------

gpuPlume now reads QUIC files directly using the .proj file associated
with the QUIC project as the starting point.  The previous mechanism
of editing a specific gpuPlume settings file so that it is associated
with a QUIC project and then copying wind field, turbulence, and cell
type files is no longer supported.

Building the code:
------------------
On Linux, you will need to type the following:

  cd fileParser
  make
  cd ..
  make

On Mac OS X, there is a bug in how the shaders work on NVIDIA graphics
cards.  We are working on this.  However, if you wanted to build the
code on OS X, the following will work:

  cd fileParser
  make
  cd ..
  make -f Makefile.osx 

On Windows, you'll need to load up the Visual Studio solution (in the
directory gpuQUIPLUME) using either VS 2005 or VS 2008.  The file
provided with gpuPlume is for VS 2005, but the code has been built
with VS 2008.  We have also tried the code on Windows XP as well as
variants of Vista.

Loading up an Urban Environment:
--------------------------------
To load an environment with gpuPlume, you need to start gpuPlume from
the command line and provide a QUIC .proj file as the first argument.
For example, several QUIC project files come with the gpuPlume source
and are located in the WindDomains directory.  To load one of these
environments, TestCase4, for instance, you'd do the following:

  ./gpuPlume WindDomains/TestCase4/TestCase4.proj

Note that the files and directory structure associated with each of
the projects in the WindDomains directories is the actual file and
directory structure associated with the QUIC system.  These same files
can be loaded into the QUIC GUI Matlab toolkit.

As of this writing most of the information needed to run a simulation
in gpuPlume is loaded from the QUIC files.  However, a few parameters
have not yet been converted to read from the QUIC files.  These
parameters can be set in the gpuPlume settings file.  Things like
number of particles, collection/concentration box information, and
various simulation options are currently set in this file.  By default
the Settings/input.txt file will be loaded.  If that file does not
exist, a default set of parameters is hard-coded into the system.  

However, if you have project specific settings files that you would
like to use, you can place a 

   gpuPlumeSettings.txt

into the top level directory of your project.  For instance, if you
edited Settings/input.txt to be specific for TestCase4 and you'd like
for that to load whenever TestCase4 is loaded, then you can copy that
file thusly:

  cp Settings/input.txt WindDomains/TestCase4/gpuPlumeSettings.txt

When setting the number of particles, you need to set the twidth and
theight parameters.  They are multiplied together to determine the
total number of particles; this is an artifact of using the OpenGL
system to do the computation and will eventually go away so you can
simply enter the total number of particles.


For Windows users, you will need to build the code with Visual Studio
and run it from the command prompt.  A Visual Studio solution file is
located in the directory gpuQUICPLUME.  Once the solution is built,
copy the executable to the main gpuPlume directory so that it can
execute as described above.  For instance,

  copy gpuQUICPLUME/Debug/gpuQUICPLUME.exe .

Once that is done, you can execute gpuplume as described above.


Navigation System:
------------------

  Keyboard Keys:

    w - Moves forward along view direction
    s - Moves backward along view direction
    a - Strafes left
    d - Strafes right

  Mouse Movement: Click button and hold for movement

    left button - turns left and right
    right button - looks up and down
    middle button - moves height of view up and down


Visual Controls:

  Layers:
  
    - - toggles between control axes, x,y, and z  

    K - increases layer along current axis
    k - decreases layer along current axis
    l - toggles between showing wind field, tau11,tau22,tau33,and tau13	
    
    1 - draws contours
    2 - draws transparent layers
    3 - toggles axis aligned and rotational planes
    
    Rotation plane:
      
      X - increases pitch
      x - decreases pitch
      Y - increases yaw
      y - decreases yaw
      R - increases roll
      r - decreases roll

  IsoSurfaces:

    i - toggles drawing the isosurface
    u - toggles fill mode for isosurface

  Particle Colors:
    
    v - toggles coloring particles using advect terms
        and using the oponnent color scheme to color directions

Simulation Controls:

  c - clears pathlines
  g - prints out prime values
  m - prints out mean velocities
  f - prints out positions
  z - toggles pause mode
  spacebar - increases one time step when in pause mode
  t - toggles visuals on or off (visuals off is faster)    

  Emitter controls:
  
    + - incrementally selects emitter
    e - toggles emission of selected emitter
    o - switches to release one particle per key press 
      	(press e to release one particle)
    p - switches to release particles per second
  
    W - move the emitter in negative Y
    S - move the emitter in positive Y
    A - move the emitter in negative X
    D - move the emitter in positive X
