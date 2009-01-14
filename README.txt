GPUPLUME

GETTING STARTED
---------------
The primary file to edit to change the simulation parameters is 

  Settings/input.txt

This file is mostly self-explanatory, but is used to change the number
of particles, the size of the domain, changing the timestep, setting
emitters within the domain, and setting up collection boxes for
concentration calculations.  One of the parameters, "quicFilesPath" is
used to provide a path in which the various QU and or QP files are
located.  Be default, the simulation will look in Settings/ for these
files, but this can be changed by setting the quicFilesPath
appropriately.

When setting the number of particles, you need to set the twidth and
theight parameters.  They are multiplied together to determine the
total number of particles; this is an artifact of using the OpenGL
system to do the computation and will eventually go away so you can
simply enter the total number of particles.

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
