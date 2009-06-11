This document describes the steps for creating a domain for use with
gpuPlume.

Notes:
------------------

1. gpuPlume was initially written with support for QUICv4.6_04-26-2007.
We have not seriously tested it with other versions of QUIC.

2. The only files gpuPlume requires are

      QU_velocity.dat
      QU_celltype.dat
      QP_turbfield.dat

Note that these are the ASCII versions of these files.  To have
quicurb or quicplume output the ASCII files, you need to edit
QU_fileoptions.inp and QP_fileoptions.inp to select ASCII.  Also note
that since we require the turbulence field, QP_fileoptions.inp must be
edited to force the output of the turbulence field.  You can also set
these from the GUI through the options menu.  We've generally edited
the files once the files are created.

Steps:
---------------------

1. Create your domain with QUICv4.6_04-26-2007 or a related version.

When you create your domain, you need to note the domain size as this
will need to be entered into the gpuplume settings file.  When
gpuplume was built, we did not build in support to read the simulation
parameters or domain parameters from the QUIC file set.

Set your wind field parameters and run QUICURB.

This should have created the celltype and velocity ASCII file.

2. The next step is to run a very quick simulation of QUICPLUME.  The
reason for this is that we need QUICPLUME to generate the turbulence
field.  So, with that in mind, you only need to set the simulation to
run very minimally.  By that, we mean that releasing 2 particles for 2
seconds would be quite sufficient.  The full release of particles will
happen with gpuPlume and is controlled by the gpuPlume settings file.

After running QUICPLUME, you may exit the matlab QUIC GUI.

3. 
