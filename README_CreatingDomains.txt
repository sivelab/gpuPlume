This document describes the steps for creating a domain for use with
gpuPlume.

Notes:
------------------

1. gpuPlume has worked with any version of QUIC as of QUICv4.6_04-26-2007.

2. gpuPlume currently reads the ASCII versions of the following files:

      QU_velocity.dat
      QU_celltype.dat
      QP_turbfield.dat

To have quicurb or quicplume output the ASCII files, you need to check
the appropriate boxes from the Options menu in the QUIC GUI to enable
ASCII output and to also output the turbulence field.

Steps:
---------------------

1. Create your domain and place buildings with the city builder.

Set your wind field parameters and run QUICURB.

This should have created the celltype and velocity ASCII file.

2. The next step is to run a very quick simulation of QUICPLUME.  The
reason for this is that we currently need QUICPLUME to generate the
turbulence field.  So, with that in mind, you only need to set the
simulation to run very minimally.  By that, we mean that releasing 2
particles for 2 seconds would be quite sufficient.  The full release
of particles will happen with gpuPlume and is controlled by the
gpuPlume settings file.

After running QUICPLUME, you may exit the matlab QUIC GUI.

3. At this point, you should be able to run your simulation in
gpuPlume by providing the proj file as the first argument to the
gpuPlume executable.

   ./gpuPlume /path/to/projects/MyProject/MyProject.proj


