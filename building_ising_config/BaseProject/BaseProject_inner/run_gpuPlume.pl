#!/usr/bin/perl

use Cwd;

my $curr_dir = getcwd;

@windangle = ('000','045','090','135','180','225','270','315');

$xfoStart = 78;
$xfoMax = 82;

$yfoStart = 78;
$yfoMax = 82;

# Need to make sure the following files and directories are copied
# into the local space to gpuPlume has the files it needs.
#
#   building.ppm
#   buildingRoof.ppm
#   case_to_numpolys
#   concrete.ppm
#   edge_connect_list
#   Shaders/
#   SkyBox/
#
# This copying is automated by this script.

$command = sprintf("XCOPY /Y /D \"%s\\..\\..\\..\\building*.ppm\" .", $curr_dir);
print "Issuing --> $command\n";
system($command);

$command = sprintf("XCOPY /Y /D \"%s\\..\\..\\..\\case_to_numpolys\" .", $curr_dir);
print "Issuing --> $command\n";
system($command);

$command = sprintf("XCOPY /Y /D \"%s\\..\\..\\..\\concrete.ppm\" .", $curr_dir);
print "Issuing --> $command\n";
system($command);

$command = sprintf("XCOPY /Y /D \"%s\\..\\..\\..\\edge_connect_list\" .", $curr_dir);
print "Issuing --> $command\n";
system($command);

$command = sprintf("XCOPY /Y /D /E /I \"%s\\..\\..\\..\\Shaders\" Shaders", $curr_dir);
print "Issuing --> $command\n";
system($command);

$command = sprintf("XCOPY /Y /D /E /I \"%s\\..\\..\\..\\SkyBox\" SkyBox", $curr_dir);
print "Issuing --> $command\n";
system($command);

# Copy the gpuPlume.exe file here.
$command = sprintf("XCOPY /Y /D /E /I \"%s\\..\\..\\..\\gpuQUICPLUME\\Debug\\gpuQUICPLUME.exe\" .", $curr_dir);
print "Issuing --> $command\n";
system($command);

#
# Once files are copied, execute gpuPlume on the various scenarios.
#
for ($xfo=$xfoStart; $xfo<=$xfoMax; $xfo++)
{
    for ($yfo=$yfoStart; $yfo<=$yfoMax; $yfo++)
    {
	for ($wa=0; $wa<8; $wa++)
	{
	    # The windows time command isn't useful and does not do
	    # what the Unix time command does (namely tell us how long
	    # the execution of the 2nd argument took.
	    #
	    # So, don't use the time command here.
	    $command = sprintf("gpuQUICPLUME oae_%d_%d_$windangle[$wa]/gpuplume_input.txt", $xfo, $yfo);
	    print "Issuing --> $command\n";
	    system($command);
	}
    }
}

# When all done, remove the files we initially copied

;


