#!/usr/bin/perl

@windangle = ('000');

# For One Building 
$xfoStart = 78;
$xfoMax = 78;

$yfoStart = 78;  # half width of building
$yfoMax = 78;    # half width of building

#
# Currently, handle increments of 1 meter
#
#for ($xfo=$xfoStart; $xfo<=$xfoMax; $xfo++)
{
    #for ($yfo=$yfoStart; $yfo<=$yfoMax; $yfo++)
    {
	#for ($wa=0; $wa<8; $wa++)
	{
	    #$command = sprintf("cd F:\VAIO\Documents\URBAN_OPTIMIZATION\OptimizationAirEnergy_Data_mod\");
	    #system("cd F:\\VAIO\\Documents\\URBAN_OPTIMIZATION\\OptimizationAirEnergy_Data_mod\\OptimizationAirEnergy_Data\\oae_78_78_000\\oae_78_78_000_inner\\");
	    `quicurb_WIN.exe`;
          #system("dir");
	    #$command = sprintf("dir");
          #system("dir");
	    #print "Issuing --> $command\n";
	    #system($command);
	    
	    #$command = sprintf("cd oae_78_78_000/oae_78_78_00_inner; F:/VAIO/Documents/URBAN_OPTIMIZATION/OptimizationAirEnergy_Data_mod/OptimizationAirEnergy_Data/quicplume_WIN.exe", $xfo, $yfo, $xfo, $yfo);

	    #print "Issuing --> $command\n";
	    #system($command);
	}
    }
}
;
