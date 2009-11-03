/* util.cpp
   
   Last Author: $Author$
   Last Date Changed in Repository: $Date$

 */
#ifdef WIN32
#include <windows.h>
#include <stdio.h>
#include <conio.h>
#include <tchar.h>
#endif

#ifndef WIN32
#include <unistd.h>
#endif

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <limits.h>

#include "util.h"
#include "gpuPlume.h"

#include "legacyFileParser.h"

Util::Util(){
  num = 0;
  numb = 0;
  bounds = new float[6];
  output_id = "";
  reuse_particles = false;
  fullscreen = false;

  // Set some default values so that if they don't 
  // appear in the configuation file the correct
  // defaults are used.
  network_mode = -1;
  viewing_mode = 0;
  treadport_view = 'c';
  static_treadport_frustum = 1;

  sun_azimuth = 0;
  sun_altitude = 0;

  onlyCalcShadows = false;
}

bool Util::isPathAbsolute(const std::string &filename)
{
  // An absolute path would have a / that begins the string.
  // Likewise, on Windows, it would have a C:\ style beginning.
  if (filename[0] == '/')
    return true;
  
  return false;
}

bool Util::isQUICProjFile(std::ifstream& inputStream)
{
  // Simply check to see if this looks like a QUIC proj file.  If it
  // does, then we will attempt to load the data for the simulation
  // from the files in the <PROJ_FILENAME_PREFIX>_inner directory.
  
  // should probably reset the stream just to make sure we're at the beginning...

  int likelihood = 0;

  char line[1024];
  std::string s1;
  while( !inputStream.eof() )
    {
      inputStream.getline(line, 1024);
      if( line[ strlen(line)] == '\n' ){
	line[ strlen(line)] = '\0';
      }

      // comments "# ...." are ok so pass over them... look only for real keywords
      if(readComment(line))
	continue;

      if(read1String(line, "Creator:", &s1)) {
	likelihood++;
      }

      else if(read1String(line, "Date:", &s1)) {
	likelihood++;
      }

      else if(read1String(line, "Notes:", &s1)) {
	likelihood++;
      }

      else if(read1String(line, "Roof", &s1)) {
	likelihood++;
      }

      else if(read1String(line, "Nested", &s1)) {
	likelihood++;
      }

      else if(read1String(line, "Wind", &s1)) {
	likelihood++;
      }
    }

  if (likelihood >= 6)
    {
      // this is quite likely to be a PROJ file, return true.
      return true;
    }
  else 
    {
      return false;
    }
}

bool Util::readInput(std::string file){

  // Is the path relative or absolute.
  hasAbsolutePath = isPathAbsolute(file);
  std::cout << "Absolute Path = " << hasAbsolutePath << std::endl;

  std::ifstream in;
  in.open(file.c_str(),std::ios::in);

  if(in == NULL) 
    {
      // bad input file
      return false;
    }

  // Check to see if the input file is a project file.  If it is, then
  // attempt to load all of the information from the standard QUIC
  // files.  Other parameters can currently be set to defaults and
  // eventually handled from the command line or from a preferences
  // file.
  bool isProjFile = isQUICProjFile( in );

  // Close the file, since we will have worked through it...
  in.close();

  if (isProjFile)
    {
      // extract the path to the file name
      // find the current working directory
      std::string cwdStr;
     
#ifdef WIN32
      const size_t bufferSz = MAX_PATH;
      TCHAR buffer[bufferSz];
      GetCurrentDirectory(bufferSz, buffer);
      cwdStr = std::string(buffer);
#else
      char *cwd = new char[PATH_MAX];
      getcwd(cwd, PATH_MAX);
      cwdStr = cwd;
#endif

      // next, extract base name from the file name, and attempt open the files we need...
      size_t lastSlashPos = file.find_last_of( "/" );
      size_t lastDotPos = file.find_last_of( "." );

      int prefixLength = lastDotPos - lastSlashPos - 1;

      std::string filePrefix, pathPrefix;
      pathPrefix = file.substr(0, lastSlashPos);
      filePrefix = file.substr(lastSlashPos+1, prefixLength);
      // std::cout << "Path Prefix: " << pathPrefix << std::endl;
      // std::cout << "File Prefix: " << filePrefix << std::endl;

      // attempt to get the path to the QU_* and QP_* files
      std::string localQuicFilePath = "";
      if (hasAbsolutePath == false)
	localQuicFilePath = cwdStr;

      std::string slash;
#ifdef WIN32
      slash = "\\";
      localQuicFilePath += slash + filePrefix + "_inner" + slash;
#else
      slash = "/";
      if (hasAbsolutePath)
	localQuicFilePath += pathPrefix + slash + filePrefix + "_inner" + slash;
      else
	localQuicFilePath += slash + pathPrefix + slash + filePrefix + "_inner" + slash;
#endif

      std::cout << "QUIC Files Path: " << localQuicFilePath << std::endl;

      readQUICBaseFiles( localQuicFilePath );

      // Now, load the base settings file which contains parameters
      // such as number of particles, and related, gpuPlume specific
      // items.  The settings file also contains some information that
      // will eventually get picked up by reading various other QUIC
      // QU_* and QP_* files.

      // Attempt to find the settings file in the
      // localQuicFilePath... if not, read the default file.
      bool foundSettingsFile = false;

      std::ifstream gpuPlumeSettings_in;
      std::string gpuPlumeSettings_filename;
#ifdef WIN32
      gpuPlumeSettings_filename = localQuicFilePath + "..\\" + "gpuPlumeSettings.txt";
#else
      gpuPlumeSettings_filename = localQuicFilePath + "../" + "gpuPlumeSettings.txt";
#endif
     
      gpuPlumeSettings_in.open(gpuPlumeSettings_filename.c_str(),std::ios::in);
      if (gpuPlumeSettings_in == NULL) 
	{
	  // bad input filename
#ifdef WIN32
	  gpuPlumeSettings_filename = "Settings\input.txt";
#else
	  gpuPlumeSettings_filename = "Settings/input.txt";
#endif
	  gpuPlumeSettings_in.open(gpuPlumeSettings_filename.c_str(),std::ios::in);
	  if (gpuPlumeSettings_in)
	    {
	      std::cout << "Opening \"" << gpuPlumeSettings_filename << "\" for parsing of additional gpuPlume settings." << std::endl;
	      foundSettingsFile = true;
	    }
	}
      else
	{
	  std::cout << "Opening \"" << gpuPlumeSettings_filename << "\" for parsing of additional gpuPlume settings." << std::endl;
	  foundSettingsFile = true;
	}

      // if there is a settings file, we'll attempt to override any
      // values set by default.  Otherwise, the default values
      // assigned when the QUIC files were parsed will be used.
      if (foundSettingsFile)
	{
	  // re-open the file for parsing
	  char line[1024];
	  while(  !gpuPlumeSettings_in.eof() )
	    {
	      gpuPlumeSettings_in.getline(line, 1024);
	      if( line[ strlen(line)] == '\n' ){
		line[ strlen(line)] = '\0';
	      }
	  
	      parseLine(line);
	    }
	  gpuPlumeSettings_in.close();
	}

      // We've successfully loaded the QUIC files, so we can move on and return true
      return true;
    }
  else
    {
      std::cout << "***************************************" << std::endl;
      std::cout << "Base input file is NOT a QUIC project file.  The code\n";
      std::cout << "has changed significantly and loading from our original\n";
      std::cout << "setup input is not encouraged since the building reflection\n";
      std::cout << "has changed significantly and is now based solely on the celltype\n";
      std::cout << "structure." << std::endl;
      std::cout << "***************************************" << std::endl;
      return false;
    }

#if 0
  // ++++++++++++++++++++++++++++++++++++++++++++
  // This was the original code for parsing our base input file... gpuPlume now loads off the 
  // QUIC files directly and only uses the input.txt settings file to pull in parameters specific
  // to gpuPlume.
  // ++++++++++++++++++++++++++++++++++++++++++++

  // re-open the file for parsing
  in.open(file.c_str(),std::ios::in);
  char line[1024];
  while(  !in.eof() )
  {
    in.getline(line, 1024);
    if( line[ strlen(line)] == '\n' ){
      line[ strlen(line)] = '\0';
    }
    
    parseLine(line);
  }

  in.close();
#endif
}
void Util::parseLine(char* line){

  float f1;
  float* b1 = new float[7];
  float* b = new float[6];
  int s_type;
  float* s = new float[8];
  float* c = new float[3];
  std::string s1;

  if(readComment(line))
    return;

  if(read1Float(line, "twidth", &f1)){
	 twidth = (int)f1;
  }
  if(read1Float(line, "theight", &f1)){
	 theight = (int)f1;
  }
  if(read1Float(line, "pwidth", &f1)){
	 pwidth = (int)f1;
  }
  if(read1Float(line, "pheight", &f1)){
	 pheight = (int)f1;
  }
  if(read1Float(line, "time_step", &f1)){
    time_step = f1;
  }
  if(read1Float(line, "nx", &f1)){
    nx = (int)f1;
  }
  if(read1Float(line, "ny", &f1)){
    ny = (int)f1;
  }
  if(read1Float(line, "nz", &f1)){
    nz = (int)f1;
  }
  if(read1Float(line, "dx", &f1)){
    dx = f1;
  }
  if(read1Float(line, "dy", &f1)){
    dy = f1;
  }
  if(read1Float(line, "dz", &f1)){
    dz = f1;
  }
  if(read1Float(line, "windFieldData", &f1)){
    windFieldData = (int)f1;
  }
  if(read1Float(line, "useRealTime", &f1)){
    if(f1 == 0)
      useRealTime = false;
    else
      useRealTime = true;
  }
  if(read1String(line, "output_file", &s1)){
    output_file = s1;
  }
  if(read1String(line, "output_id", &s1)){
    output_id = s1;
  }

  if(read1Float(line, "duration", &f1)){
    duration = (double)f1;
  }
  if(read1Float(line, "startCBoxTime", &f1)){
    startCBoxTime = (double)f1;
  }
  if(read1Float(line, "endCBoxTime", &f1)){
    endCBoxTime = (double)f1;
  }
  if(read1Float(line, "averagingTime", &f1)){
    averagingTime = (double)f1;
  }
  if(read6Float(line, "bounds", b)){
    bounds = b;
  }
  if(read1Float(line, "numBox_x", &f1)){
    numBox_x = (int)f1;
  }
  if(read1Float(line, "numBox_y", &f1)){
    numBox_y = (int)f1;
  }
  if(read1Float(line, "numBox_z", &f1)){
    numBox_z = (int)f1;
	volumeBox();
  }
  if(read1Float(line, "ustar", &f1)){
    ustar = f1;
    sigU = 2.0f*f1;
    sigV = 2.0f*f1;
    sigW = 1.3f*f1;
  }
  if(read1Float(line, "reuse_particles", &f1)){
    reuse_particles = false;
    if((int)f1 == 0)
      reuse_particles = false;
    else
      reuse_particles = true;
  }
  if(read1Float(line, "show_particle_visuals", &f1)){
    if(f1 == 0)
      show_particle_visuals = false;
    else
      show_particle_visuals = true;
  }
  if(read1Float(line, "show_collectionBox_visuals", &f1)){
    if(f1 == 0)
      show_collectionBox_visuals = false;
    else
      show_collectionBox_visuals = true;
  }
  if(read1Float(line, "advectChoice", &f1)){
    advectChoice = (int)f1;
  }
  if(read1Float(line, "num_of_sources", &f1)){
    numOfPE = (int)f1;
    petype = new int[numOfPE];
    xpos = new float[numOfPE];
    ypos = new float[numOfPE];
    zpos = new float[numOfPE];
    xpos_e = new float[numOfPE];
    ypos_e = new float[numOfPE];
    zpos_e = new float[numOfPE];
    radius = new float[numOfPE];
    rate = new float[numOfPE];
  }
  if(readSourceInfo(line, "source_info", s_type, s)){
    if (s_type == 1) // POINT
      {
	petype[num] = 1;
	xpos[num] = s[0];
	ypos[num] = s[1];
	zpos[num] = s[2];
	radius[num] = 0.0;
	rate[num] = s[3];
      }
    else if (s_type == 2) // LINE
      {
	petype[num] = 2;
	xpos[num] = s[0];
	ypos[num] = s[1];
	zpos[num] = s[2];
	xpos_e[num] = s[3];
	ypos_e[num] = s[4];
	zpos_e[num] = s[5];
	radius[num] = 0.0;
	rate[num] = s[6];
      }
    else if (s_type == 3) // SPHERE
      {
	petype[num] = 3;
	xpos[num] = s[0];
	ypos[num] = s[1];
	zpos[num] = s[2];
	radius[num] = s[3];
	rate[num] = s[4];
      }
    else if (s_type == 4) // PLANE
      {
	petype[num] = 4;
	xpos[num] = s[0];
	ypos[num] = s[1];
	zpos[num] = s[2];
	xpos_e[num] = s[3];  // use xpos_e to store width
	ypos_e[num] = s[4];  // use ypos_e to store height
	rate[num] = s[5];
      }
    num++;
  }
  if(read1Float(line, "release_type",&f1)){
    releaseType = (int)f1;
  }
  if(read1Float(line, "emit_method", &f1)){
    emit_method = (int)f1;
  }
  if(read1Float(line, "numBuild", &f1)){
    numBuild = (int)f1;
    //introducing number of sides of the building for drawing purpose
    numSides = new int[numBuild];
    xfo = new float[numBuild];
    yfo = new float[numBuild];
    zfo = new float[numBuild];
    ht = new float[numBuild];
    wti = new float[numBuild];
    lti = new float[numBuild];
  }
   if(read7Float(line, "build_param", b1)){
	//printf("\n %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n", b1[0], b1[1], b1[2], b1[3], b1[4], b1[5], b1[6]);
    numSides[numb] = b1[0];
    xfo[numb] = b1[1];
    yfo[numb] = b1[2];
    zfo[numb] = b1[3];
    ht[numb]  = b1[4];
    wti[numb] = b1[5];
    lti[numb] = b1[6];
    printf("\n %d %f %f %f %f %f %f \n",numSides[numb], xfo[numb], yfo[numb], zfo[numb], ht[numb], wti[numb], b1[6]);
    numb++;
  }
  if(read1String(line, "quicFilesPath", &s1)){
    quicFilesPath = s1;
  }
  if(read1Float(line, "pauseMode", &f1)){
    if(f1 == 0)
      pauseMode = true;
    else
      pauseMode = false;
  }
  if(read1Float(line, "calculateMeanVel", &f1)){
    if(f1 == 0)
      calculateMeanVel = false;
    else
      calculateMeanVel = true;
  }

  if(read1Float(line, "updateParticleColors", &f1)){
    if(f1 == 0)
      updateParticleColors = false;
    else
      updateParticleColors = true;
  }
  else
    // if there's no updateParticles line... make sure we turn it on
    updateParticleColors = true;    

  if(read3Float(line, "back_color", c)){
    bcolor[0] = c[0];
    bcolor[1] = c[1];
    bcolor[2] = c[2];
  }
  if(read1Float(line, "contour_regions", &f1)){
    num_contour_regions = (int)f1;
  }
  if(read1Float(line, "network_mode", &f1)) {
    network_mode = int(f1);
  }
  if(read1Float(line, "viewing_mode", &f1)) {
    viewing_mode = (int)f1;
  }
  if(read1String(line, "treadport_view", &s1)) {
    treadport_view = s1.c_str()[0];
  }
  if(read1Float(line, "static_treadport_frustum", &f1)) {
    static_treadport_frustum = (int)f1;
  }
  if(read1Float(line, "sun_azimuth", &f1)) {
    sun_azimuth = (float)f1;
  }
  if(read1Float(line, "sun_altitude", &f1)) {
    sun_altitude = (float)f1;
  }
}

bool Util::readQUICBaseFiles( std::string& QUICFilesPath )
{
  float f1;
  float* b = new float[6];
  int s_type;
  float* s = new float[8];
  float* c = new float[3];
  std::string s1;

  pwidth = 10;
  pheight = 10;

  show_particle_visuals = true;
  show_collectionBox_visuals = false;

  // background color
  bcolor[0] = 0.2;
  bcolor[1] = 0.2;
  bcolor[2] = 0.2;

  advectChoice = 4;  // #4 = one shader for the multiple buildings model
  windFieldData = 5; // #use value of 5 for reading in from QUIC-FILES

  emit_method = 0;

  ustar = 0.084;
  sigU = 2.0f*ustar;
  sigV = 2.0f*ustar;
  sigW = 1.3f*ustar;

  pauseMode = true;
  calculateMeanVel = false;
  updateParticleColors = true;

  num_contour_regions = 5;

  //
  // can get these from QUIC
  //
  quicFilesPath = QUICFilesPath;

  output_file = "/tmp/junkGPUPlume.txt";
  output_id = 42;

  // create the legacy file parser to parse the QU_simparams.inp file.
  legacyFileParser* lfp = new legacyFileParser();  

  intElement ie_nx = intElement("nx - Domain Length(X) Grid Cells");
  intElement ie_ny = intElement("ny - Domain Width(Y) Grid Cells");
  intElement ie_nz = intElement("nz - Domain Height(Z) Grid Cells");
  lfp->commit(ie_nx);
  lfp->commit(ie_ny);
  lfp->commit(ie_nz);

  floatElement fe_dx = floatElement("dx (meters)");
  floatElement fe_dy = floatElement("dy (meters)");
  floatElement fe_dz = floatElement("dz (meters)");
  lfp->commit(fe_dx);
  lfp->commit(fe_dy);
  lfp->commit(fe_dz);
		
  floatElement fe_start_time   = floatElement("decimal start time (hr)");
  floatElement fe_time_incr    = floatElement("time increment (hr)");
  intElement ie_num_time_steps = intElement("total time increments");
  lfp->commit(fe_start_time);
  lfp->commit(fe_time_incr);
  lfp->commit(ie_num_time_steps);
		
  intElement ie_roof_type   = intElement("rooftop flag (0-none, 1-log profile, 2-vortex)");
  intElement ie_upwind_type = intElement("upwind cavity flag (0-none, 1-Rockle, 2-MVP, 3-HMVP)");
  intElement ie_canyon_type = intElement("street canyon flag (0-none, 1-Roeckle, 2-CPB, 3-exp. param. PKK, 4-Roeckle w/ Fackrel)");
  boolElement be_intersection_flag = boolElement("street intersection flag (0-off, 1-on)");
  lfp->commit(ie_roof_type);
  lfp->commit(ie_upwind_type);
  lfp->commit(ie_canyon_type);
  lfp->commit(be_intersection_flag);
		
  intElement ie_max_iterations     = intElement("Maximum number of iterations");
  intElement ie_residual_reduction = intElement("Residual Reduction (Orders of Magnitude)");
  boolElement be_diffusion_flag    = boolElement("Use Diffusion Algorithm (1 = on)");
  intElement ie_diffusion_step     = intElement("Number of Diffusion iterations");
  lfp->commit(ie_max_iterations);
  lfp->commit(ie_residual_reduction);
  lfp->commit(be_diffusion_flag);
  lfp->commit(ie_diffusion_step);
		
  floatElement fe_domain_rotation = floatElement("Domain rotation relative to true north (cw = +)");
  intElement ie_utmx              = intElement("UTMX of domain origin (m)");
  intElement ie_utmy              = intElement("UTMY of domain origin (m)");
  lfp->commit(fe_domain_rotation);
  lfp->commit(ie_utmx);
  lfp->commit(ie_utmy);
		
  intElement ie_utm_zone      = intElement("UTM zone");
  intElement be_quic_cfd_type = intElement("QUIC-CFD Flag");
  intElement ie_wake_type     = intElement("wake flag (0-none, 1-Rockle, 2-Modified Rockle)");
  lfp->commit(ie_utm_zone);
  lfp->commit(be_quic_cfd_type);
  lfp->commit(ie_wake_type);
		
  std::cout << "\tParsing: " << "QU_simparams.inp" << std::endl;

  lfp->study(quicFilesPath + "QU_simparams.inp");

  // Check for discovery and default if necessary.		
  nx = (lfp->recall(ie_nx)) ? ie_nx.value : 0 ;
  ny = (lfp->recall(ie_ny)) ? ie_ny.value : 0 ;
  nz = (lfp->recall(ie_nz)) ? ie_nz.value : 0 ; nz++;

  if(nx == 0 || ny == 0) {std::cerr << "Error::urbSetup::one or more dimensions is zero." << std::endl; exit(EXIT_FAILURE);}

  dx = (lfp->recall(fe_dx)) ? fe_dx.value : 1. ;
  dy = (lfp->recall(fe_dy)) ? fe_dy.value : 1. ;
  dz = (lfp->recall(fe_dz)) ? fe_dz.value : 1. ;
		
  float start_time = (lfp->recall(fe_start_time))     ? fe_start_time.value     : 0. ;
  float QU_time_step = (lfp->recall(fe_time_incr))      ? fe_time_incr.value      : 0. ;
  int num_time_steps = (lfp->recall(ie_num_time_steps)) ? ie_num_time_steps.value : 1 ;
		
  int roof_type   = (lfp->recall(ie_roof_type))   ? ie_roof_type.value   : 0 ;
  int upwind_type = (lfp->recall(ie_upwind_type)) ? ie_upwind_type.value : 0 ;
  int canyon_type = (lfp->recall(ie_canyon_type)) ? ie_canyon_type.value : 0 ;
  int intersection_flag = (lfp->recall(be_intersection_flag)) ? be_intersection_flag.value : false ;
		
  int max_iterations     = (lfp->recall(ie_max_iterations))     ? ie_max_iterations.value     : 10000 ;
  int residual_reduction = (lfp->recall(ie_residual_reduction)) ? ie_residual_reduction.value :     3 ;
  int diffusion_flag     = (lfp->recall(be_diffusion_flag))     ? be_diffusion_flag.value     : false ;
  int diffusion_step     = (lfp->recall(ie_diffusion_step))     ? ie_diffusion_step.value     :     1 ;
		
  int domain_rotation = (lfp->recall(fe_domain_rotation)) ? fe_domain_rotation.value : 0. ;
  float utmx = (lfp->recall(ie_utmx)) ? ie_utmx.value : 0 ;
  float utmy = (lfp->recall(ie_utmy)) ? ie_utmy.value : 0 ;
		
  int utm_zone      = (lfp->recall(ie_utm_zone))      ? ie_utm_zone.value      :     0 ;
  int quic_cfd_type = (lfp->recall(be_quic_cfd_type)) ? be_quic_cfd_type.value : false ;
  int wake_type     = (lfp->recall(ie_wake_type))     ? ie_wake_type.value     :     0 ;
		
  delete lfp;



  // ///////////////////////////////////////////////////////////
  // 
  // create the legacy file parser to parse the QP_params.inp file.
  // 
  lfp = new legacyFileParser();  

  intElement ie_sourceTypeFlag = intElement("Source type flag (1 = Basic, 2 = Dense Gas, 3 = Distributed Particle Size, 4 = Explosive, 5 = ERAD source, 6 = Bio Slurry, 7 = 2-phase, 8 = Exfiltration)");
  boolElement ie_isiteflag = boolElement("normal QUIC (isitefl=0) or sensor siting (=1) mode");
  boolElement ie_iindoorflag = boolElement("indoor calculations turned off (=0) or turned on (=1)");
  intElement ie_inextgridflag = intElement("1 - inner grid, 2 - outer grid");
  floatElement ie_westernEdge = floatElement("Location of western edge of inner grid relative to outer grid (m)");
  floatElement ie_southernEdge = floatElement("Location of southern edge of inner relative to outer grid (m)");
  floatElement ie_z0 = floatElement("Wall Surface Roughness Length (m)");
  floatElement ie_rcl = floatElement("Reciprocal Monin-Obukhov length(1/m)");
  floatElement ie_boundaryLayerHeight = floatElement("Boundary Layer height (m)");
  boolElement ie_nonLocalMixing = boolElement("use 1 to enable non-local mixing");
  intElement ie_numParticles = intElement("number of particles released over entire simulation");
  intElement ie_particleDistFlag = intElement("Number of particle distribution flag (1 = by mass, 2 = by source)");
  boolElement ie_particleSplitFlag = boolElement("Particle splitting flag");
  boolElement ie_particleRecyclingFlag = boolElement("Particle recycling flag");
  intElement ie_partNumIncreaseFactor = intElement("Total particle number increase factor");
  intElement ie_numParticleSplit = intElement("Number of particles a particle is split into");
  floatElement ie_partSplittingDosage = floatElement("Particle splitting target dose (gs/m^3)");
  floatElement ie_taylorMicroscaleMin = floatElement("Enable Taylor microscale lower limit to sub-time steps");
  intElement ie_randomNumberSeed = intElement("Random number seed");
  floatElement ie_timeStep = floatElement("time step (s)");
  floatElement ie_duration = floatElement("duration (s)");
  floatElement ie_concAvgTime = floatElement("concentration averaging time (s)");
  floatElement ie_concStartTime = floatElement("starting time for concentration averaging (s)");
  floatElement ie_partOutputPeriod = floatElement("particle output period (s)");
  floatElement ie_nbx = floatElement("in x direction, # of collecting boxes (concentration grid cells) ");
  floatElement ie_nby = floatElement("in y direction, # of collecting boxes (concentration grid cells) ");
  floatElement ie_nbz = floatElement("in z direction, # of collecting boxes (concentration grid cells) ");
  floatElement ie_xbl = floatElement("lower limits for collecting boxes in x in meters");
  floatElement ie_xbu = floatElement("upper limits for collecting boxes in x direction in meters");
  floatElement ie_ybl = floatElement("lower limits for collecting boxes in y in meters");
  floatElement ie_ybu = floatElement("upper limits for collecting boxes in y direction in meters");
  floatElement ie_zbl = floatElement("lower limits for collecting boxes in z in meters");
  floatElement ie_zbu = floatElement("upper limits for collecting boxes in z direction in meters");

  lfp->commit(ie_sourceTypeFlag);
  lfp->commit(ie_isiteflag);
  lfp->commit(ie_iindoorflag);
  lfp->commit(ie_inextgridflag);
  lfp->commit(ie_westernEdge);
  lfp->commit(ie_southernEdge);
  lfp->commit(ie_z0);
  lfp->commit(ie_rcl);
  lfp->commit(ie_boundaryLayerHeight);
  lfp->commit(ie_nonLocalMixing);
  lfp->commit(ie_numParticles);
  lfp->commit(ie_particleDistFlag);
  lfp->commit(ie_particleSplitFlag);
  lfp->commit(ie_particleRecyclingFlag);
  lfp->commit(ie_partNumIncreaseFactor);
  lfp->commit(ie_numParticleSplit);
  lfp->commit(ie_partSplittingDosage);
  lfp->commit(ie_taylorMicroscaleMin);
  lfp->commit(ie_randomNumberSeed);
  lfp->commit(ie_timeStep);
  lfp->commit(ie_duration);
  lfp->commit(ie_concAvgTime);
  lfp->commit(ie_concStartTime);
  lfp->commit(ie_partOutputPeriod);
  lfp->commit(ie_nbx);
  lfp->commit(ie_nby);
  lfp->commit(ie_nbz);
  lfp->commit(ie_xbl);
  lfp->commit(ie_xbu);
  lfp->commit(ie_ybl);
  lfp->commit(ie_ybu);
  lfp->commit(ie_zbl);
  lfp->commit(ie_zbu);

  std::cout << "\tParsing: " << "QP_params.inp" << std::endl;
  lfp->study(quicFilesPath + "QP_params.inp");

  // Check for discovery and default if necessary.		
  if (lfp->recall(ie_sourceTypeFlag))
    {
      if (ie_sourceTypeFlag.value == 1) qpParamData.sourceFlag = qpParams::BASIC;
      else if (ie_sourceTypeFlag.value == 2) qpParamData.sourceFlag = qpParams::DENSEGAS;
      else if (ie_sourceTypeFlag.value == 3) qpParamData.sourceFlag = qpParams::DISTPARTSIZE;
      else if (ie_sourceTypeFlag.value == 4) qpParamData.sourceFlag = qpParams::EXPLOSIVE;
      else if (ie_sourceTypeFlag.value == 5) qpParamData.sourceFlag = qpParams::ERADSOURCE;
      else if (ie_sourceTypeFlag.value == 6) qpParamData.sourceFlag = qpParams::BIOSLURRY;
      else if (ie_sourceTypeFlag.value == 7) qpParamData.sourceFlag = qpParams::TWOPHASE;
      else if (ie_sourceTypeFlag.value == 8) qpParamData.sourceFlag = qpParams::EXFILTRATION;
      else 
	qpParamData.sourceFlag = qpParams::BASIC;
    }

  qpParamData.isiteflag = (lfp->recall(ie_isiteflag)) ? ie_isiteflag.value : 0;
  qpParamData.iindoorflag = (lfp->recall(ie_iindoorflag)) ? ie_iindoorflag.value : false;
  qpParamData.inextgridflag = (lfp->recall(ie_inextgridflag)) ? ie_inextgridflag.value : 1;
  qpParamData.westernEdge = (lfp->recall(ie_westernEdge)) ? ie_westernEdge.value : 0;
  qpParamData.southernEdge = (lfp->recall(ie_southernEdge)) ? ie_southernEdge.value : 0;
  qpParamData.z0 = (lfp->recall(ie_z0)) ? ie_z0.value : 0;
  qpParamData.rcl = (lfp->recall(ie_rcl)) ? ie_rcl.value : 0;
  qpParamData.boundaryLayerHeight = (lfp->recall(ie_boundaryLayerHeight)) ? ie_boundaryLayerHeight.value : 0;
  qpParamData.nonLocalMixing = (lfp->recall(ie_nonLocalMixing)) ? ie_nonLocalMixing.value : 0;
  qpParamData.numParticles = (lfp->recall(ie_numParticles)) ? ie_numParticles.value : 0;
  qpParamData.particleDistFlag = (lfp->recall(ie_particleDistFlag)) ? ie_particleDistFlag.value : 0;
  qpParamData.particleSplitFlag = (lfp->recall(ie_particleSplitFlag)) ? ie_particleSplitFlag.value : 0;
  qpParamData.particleRecyclingFlag = (lfp->recall(ie_particleRecyclingFlag)) ? ie_particleRecyclingFlag.value : 0;
  qpParamData.partNumIncreaseFactor = (lfp->recall(ie_partNumIncreaseFactor)) ? ie_partNumIncreaseFactor.value : 0;
  qpParamData.numParticleSplit = (lfp->recall(ie_numParticleSplit)) ? ie_numParticleSplit.value : 0;
  qpParamData.partSplittingDosage = (lfp->recall(ie_partSplittingDosage)) ? ie_partSplittingDosage.value : 0;
  qpParamData.taylorMicroscaleMin = (lfp->recall(ie_taylorMicroscaleMin)) ? ie_taylorMicroscaleMin.value : 0;
  qpParamData.randomNumberSeed = (lfp->recall(ie_randomNumberSeed)) ? ie_randomNumberSeed.value : 0;
  qpParamData.timeStep = (lfp->recall(ie_timeStep)) ? ie_timeStep.value : 0;
  qpParamData.duration = (lfp->recall(ie_duration)) ? ie_duration.value : 0;
  qpParamData.concAvgTime = (lfp->recall(ie_concAvgTime)) ? ie_concAvgTime.value : 0;
  qpParamData.concStartTime = (lfp->recall(ie_concStartTime)) ? ie_concStartTime.value : 0;
  qpParamData.partOutputPeriod = (lfp->recall(ie_partOutputPeriod)) ? ie_partOutputPeriod.value : 0;
  qpParamData.nbx = (lfp->recall(ie_nbx)) ? ie_nbx.value : 0;
  qpParamData.nby = (lfp->recall(ie_nby)) ? ie_nby.value : 0;
  qpParamData.nbz = (lfp->recall(ie_nbz)) ? ie_nbz.value : 0;
  qpParamData.xbl = (lfp->recall(ie_xbl)) ? ie_xbl.value : 0;
  qpParamData.xbu = (lfp->recall(ie_xbu)) ? ie_xbu.value : 0;
  qpParamData.ybl = (lfp->recall(ie_ybl)) ? ie_ybl.value : 0;
  qpParamData.ybu = (lfp->recall(ie_ybu)) ? ie_ybu.value : 0;
  qpParamData.zbl = (lfp->recall(ie_zbl)) ? ie_zbl.value : 0;
  qpParamData.zbu = (lfp->recall(ie_zbu)) ? ie_zbu.value : 0;

  delete lfp;

  // set the various components used by the previous incantations of
  // the software.  In particular, use the particle number here unless
  // overwritten by the command line... set collection box info,
  // particle reuse, time_step, etc...

  // for now, something basic...
  twidth = (int)sqrt(qpParamData.numParticles);
  theight = (int)sqrt(qpParamData.numParticles);
  std::cout << "Requested " << qpParamData.numParticles << " particles.  Actually using " << twidth * theight << " particles!" << std::endl;

  reuse_particles = true; // qpParamData.particleRecyclingFlag;

  useRealTime = false;
  time_step = qpParamData.timeStep;

  // if the duration is close to zero... don't capture concentrations...
  duration = qpParamData.duration;
  if (duration >= -0.00001 && duration <= 0.00001)
    releaseType = 1;
  else 
    releaseType = 0;
  startCBoxTime = qpParamData.concStartTime;
  averagingTime = qpParamData.concAvgTime;
  endCBoxTime = qpParamData.concStartTime + qpParamData.concAvgTime;

  bounds[0] = qpParamData.xbl;
  bounds[1] = qpParamData.ybl;
  bounds[2] = qpParamData.zbl;
  
  bounds[3] = qpParamData.xbu;
  bounds[4] = qpParamData.ybu;
  bounds[5] = qpParamData.zbu;

  numBox_x = qpParamData.nbx;
  numBox_y = qpParamData.nby;
  numBox_z = qpParamData.nbz;
  volumeBox();

  std::cout << "QP_Params file: " << std::endl;
  std::cout << "qpParamData.sourceFlag = " << qpParamData.sourceFlag << std::endl;
  std::cout << "qpParamData.isiteflag = " << qpParamData.isiteflag << std::endl;
  std::cout << "qpParamData.iindoorflag = " << qpParamData.iindoorflag << std::endl;
  std::cout << "qpParamData.inextgridflag = " << qpParamData.inextgridflag << std::endl;
  std::cout << "qpParamData.westernEdge = " << qpParamData.westernEdge << std::endl;
  std::cout << "qpParamData.southernEdge = " << qpParamData.southernEdge << std::endl;
  std::cout << "qpParamData.z0 = " << qpParamData.z0 << std::endl;
  std::cout << "qpParamData.rcl = " << qpParamData.rcl << std::endl;
  std::cout << "qpParamData.boundaryLayerHeight = " << qpParamData.boundaryLayerHeight << std::endl;
  std::cout << "qpParamData.nonLocalMixing = " << qpParamData.nonLocalMixing << std::endl;
  std::cout << "qpParamData.numParticles = " << qpParamData.numParticles << std::endl;
  std::cout << "qpParamData.particleDistFlag = " << qpParamData.particleDistFlag << std::endl;
  std::cout << "qpParamData.particleSplitFlag = " << qpParamData.particleSplitFlag << std::endl;
  std::cout << "qpParamData.particleRecyclingFlag = " << qpParamData.particleRecyclingFlag << std::endl;
  std::cout << "qpParamData.partNumIncreaseFactor = " << qpParamData.partNumIncreaseFactor << std::endl;
  std::cout << "qpParamData.numParticleSplit = " << qpParamData.numParticleSplit << std::endl;
  std::cout << "qpParamData.partSplittingDosage = " << qpParamData.partSplittingDosage << std::endl;
  std::cout << "qpParamData.taylorMicroscaleMin = " << qpParamData.taylorMicroscaleMin << std::endl;
  std::cout << "qpParamData.randomNumberSeed = " << qpParamData.randomNumberSeed << std::endl;
  std::cout << "qpParamData.timeStep = " << qpParamData.timeStep << std::endl;
  std::cout << "qpParamData.duration = " << qpParamData.duration << std::endl;
  std::cout << "qpParamData.concAvgTime = " << qpParamData.concAvgTime << std::endl;
  std::cout << "qpParamData.concStartTime = " << qpParamData.concStartTime << std::endl;
  std::cout << "qpParamData.partOutputPeriod = " << qpParamData.partOutputPeriod << std::endl;
  std::cout << "qpParamData.nbx = " << qpParamData.nbx << std::endl;
  std::cout << "qpParamData.nby = " << qpParamData.nby << std::endl;
  std::cout << "qpParamData.nbz = " << qpParamData.nbz << std::endl;
  std::cout << "qpParamData.xbl = " << qpParamData.xbl << std::endl;
  std::cout << "qpParamData.xbu = " << qpParamData.xbu << std::endl;
  std::cout << "qpParamData.ybl = " << qpParamData.ybl << std::endl;
  std::cout << "qpParamData.ybu = " << qpParamData.ybu << std::endl;
  std::cout << "qpParamData.zbl = " << qpParamData.zbl << std::endl;
  std::cout << "qpParamData.zbu = " << qpParamData.zbu << std::endl;

  std::cout << "\tParsing: " << "QU_buildings.inp" << std::endl;

  // It's special...a common format is needed.
  std::string bld_filepath = quicFilesPath + "QU_buildings.inp";
  std::ifstream bldFile(bld_filepath.c_str(), std::ifstream::in);
  if(!bldFile.is_open())
    {
      std::cerr << "urbSetup could not open :: " << bld_filepath << "." << std::endl;
      exit(EXIT_FAILURE);
    }
		
  std::string line;
  std::stringstream ss(line, std::stringstream::in | std::stringstream::out);

  // first thing in these files is now a comment 
  getline(bldFile, line);

  int x_subdomain_start, y_subdomain_start, x_subdomain_end, y_subdomain_end;
  float zo;

  // x subdomain (southwest corner)
  getline(bldFile, line);
  ss.str(line);
  ss >> x_subdomain_start;
		
  // y subdomain (southwest corner)
  getline(bldFile, line);
  ss.str(line);
  ss >> y_subdomain_start;

  // x subdomain (northeast corner)
  getline(bldFile, line);
  ss.str(line);
  ss >> x_subdomain_end;
		
  // y subdomain (northeast corner)
  getline(bldFile, line);
  ss.str(line);
  ss >> y_subdomain_end;
		
  // wall roughness
  getline(bldFile, line);
  ss.str(line);
  ss >> zo;
		
  // number of buildings
  getline(bldFile, line);
  ss.str(line);
  int numbuilds = 0;
  ss >> numbuilds;

  numBuild = numbuilds;
  xfo = new float[numBuild];
  yfo = new float[numBuild];
  zfo = new float[numBuild];
  ht = new float[numBuild];
  wti = new float[numBuild];
  lti = new float[numBuild];
  numSides = new int[numBuild];
  gamma = new float[numBuild];
		
  // building description !Bld #	Group	Type	Height	Width	Length	Xfo	Yfo	Zfo	Gamma	Attenuation	Values in grid cell units
  //						!1	1	1	10	48	49	37	63	0	0	0
  getline(bldFile, line);
		
  // buildings
  int num = 0;
  int group = 0;
  int type = 0;
		
  float gamma_degrees=0.0,attenuation = 0.0;
	    
  float h,w,l;
  float x,y,z;

  for(int i = 0; i < numbuilds; i++)
    {
      getline(bldFile, line);
      ss.str(line);
      ss >> num	>> group >> type;
      ss >> h >> w >> l >> x >> y >> z >> gamma_degrees >> attenuation;
      ss.clear();
      xfo[i] = x;
      yfo[i] = y;
      zfo[i] = z;      
      ht[i] = h;
      wti[i] = w;
      lti[i] = l;
      gamma[i] = gamma_degrees;			
      
      switch(type)
      {
	  case 1:
            numSides[i]=4;
	    break;
	    
	case 2:    // building::CYLINDICAL
            
            numSides[i]=1;
	    break;
	    
	case 3:         // building::PENTAGON:
           
            numSides[i]=5;
	    break;

	    //case building::VEGETATION:	b = new vegetation(); break;

	  default:
	    std::cerr << "I don't know what kind of building " << type << " is." << std::endl;
	    break;
	}
    }
		
  bldFile.close();

  //
  // Read the QP_source file
  //
  std::cout << "\tParsing: QP_source.inp" << std::endl;;

  std::string source_filepath = quicFilesPath + "QP_source.inp";
  std::ifstream sourceFile(source_filepath.c_str(), std::ifstream::in);
  if(!sourceFile.is_open())
    {
      std::cerr << "gpuPlume could not open :: " << source_filepath << "." << std::endl;
      exit(EXIT_FAILURE);
    }
		
  // first thing in these files is now a comment 
  getline(sourceFile, line);

  int numberOfSources, numberOfSourceNodes;

  // Number of sources
  getline(sourceFile, line);
  ss.str(line);
  ss >> numberOfSources;
  ss.clear();
		
  // Number of source nodes
  getline(sourceFile, line);
  ss.str(line);
  ss >> numberOfSourceNodes;
  ss.clear();

  //
  // Allocate space for the sources
  //
  numOfPE = numberOfSources;;
  petype = new int[numOfPE];
  xpos = new float[numOfPE];
  ypos = new float[numOfPE];
  zpos = new float[numOfPE];
  xpos_e = new float[numOfPE];
  ypos_e = new float[numOfPE];
  zpos_e = new float[numOfPE];
  radius = new float[numOfPE];
  rate = new float[numOfPE];

  // read over the remainder of the source file and pull out the respective parts
  for(int i = 0; i < numOfPE; i++)
    {
      // First line in the source info is a comment like this: !Start of source number 1
      getline(sourceFile, line);

      // next is source name, which we don't use yet...
      getline(sourceFile, line);

      // source strength units
      int strengthUnits = -1;
      getline(sourceFile, line);
      ss.str(line);
      ss >> strengthUnits;
      ss.clear();

      // source strength 
      int sourceStr;
      getline(sourceFile, line);
      ss.str(line);
      ss >> sourceStr;
      ss.clear();

      // source density
      int sourceDensity;
      getline(sourceFile, line);
      ss.str(line);
      ss >> sourceDensity;
      ss.clear();

      // release type
      int rType;
      getline(sourceFile, line);
      ss.str(line);
      ss >> rType;
      ss.clear();

      // Release Type: 1 for instantaneous
      //               2 for continuous
      //               3 for discrete continous
      // 
      // Need to relate these values to our values, which are, of
      // course, different.  Ugh.  Needs to be reworked.
      //
      if (rType == 1) // IR
	releaseType = 1;
      else if (rType == 2)
	releaseType = 1;
      else if (rType == 3)
	releaseType = 0;

      releaseType = 1;

      // source start time
      int sourceStart;
      getline(sourceFile, line);
      ss.str(line);
      ss >> sourceStart;
      ss.clear();

      // source duration
      int sourceDuration;
      getline(sourceFile, line);
      ss.str(line);
      ss >> sourceDuration;
      ss.clear();

      // source geometry
      int geomType;
      getline(sourceFile, line);
      ss.str(line);
      ss >> geomType;
      ss.clear();

      // Source geometry (1 = spherical shell, 2 = line, 3 = cylinder,
      // 4 = Explosive,5 = Area, 6 = Moving Point, 7 = spherical
      // volume, 8 = Submunitions)
      switch(geomType)
	{
	  case 1:  // spherical shell
	  case 7:  // spherical volume
	    // spherical shell
	    petype[i] = 3;
	    
	    // x coord of sphere
	    getline(sourceFile, line);  
	    ss.str(line);
	    ss >> xpos[i];
	    ss.clear();

	    // y coord of sphere
	    getline(sourceFile, line);  
	    ss.str(line);
	    ss >> ypos[i];
	    ss.clear();

	    // z coord of sphere
	    getline(sourceFile, line);  
	    ss.str(line);
	    ss >> zpos[i];
	    ss.clear();

	    // radius
	    getline(sourceFile, line);  
	    ss.str(line);
	    ss >> radius[i];
	    ss.clear();

	    rate[i] = 800.0;

	    // Adding sphere source
	    std::cout << "\t\tSphere Source: " << xpos[i] << ',' << ypos[i] << ',' << zpos[i] << std::endl;
	    break;
	    
	  case 2: // line
	    petype[i] = 2;
	    
	    // !Numnber of data points
	    int numPts;
	    getline(sourceFile, line);
	    ss.str(line);
	    ss >> numPts;
	    ss.clear();

	    // !x (m)   y (m)   z (m)
	    getline(sourceFile, line);

	    // for (nPts = 0; nPts < numPts; nPts++)
	    // {
	    getline(sourceFile, line);  
	    ss.str(line);
	    ss >> xpos[i] >> ypos[i] >> zpos[i];
	    ss.clear();

	    getline(sourceFile, line);  
	    ss.str(line);
	    ss >> xpos_e[i] >> ypos_e[i] >> zpos_e[i];
	    ss.clear();

	    radius[i] = 0.0;
	    rate[i] = 800.0;

	    // Adding line source
	    std::cout << "\t\tLine Source: " << xpos[i] << ',' << ypos[i] << ',' << zpos[i] << " <---> " << xpos_e[i] << ',' << ypos_e[i] << ',' << zpos_e[i] << std::endl;
	    break;

	  case 3: // cylinder
	    // don't suppot cylinder yet, so stick a sphere there...
	    petype[i] = 3;

	    // !x coord of center of cylinder base (m)
	    getline(sourceFile, line);
	    ss.str(line);
	    ss >> xpos[i];
	    ss.clear();

	    // !y coord of center of cylinder base (m)
	    getline(sourceFile, line);
	    ss.str(line);
	    ss >> ypos[i];
	    ss.clear();

	    // !z coord of cylinder base (m)
	    getline(sourceFile, line);
	    ss.str(line);
	    ss >> zpos[i];
	    ss.clear();

	    // !radius of cylinder (m)
	    getline(sourceFile, line);
	    ss.str(line);
	    ss >> radius[i];
	    ss.clear();

	    // !height of cylinder (m)
	    getline(sourceFile, line);

	    rate[i] = 800.0;

	    std::cout << "\t\tCylinder Source: not added as cylinder... but as sphere." << std::endl;
	    break;

	  case 5: // area
	    float a_xfo, a_yfo, a_zfo, a_w, a_h, a_l, a_rot;

	    // !Area source xfo (m)
	    getline(sourceFile, line);
	    ss.str(line);
	    ss >> a_xfo;
	    ss.clear();

	    // !Area source yfo (m)
	    getline(sourceFile, line);
	    ss.str(line);
	    ss >> a_yfo;
	    ss.clear();

	    // !Area source zfo (m)
	    getline(sourceFile, line);
	    ss.str(line);
	    ss >> a_zfo;
	    ss.clear();

	    // !Area source length (m)
	    getline(sourceFile, line);
	    ss.str(line);
	    ss >> a_l;
	    ss.clear();

	    // !Area source width (m)
	    getline(sourceFile, line);
	    ss.str(line);
	    ss >> a_w;
	    ss.clear();

	    // !Area source height (m)
	    getline(sourceFile, line);
	    ss.str(line);
	    ss >> a_h;
	    ss.clear();

	    // !Area source rotation angle (o)
	    getline(sourceFile, line);
	    ss.str(line);
	    ss >> a_rot;
	    ss.clear();

	    //
	    // don't suppot area yet, so stick a sphere there at a reasonable location...
	    //
	    petype[i] = 3;
	    
	    xpos[i] = a_xfo;
	    ypos[i] = a_yfo;
	    zpos[i] = a_zfo;
	    radius[i] = a_h;

	    rate[i] = 800.0;

	    std::cout << "\t\tArea Source: not added directly, but represented as sphere." << std::endl;
	    break;

	    // case 4: // explosive
	    // case 6: // moving point
	    // case 8: // submunitions
	  default:
	    std::cout << "\t\tEmitter Type " << geomType << " not yet supported." << std::endl;
	    exit(EXIT_FAILURE);
	    break;
	}

      // After this, we again have a comment
      getline(sourceFile, line);        
    }
  sourceFile.close();

  return true;
}


bool Util::readSourceInfo(char *line, std::string settingName, int &source_type, float *f)
{
	std::istringstream ist(line);

	std::string w, source_typename;

	ist >> w;  // in other words, "source_info"
	if(w == settingName){

	  // check the source type, which will determine the remaining arguments
	  ist >> source_typename;
	  if (source_typename == "point")
	    {
	      source_type = 1;
	      ist >> f[0];
	      ist >> f[1];
	      ist >> f[2];
	      ist >> f[3];
	    }
	  else if (source_typename == "line")
	    {
	      source_type = 2;
	      ist >> f[0];
	      ist >> f[1];
	      ist >> f[2];
	      ist >> f[3];
	      ist >> f[4];
	      ist >> f[5];
	      ist >> f[6];
	    }
	  else if (source_typename == "sphere")
	    {
	      source_type = 3;
	      ist >> f[0];
	      ist >> f[1];
	      ist >> f[2];
	      ist >> f[3];
	      ist >> f[4];
	    }
	  else if (source_typename == "plane")
	    {
	      source_type = 4;
	      ist >> f[0];
	      ist >> f[1];
	      ist >> f[2];
	      ist >> f[3];
	      ist >> f[4];
	      ist >> f[5];
	    }
	  else
	    {
	      std::cerr << "\n*********************\nUnknown source type in settings file!\n*********************" << std::endl;
	      return false;
	    }

	  return true;
	}

    return false;
}
bool Util::read3Float(char *line, std::string settingName, float *f)
{
	std::istringstream ist(line);

	std::string w;
	ist >> w;
	if(w == settingName){
		ist >> f[0];
		ist >> f[1];
		ist >> f[2];
		return true;
	}

    return false;
}
bool Util::read7Float(char *line, std::string settingName, float *f)
{
	std::istringstream ist(line);

	std::string w;
	ist >> w;
	if(w == settingName){
		ist >> f[0];
		ist >> f[1];
		ist >> f[2];
		ist >> f[3];
		ist >> f[4];
		ist >> f[5];
                ist >> f[6];
		return true;
	}
	
    return false;
}
bool Util::read6Float(char *line, std::string settingName, float *f)
{
	std::istringstream ist(line);

	std::string w;
	ist >> w;
	if(w == settingName){
		ist >> f[0];
		ist >> f[1];
		ist >> f[2];
		ist >> f[3];
		ist >> f[4];
		ist >> f[5];
		return true;
	}

    return false;
}

bool Util::read1Float(char *line, std::string settingName, float *f)
{
	std::istringstream ist(line);

	std::string w;
	ist >> w;
	if(w == settingName){
		ist >> *f;
		return true;
	}

    return false;
}
bool Util::readComment(const char *line)
{

    if(strlen(line)==0)
    {
        return true;
    }

    int i=0;
    while(line[i] == ' ')
        i++;

    if(line[i] == '#')
    {
        return true;
    }

    return false;

}
bool Util::read1String(const char *line, const char *settingName, std::string *s)
{
	std::istringstream ist(line);

	std::string w;
	ist >> w;
	if(w == settingName){
		ist >> *s;
		return true;
	}

    return false;
}

void Util::volumeBox(){
	double xBoxlen = (bounds[3]-bounds[0])/numBox_x;
	double yBoxlen = (bounds[4]-bounds[1])/numBox_y;
	double zBoxlen = (bounds[5]-bounds[2])/numBox_z;

	volume= xBoxlen*yBoxlen*zBoxlen;
}

