//

//
// Read QP_params and all QU_simparams
// and all variable in QP_buildout, and 
//

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

// #include "legacyFileParser.h"

Util::Util(){
  num = 0;
  numb = 0;
  bounds = new float[6];
  output_id = "";
  reuse_particles = false;
  fullscreen = false;
  ignoreSignal = false;
  offscreenRender = false;

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

  winWidth = 1;
  winHeight = 1;
  problemID = 0;
  problemInstanceID = 0;
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
  bcolor[0] = 0.2f;
  bcolor[1] = 0.2f;
  bcolor[2] = 0.2f;

  advectChoice = 4;  // #4 = one shader for the multiple buildings model
  windFieldData = 5; // #use value of 5 for reading in from QUIC-FILES

  emit_method = 0;

  ustar = 0.084;
  sigU = 2.0f*ustar;
  sigV = 2.0f*ustar;
  sigW = 1.3f*ustar;

  pauseMode = false;
  calculateMeanVel = false;
  updateParticleColors = true;

  num_contour_regions = 5;

  //
  // can get these from QUIC
  //
  quicFilesPath = QUICFilesPath;

  output_file = "/tmp/junkGPUPlume.txt";
  output_id = 42;

  // ///////////////////////////////////////////////////////////
  // 
  // Parse and Read QU_metparams.inp file.
  // ///////////////////////////////////////////////////////////
  quMetParamData.readQUICFile(quicFilesPath + "QU_metparams.inp");

  // ///////////////////////////////////////////////////////////
  // 
  // Parse and Read QU_simparams.inp file.
  // ///////////////////////////////////////////////////////////
  quSimParamData.readQUICFile(quicFilesPath + "QU_simparams.inp");

  // Check for discovery and default if necessary.		
  nx = quSimParamData.nx;
  ny = quSimParamData.ny;
  nz = quSimParamData.nz;

  dx = quSimParamData.dx;
  dy = quSimParamData.dy;
  dz = quSimParamData.dz;
		
  // ///////////////////////////////////////////////////////////
  // 
  // Parse and Read QP_buildout.inp file.
  // ///////////////////////////////////////////////////////////
  qpBuildoutData.readQUICFile(quicFilesPath + "QP_buildout.inp");

  // ///////////////////////////////////////////////////////////
  // 
  // Parse and Read QP_params.inp file.
  // ///////////////////////////////////////////////////////////
  qpParamData.readQUICFile(quicFilesPath + "QP_params.inp");

  // set the various components used by the previous incantations of
  // the software.  In particular, use the particle number here unless
  // overwritten by the command line... set collection box info,
  // particle reuse, time_step, etc...

  // for now, something basic...
  twidth = (int)sqrt(static_cast<float>(qpParamData.numParticles));
  theight = (int)sqrt(static_cast<float>(qpParamData.numParticles));
  std::cout << "\t\tRequested " << qpParamData.numParticles << " particles.  Actually using " << twidth * theight << " particles!" << std::endl;

  reuse_particles = qpParamData.particleRecyclingFlag;

  useRealTime = false;
  time_step = static_cast<float>(qpParamData.timeStep);

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



  // ///////////////////////////////////////////////////////////
  // 
  // Parse and Read QU_buildings.inp file.
  // ///////////////////////////////////////////////////////////
  quBuildingData.readQUICFile(quicFilesPath + "QU_buildings.inp");

  numBuild = quBuildingData.buildings.size();

  // keep the old structures... for now... and fill in with these values...
  xfo = new float[numBuild];
  yfo = new float[numBuild];
  zfo = new float[numBuild];
  ht = new float[numBuild];
  wti = new float[numBuild];
  lti = new float[numBuild];
  numSides = new int[numBuild];
  gamma = new float[numBuild];
		
  for (int i=0; i<numBuild; i++)
    {
      xfo[i] = quBuildingData.buildings[i].xfo;
      yfo[i] = quBuildingData.buildings[i].yfo;
      zfo[i] = quBuildingData.buildings[i].zfo;
      ht[i] = quBuildingData.buildings[i].height;
      wti[i] = quBuildingData.buildings[i].width;
      lti[i] = quBuildingData.buildings[i].length;
      gamma[i] = quBuildingData.buildings[i].gamma;

      switch(quBuildingData.buildings[i].type)
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
	    std::cerr << "I don't know what kind of building " << quBuildingData.buildings[i].type << " is." << std::endl;
	    break;
	}
    }


  // ///////////////////////////////////////////////////////////
  // 
  // Parse and Read QP_params.inp file.
  // ///////////////////////////////////////////////////////////
  qpSourceData.readQUICFile(quicFilesPath + "QP_source.inp");

  // make it work with the other stuff...
  //
  // Allocate space for the sources
  //
  numOfPE = qpSourceData.sources.size();

  petype = new int[numOfPE];
  xpos = new float[numOfPE];
  ypos = new float[numOfPE];
  zpos = new float[numOfPE];
  xpos_e = new float[numOfPE];
  ypos_e = new float[numOfPE];
  zpos_e = new float[numOfPE];
  radius = new float[numOfPE];
  rate = new float[numOfPE];

  // if the releaseType was previously set to 0, it means that the
  // simulation is set to run for a specific duration, likely for an
  // experiment.  As such, we will keep the releaseType at 0.
  // Otherwise, set it to something else.  Currently, we are not doing
  // a good job of dealing with the release types as they pertain to
  // QUIC.  -Pete

  // Also!!!! we are not dealing with release type well on a source by
  // source basis.

  std::cerr << "DEV NOTE: source by source release types are not supported!" << std::endl;

  if (releaseType != 0)
    {
      if (qpSourceData.sources[0].release == 1)       // IR
	releaseType = 2;
      else if (qpSourceData.sources[0].release == 2)  // CONTINUOUS
	releaseType = 1;
      else if (qpSourceData.sources[0].release == 3)  // DISCRETE CONTINUOUS
	releaseType = 1;
    }

  for (int sourceId = 0; sourceId < qpSourceData.sources.size(); sourceId++)
    {
      switch(qpSourceData.sources[sourceId].geometry)
	{
	  // these are all treated as a sphere...
	  case qpSource::SPHERICAL_SHELL:
	  case qpSource::CYLINDER:
	  case qpSource::EXPLOSIVE:
	  case qpSource::AREA:
	    petype[sourceId] = 3;
	    
	    xpos[sourceId] = qpSourceData.sources[sourceId].points[0].x;
	    ypos[sourceId] = qpSourceData.sources[sourceId].points[0].y;
	    zpos[sourceId] = qpSourceData.sources[sourceId].points[0].z;
	    radius[sourceId] = qpSourceData.sources[sourceId].radius;
	    
	    rate[sourceId] = 800.0;
	    break;

	  case qpSource::LINE:
	    petype[sourceId] = 2;
	    
	    xpos[sourceId] = qpSourceData.sources[sourceId].points[0].x;
	    ypos[sourceId] = qpSourceData.sources[sourceId].points[0].y;
	    zpos[sourceId] = qpSourceData.sources[sourceId].points[0].z;

	    xpos_e[sourceId] = qpSourceData.sources[sourceId].points[1].x;
	    ypos_e[sourceId] = qpSourceData.sources[sourceId].points[1].y;
	    zpos_e[sourceId] = qpSourceData.sources[sourceId].points[1].z;

	    radius[sourceId] = 0.0;
	    rate[sourceId] = 800.0;
	    break;

	  default:
	    // case 4: // explosive
	    // case 6: // moving point
	    // case 8: // submunitions
	    std::cout << "\t\tEmitter Type " << qpSourceData.sources[sourceId].geometry << " not yet supported." << std::endl;
	    exit(EXIT_FAILURE);
	    break;
	}
    }
  
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

