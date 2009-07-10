/* util.cpp
   
   Last Author: $Author$
   Last Date Changed in Repository: $Date$

 */
#ifndef WIN32
#include <unistd.h>
#endif

#ifdef WIN32
#include <windows.h>
#include <stdio.h>
#include <conio.h>
#include <tchar.h>
#endif

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "util.h"
#include "gpuPlume.h"

#include "legacyFileParser.h"

Util::Util(){
  num = 0;
  numb = 0;
  bounds = new float[6];
  output_id = "";
  reuse_particles = false;
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
		char *cwd = new char[MAX_PATH];
		getcwd(cwd, 1024);
		cwdStr = cwd;
#endif
      // next, extract base name from the file name, and attempt open the files we need...
      size_t lastSlashPos = file.find_last_of( "/" );
      size_t lastDotPos = file.find_last_of( "." );

      int prefixLength = lastDotPos - lastSlashPos - 1;

      std::string filePrefix, pathPrefix;
      pathPrefix = file.substr(0, lastSlashPos);
      filePrefix = file.substr(lastSlashPos+1, prefixLength);
      std::cout << "Path Prefix: " << pathPrefix << std::endl;
      std::cout << "File Prefix: " << filePrefix << std::endl;

      // attempt to get the path to the QU_* and QP_* files
#ifdef WIN32
	  std::string localQuicFilePath = cwdStr + "\\" + filePrefix + "_inner\\";
#else
      std::string localQuicFilePath = cwdStr + "/" + pathPrefix + "/" + filePrefix + "_inner/";
#endif
      std::cout << "QUIC Files Path: " << localQuicFilePath << std::endl;

      readQUICBaseFiles( localQuicFilePath );

      // We've successfully loaded the QUIC files, so we can move on and return true
      return true;
    }
  else
    {
      std::cout << "Input file is NOT a QUIC project file.  Attempting to use original settings file parser." << std::endl;
    }

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
  return true;
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
}

bool Util::readQUICBaseFiles( std::string& QUICFilesPath )
{
  float f1;
  float* b = new float[6];
  int s_type;
  float* s = new float[8];
  float* c = new float[3];
  std::string s1;

  // set some defaults for now...
  twidth = 1024;
  theight = 1024;

  reuse_particles = true;

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

  useRealTime = false;
  time_step = 0.01;   // should read this from QP_params.inp

  emit_method = 0;

  // eventually get these from plume...
  duration = 0.0;
  startCBoxTime = 0.0;
  endCBoxTime = 0.0;
  averagingTime = 0.0;

  bounds[0] = bounds[1] = bounds[2] = 0.0;
  bounds[3] = bounds[4] = bounds[5] = 1.0;

  numBox_x = 1;
  numBox_y = 1;
  numBox_z = 1;
  volumeBox();

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
