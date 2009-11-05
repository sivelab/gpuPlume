#ifndef __QUICDATA_H__
#define __QUICDATA_H__ 1

// classes to hold the QUIC data...

class quicDataFile
{
public:
  quicDataFile() {}
  virtual ~quicDataFile() {}
  
  virtual bool readQUICFile(const std::string &filename) = 0;
  virtual bool writeQUICFile(const std::string &filename) = 0;

protected:
private:
};

// //////////////////////////////////////////////////////////////////
// 
// Class for holding the QP_Params.inp file
// 
// //////////////////////////////////////////////////////////////////
class qpParams : public quicDataFile
{
public:
  enum SourceType {
    BASIC = 1,
    DENSEGAS = 2,
    DISTPARTSIZE = 3,
    EXPLOSIVE = 4,
    ERADSOURCE = 5,
    BIOSLURRY = 6,
    TWOPHASE = 7,
    EXFILTRATION = 8    
  };
  
  qpParams() {}
  ~qpParams() {}

  SourceType sourceFlag;  
  short isiteflag;   // !normal QUIC (isitefl=0) or sensor siting (=1) mode
  bool iindoorflag; // !indoor calculations turned off (=0) or turned on (=1)
  short inextgridflag;      // !1 - inner grid, 2 - outer grid
  float westernEdge;  // !Location of western edge of inner grid relative to outer grid (m)
  float southernEdge; // !Location of southern edge of inner relative to outer grid (m)
  float z0;  // wallSurfRoughness;    // !Wall Surface Roughness Length (m)
  float rcl;  // !Reciprocal Monin-Obukhov length(1/m)
  float boundaryLayerHeight;  // !Boundary Layer height (m)
  bool nonLocalMixing;        // !use 1 to enable non-local mixing
  int numParticles;          // !number of particles released over entire simulation
  short particleDistFlag;     // !Number of particle distribution flag (1 = by mass, 2 = by source)
  bool particleSplitFlag;     // !Particle splitting flag
  bool particleRecyclingFlag; // !Particle recycling flag
  int partNumIncreaseFactor;  // !Total particle number increase factor
  short numParticleSplit;     // !Number of particles a particle is split into
  double partSplittingDosage; // !Particle splitting target dose (gs/m^3)
  float taylorMicroscaleMin;  // !Enable Taylor microscale lower limit to sub-time steps
  int randomNumberSeed;  // !Random number seed
  double timeStep;       // !time step (s)
  double duration;       // !duration (s)
  double concAvgTime;   // !concentration averaging time (s)
  double concStartTime; // !starting time for concentration averaging (s)
  double partOutputPeriod; // !particle output period (s)
  float nbx;  // !in x direction, # of collecting boxes (concentration grid cells) 
  float nby;  // !in y direction, # of collecting boxes (concentration grid cells) 
  float nbz;  // !in z direction, # of collecting boxes (concentration grid cells) 
  float xbl;  // !lower limits for collecting boxes in x in meters
  float xbu;  // !upper limits for collecting boxes in x direction in meters
  float ybl;  // !lower limits for collecting boxes in y in meters
  float ybu;  // !upper limits for collecting boxes in y direction in meters
  float zbl;  // !lower limits for collecting boxes in z in meters
  float zbu;  // !upper limits for collecting boxes in z direction in meters

  readQUICFile(const std::string &filename);
  writeQUICFile(const std::string &filename);

private:
};

// //////////////////////////////////////////////////////////////////
// 
// Class for holding the QP_Params.inp file
// 
// //////////////////////////////////////////////////////////////////
class quSimParams
{
public:
};

#if 0
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

#endif



// //////////////////////////////////////////////////////////////////
// 
// Class for holding the QU_buildings.inp file
// 
// //////////////////////////////////////////////////////////////////

class quBuildings
{
public:
};

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

// //////////////////////////////////////////////////////////////////
// 
// Class for holding the QP_source.inp file
// 
// //////////////////////////////////////////////////////////////////

class qpSource
{
public:
};

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


#   QU_metparams.inp

#     BaseProject.info  


#     QU_fileoptions.inp


#     sensor1.inp

#endif // #ifndef __QUICDATA_H__
