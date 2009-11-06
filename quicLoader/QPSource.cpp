#include "QPSource.h"

bool qpSource::readQUICFile(const std::string &filename)
{
  std::cout << "\tParsing: QP_source.inp" << std::endl;;

  std::string source_filepath = quicFilesPath + "QP_source.inp";
  std::ifstream sourceFile(source_filepath.c_str(), std::ifstream::in);
  if(!sourceFile.is_open())
    {
      std::cerr << "gpuPlume could not open :: " << source_filepath << "." << std::endl;
      exit(EXIT_FAILURE);
    }
		
  std::string line;
  std::stringstream ss(line, std::stringstream::in | std::stringstream::out);

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



  return true;
}

bool qpSource::writeQUICFile(const std::string &filename)
{
  std::ofstream qpfile;
  qpfile.open(filename.c_str());

  if (qpfile.is_open())
    {
      qpfile << "!QUIC 5.51" << std::endl;

      

      return true;
    }

  return true;
}

