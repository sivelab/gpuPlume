#include "QUSimparams.h"

bool quSimParams::readQUICFile(const std::string &filename)
{
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
		
  std::cout << "\tParsing QU_simparams.inp:" << filename << std::endl;
  lfp->study(filename);

  // Check for discovery and default if necessary.		
  nx = (lfp->recall(ie_nx)) ? ie_nx.value : 0 ;
  ny = (lfp->recall(ie_ny)) ? ie_ny.value : 0 ;
  nz = (lfp->recall(ie_nz)) ? ie_nz.value : 0 ; nz++;

  if(nx == 0 || ny == 0) {std::cerr << "Error::quicLoader::one or more dimensions is zero." << std::endl; exit(EXIT_FAILURE);}

  dx = (lfp->recall(fe_dx)) ? fe_dx.value : 1. ;
  dy = (lfp->recall(fe_dy)) ? fe_dy.value : 1. ;
  dz = (lfp->recall(fe_dz)) ? fe_dz.value : 1. ;
		
  start_time = (lfp->recall(fe_start_time))     ? fe_start_time.value     : 0. ;
  time_incr = (lfp->recall(fe_time_incr))      ? fe_time_incr.value      : 0. ;
  num_time_steps = (lfp->recall(ie_num_time_steps)) ? ie_num_time_steps.value : 1 ;
		
  roof_type   = (lfp->recall(ie_roof_type))   ? ie_roof_type.value   : 0 ;
  upwind_type = (lfp->recall(ie_upwind_type)) ? ie_upwind_type.value : 0 ;
  canyon_type = (lfp->recall(ie_canyon_type)) ? ie_canyon_type.value : 0 ;
  intersection_flag = (lfp->recall(be_intersection_flag)) ? be_intersection_flag.value : false ;
		
  max_iterations     = (lfp->recall(ie_max_iterations))     ? ie_max_iterations.value     : 10000 ;
  residual_reduction = (lfp->recall(ie_residual_reduction)) ? ie_residual_reduction.value :     3 ;
  diffusion_flag     = (lfp->recall(be_diffusion_flag))     ? be_diffusion_flag.value     : false ;
  diffusion_step     = (lfp->recall(ie_diffusion_step))     ? ie_diffusion_step.value     :     1 ;
		
  domain_rotation = (lfp->recall(fe_domain_rotation)) ? fe_domain_rotation.value : 0. ;
  utmx = (lfp->recall(ie_utmx)) ? ie_utmx.value : 0 ;
  utmy = (lfp->recall(ie_utmy)) ? ie_utmy.value : 0 ;
		
  utm_zone      = (lfp->recall(ie_utm_zone))      ? ie_utm_zone.value      :     0 ;
  quic_cfd_type = (lfp->recall(be_quic_cfd_type)) ? be_quic_cfd_type.value : false ;
  wake_type     = (lfp->recall(ie_wake_type))     ? ie_wake_type.value     :     0 ;
		
  delete lfp;
}

bool quSimParams::writeQUICFile(const std::string &filename)
{
  std::ofstream qufile;
  qufile.open(filename.c_str());
  qufile << "!QUIC 5.51" << std::endl;

  if (qufile.is_open())
    {
      qufile << nx << "\t\t\t!nx - Domain Length(X) Grid Cells" << std::endl;
      qufile << ny << "\t\t\t!ny - Domain Width(Y) Grid Cells" << std::endl;
      qufile << nz << "\t\t\t!nz - Domain Height(Z) Grid Cells" << std::endl;
      qufile << dx << "\t\t\t!dx (meters)" << std::endl;
      qufile << dy << "\t\t\t!dy (meters)" << std::endl;
      qufile << dz << "\t\t\t!dz (meters)" << std::endl;
      qufile << start_time << "\t\t\t!decimal start time (hr)" << std::endl;
      qufile << time_incr << "\t\t\t!time increment (hr)" << std::endl;
      qufile << num_time_steps << "\t\t\t!total time increments" << std::endl;
      qufile << roof_type << "\t\t\t!rooftop flag (0-none, 1-log profile, 2-vortex)" << std::endl;
      qufile << upwind_type << "\t\t\t!upwind cavity flag (0-none, 1-Rockle, 2-MVP, 3-HMVP)" << std::endl;
      qufile << canyon_type << "\t\t\t!street canyon flag (0-none, 1-Roeckle, 2-CPB, 3-exp. param. PKK, 4-Roeckle w/ Fackrel)" << std::endl;
      qufile << intersection_flag << "\t\t\t!street intersection flag (0-off, 1-on)" << std::endl;
      qufile << max_iterations << "\t\t\t!Maximum number of iterations" << std::endl;
      qufile << residual_reduction << "\t\t\t!Residual Reduction (Orders of Magnitude)" << std::endl;
      qufile << diffusion_flag << "\t\t\t!Use Diffusion Algorithm (1 = on)" << std::endl;
      qufile << diffusion_step << "\t\t\t!Number of Diffusion iterations" << std::endl;
      qufile << domain_rotation << "\t\t\t!Domain rotation relative to true north (cw = +)" << std::endl;
      qufile << utmx << "\t\t\t!UTMX of domain origin (m)" << std::endl;
      qufile << utmy << "\t\t\t!UTMY of domain origin (m)" << std::endl;
      qufile << utm_zone << "\t\t\t!UTM zone" << std::endl;
      qufile << quic_cfd_type << "\t\t\t!QUIC-CFD Flag" << std::endl;
      qufile << wake_type << "\t\t\t!wake flag (0-none, 1-Rockle, 2-Modified Rockle)" << std::endl;

      return true;
    }

  return false;
}