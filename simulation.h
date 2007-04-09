#ifndef __SIMULATION_H__
#define __SIMULATION_H__

#include "Timer.h"

class Simulation{

 public:
  Simulation(bool,double,float*);

  void init();
  void setStartTime();
  bool update(float*);

  //holds value for total elapsed time in running collection boxes
  double totalTime;
  
  //Duration in time(seconds) to run simulation
  //Set to 0 to run until user decides to quit
  double simDuration;

  //Total number of time steps to run simulation
  double tts;
  //Keeps track of the current time step of simulation
  double curr_timeStep;

 private:
  bool useRealTime;
  //controls whether or not to run the simulation until
  //user quits or for a total number of time steps
  bool infiniteSim;

  Timer* sim_clock;
  Timer_t sim_time[2];
  Timer_t display_time[2];

};

#endif //__SIMULATION_H__
