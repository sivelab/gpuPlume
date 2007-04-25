#include "simulation.h"
#include <iostream>

Simulation::Simulation(bool time,double t,float* time_step){
  useRealTime = time;
  simDuration = t;
  totalTime = 0.0; 
  curr_timeStep = 0;

  //tts = floor(simDuration/(*time_step));
  
  if(simDuration == 0)
    infiniteSim = true;
  else infiniteSim = false;

}
void Simulation::init(){
  if(useRealTime){
    sim_clock = new Timer(true);
    display_time[0] = sim_clock->tic();
  }
}

void Simulation::setStartTime(){
  if(useRealTime){
    sim_time[0] = sim_clock->tic();
    display_time[0] = sim_clock->tic();
  }
}
bool Simulation::update(float* time_step){ 

  //increment the current time step of the simulation
  curr_timeStep += 1.0;

  if(useRealTime){
    display_time[1] = sim_clock->tic();
    sim_time[1] = display_time[1];

    //real-time time_step
    *time_step = sim_clock->deltas(display_time[0],display_time[1]);
    //Total time elapsed
    totalTime = sim_clock->deltas(sim_time[0],sim_time[1]);
   
    display_time[0] = sim_clock->tic();
  }
  else{
    totalTime += *time_step;
  }

  
  if(!infiniteSim){
    if(totalTime >= simDuration){
      return true;
    }
  }
  
  return false;
}
float Simulation::randVal(){
#ifdef WIN32
  // Rand functions
  return (float)(rand()/(float)RAND_MAX); 
#else
  return drand48(); 
#endif

}
