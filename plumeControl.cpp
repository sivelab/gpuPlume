#include <math.h>

#include "plumeControl.h"
#include "glErrorUtil.h"

#ifdef WIN32
#include <windows.h>
#include <stdio.h>
#include <conio.h>
#include <tchar.h>

#endif

PlumeControl::~PlumeControl(){}

void PlumeControl::init(bool){}

int PlumeControl::display(long int seedValue)  ///value being passed over to print within the concentration file 
{ 
  return 1;
}

void PlumeControl::setupEmitters()
{
  for(int i=0; i < util->numOfPE; i++)
    {
      if(util->petype[i] == 1)
	// POINT EMITTER
	pe[i] = new PointEmitter(util->xpos[i],util->ypos[i],util->zpos[i], 
				 util->rate[i], twidth, theight, &indices, &emit_shader,
				 &random_values,pc->sig,nxdx,nydy,nzdz);
      else if (util->petype[i] == 2)
	// LINE EMITTER
	pe[i] = new LineEmitter(util->xpos[i], util->ypos[i], util->zpos[i], 
				util->xpos_e[i], util->ypos_e[i], util->zpos_e[i], 
				util->rate[i], 
				twidth, theight, &indices, &emit_shader,
				&random_values, pc->sig, nxdx, nydy, nzdz);
      else if (util->petype[i] == 3)
	// SPHERE EMITTER
	// pe[i] = new SphereEmitter(util->xpos[i],util->ypos[i],util->zpos[i], 
	// util->rate[i], util->radius[i], twidth, theight, &indices, &emit_shader,
	// &random_values,pc->sig,nxdx,nydy,nzdz);

	pe[i] = new SphereEmitter(util->xpos[i],util->ypos[i],util->zpos[i], 
				  util->rate[i], util->radius[i], twidth, theight, &indices, &emit_shader,
				  &random_values,pc->sig,nxdx,nydy,nzdz,pc->alph1ij,pc->alph2ij,pc->alph3ij
                                  ,pc->bet1ij,pc->bet2ij,pc->bet3ij,pc->gam1ij,pc->gam2ij,pc->gam3ij);
      //Balli: Included rotation parameter in sphere emitter call, other source types should also be changed

      else if (util->petype[i] == 4)
	// PLANE EMITTER
	pe[i] = new PlaneEmitter(util->xpos[i],util->ypos[i],util->zpos[i], 
				 util->xpos_e[i], util->ypos_e[i], 
				 util->rate[i],
				 twidth, theight, &indices, &emit_shader,
				 &random_values,pc->sig,nxdx,nydy,nzdz);
    }
 
  for(int i=0; i < util->numOfPE; i++){
    //if(reuseParticles)
      //pe[i]->setParticleReuse(&indicesInUse, lifeTime);
    if(continuousParticleFlow)
	pe[i]->setContinuousParticleFlow(&indicesInUse, lifeTime);

    pe[i]->emit = false;
    //Set for different methods of emitting particles
    if(util->emit_method == 0){
      pe[i]->Punch_Hole = true;
    }
    else pe[i]->Punch_Hole = false;

    //Set up the ParticleEmitter method to release particles
    //Release particles per time step only if duration is defined and
    //there is a fixed time step.
    switch(util->releaseType){
    case 0:
      pe[i]->releaseType = perTimeStep;
      break;
    case 1:
      pe[i]->releaseType = perSecond;
      break;
    case 2:
      pe[i]->releaseType = instantaneous;
      break;
    default:
      std::cout << "Error in setting up particle release type" << std::endl;
    }

    // If the release type equates to a certain number of particles
    // over a fixed duration, we need to calculate the number to emit.
    if(pe[i]->releaseType == perTimeStep)
      {
	//set number of particles to emit = (number of particles/ total number of time steps);
	float number = (twidth*theight) / (util->duration/time_step);

	// Num to emit should be averaged across all emitters (for now
	// until we figure out a more controllable scheme)
	number = number / (float)util->numOfPE;

	// number should likely be rounded up
	// pe[i]->setNumToEmit( round(number) );
	pe[i]->setNumToEmit( ceilf(number) );

	std::cout << "Emitting " << ceilf(number) << " particles per time step for " 
		  << util->duration/time_step << " simulation steps." << std::endl;
    }
  }
}

void PlumeControl::initFBO(void){
  // Get the Framebuffer Object ready
  fbo = new FramebufferObject();
  fbo->Bind();
      
  //rb = new Renderbuffer();
  //rb->Set(GL_DEPTH_COMPONENT24, twidth, theight);
  //fbo->AttachRenderBuffer(GL_DEPTH_ATTACHMENT_EXT, rb->GetId() );

  //Attach textures to framebuffer object
  fbo->AttachTexture(GL_COLOR_ATTACHMENT0_EXT, texType, texid[0]);
  fbo->AttachTexture(GL_COLOR_ATTACHMENT1_EXT, texType, texid[1]);
  fbo->AttachTexture(GL_COLOR_ATTACHMENT2_EXT, texType, prime0);
  fbo->AttachTexture(GL_COLOR_ATTACHMENT3_EXT, texType, prime1);

  fbo->IsValid();
  FramebufferObject::Disable();
}

void PlumeControl::setupTextures()
{
}

void PlumeControl::particleReuse()
{
  if(frameCount == 0 && useRealTime)
    reuse_time[0] = display_clock->tic();//reuse_clock->tic();
  
  frameCount++;
  if(frameCount == 10)
    {
      if(useRealTime)
	reuse_time[1] = display_clock->tic();//reuse_clock->tic();

      // Iterate list indicesInUse; add time difference; if over
      // lifetime; remove from list in indicesInUse add that index
      // into list indices

      iter = indicesInUse.begin();
      bool exit = false;
      while(iter != indicesInUse.end() && !exit)
	{
	  pIndex &reIndex = *iter;
	  if(useRealTime)
	    reIndex.time += display_clock->deltas(reuse_time[0],reuse_time[1]);
	  else
	    reIndex.time += frameCount*time_step;
	     
	  if(reIndex.time >= lifeTime)
	    {
	      std::cout << "reIndex.time = " << reIndex.time << ", lifeTime = " << lifeTime << std::endl;
	      totalNumPar -= 1.0;

	      indicesInUse.erase(iter);
	      std::cout << "Size of indices before = " << indices.size() << std::endl;
	      indices.push_front(reIndex.id);
	      std::cout << "Size of indices = " << indices.size() << std::endl;
	      std::cout << reIndex.time << " " << reIndex.id << std::endl;
	      exit = true;
	    }
	  iter++;
	}
      frameCount = 0;
    }


}


void PlumeControl::swapPauseMode() {
  inPauseMode = !inPauseMode;
}


void PlumeControl::writeShadowMapToFile() {
  // Nothing is done here.
}
