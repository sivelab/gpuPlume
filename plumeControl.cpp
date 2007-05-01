#include <math.h>

#include "plumeControl.h"
#include "glErrorUtil.h"

#ifdef WIN32
#include <windows.h>
#include <stdio.h>
#include <conio.h>
#include <tchar.h>

#endif

// //////////////////////////////////////
// BEGIN -----> QUIC PLUME FORTRAN REFERENCES
// //////////////////////////////////////

#ifdef USE_PLUME_DATA

extern "C"
{
  void readfiles_();
}

// Domain size stored in nx, ny, and nz
extern "C" int __datamodule__nx;
extern "C" int __datamodule__ny;
extern "C" int __datamodule__nz;

extern "C" double __datamodule__dx;
extern "C" double __datamodule__dy;
extern "C" double __datamodule__dz;

// UVW contains the wind field
extern "C" double* __datamodule__u;
extern "C" double* __datamodule__v;
extern "C" double* __datamodule__w;

extern "C" int __datamodule__inumbuild;   // integer number of buildings
extern "C" double* __datamodule__xfo;
extern "C" double* __datamodule__yfo;
extern "C" double* __datamodule__zfo; 
extern "C" double* __datamodule__ht;
extern "C" double* __datamodule__wti;
extern "C" double* __datamodule__lti; 

#endif
// //////////////////////////////////////
// END ----> QUIC PLUME FORTRAN REFERENCES
// //////////////////////////////////////


PlumeControl::PlumeControl(){
 
#ifdef USE_PLUME_DATA
  // Call the PLUME code to read in the data files.
  std::cout << "Reading data using PLUME code..." << std::endl;
  readfiles_();

  //QuicPlume data for the domain
  nx = __datamodule__nx; //domain in the x direction
  ny = __datamodule__ny; //domain in the y direction
  nz = __datamodule__nz; //domain in the z direction

  //nx = (__datamodule__nx - 1) * __datamodule__dx; //domain in the x direction
  //ny = (__datamodule__nz - 1) * __datamodule__dz; //domain in the y direction
  //nz = (__datamodule__ny - 1) * __datamodule__dy; //domain in the z direction

  std::cout << "QUIC PLUME domain size: " << nx << " (in X) by " 
	    << ny << " (in Y) by " << nz << " (in Z)" << std::endl;

  //QuicPlume data for the windfield
  u = __datamodule__u;
  v = __datamodule__v;
  w = __datamodule__w;

  //QuicPlume data for the buildings
  numBuild = __datamodule__inumbuild;
  xfo = __datamodule__xfo;
  yfo = __datamodule__yfo;
  zfo = __datamodule__zfo;
  ht = __datamodule__ht;
  wti = __datamodule__wti;
  lti = __datamodule__lti;
  
#else
  nx = 60;
  ny = 60;//140;
  nz = 20;
  u=0;
  v=0;
  w=0;
#endif
  
  utility = new Util(this);
  utility->readInput("Settings/input.txt");

  //Sets up the type of simulation to run
  sim = new Simulation(useRealTime,duration,&time_step);
  
  texType = GL_TEXTURE_RECTANGLE_ARB;
  int_format = GL_RGBA32F_ARB;
  int_format_init = GL_RGBA;

  totalNumPar = 0.0;
  pos_buffer = new GLfloat[ twidth * theight * 4 ];

  //CollectionBox Settings
  avgTime = averagingTime + startCBoxTime;
  //If we want to output concentrations only at end set this
  avgTime = endCBoxTime;
  cBoxes[0] = new CollectionBox(numBox_x,numBox_y,numBox_z,bounds,averagingTime);
  num_cBoxes = 1;
  
  firstTime = true;
  endCBox = false;
  output_CollectionBox = false;
  odd = true; 
  dump_contents = false;
  quitSimulation = false;

  stream = new StreamLine(twidth,theight,nx,ny,nz);
  
  //Set whether to reuse particles or not
  //If reuseParticles is set to false: fourth coordinate of particle is -1 if emitted, else 0
  //If reuseParticles is set to true: fourth coordinate is <= lifetime if emiited, else lifetime+1
  reuseParticles = false;
  frameCount = 0;
  if(reuseParticles)
    lifeTime = 30.0;
  else lifeTime = -1.0;

  for(int i = twidth*theight-1; i >= 0; i--)
    indices.push_back(i);
   
}

void PlumeControl::init(bool OSG){
  osgPlume = OSG;

  pc = new ParticleControl(texType, twidth,theight,nx,ny,nz,u,v,w);
  pc->setUstarAndSigmas(ustar);

  dc = new DisplayControl(nx,ny,nz, texType);  
  dc->initVars(numBuild,xfo,yfo,zfo,ht,wti,lti);

  if(osgPlume){
    vp = new GLint[4];
    mvm = new GLfloat[16];
    pm = new GLfloat[16];
    dc->osgPlume = true;
  }
  
  if(testcase == 3){   
    numOfPE = 1;
    dc->draw_buildings = true;
    //Creates a point emitter with position(10,10,10) , emitting 10 particles per second
    pe[0] = new PointEmitter(10.0,10.0,10.0, 10.0, &twidth, &theight, &indices, &emit_shader);
  }
  else{
    dc->draw_buildings = false;
    for(int i=0; i < numOfPE; i++){
      if(radius[i] == 0)
	pe[i] = new PointEmitter(xpos[i],ypos[i],zpos[i], 10.0, &twidth, &theight, &indices, &emit_shader);
      else
	pe[i] = new SphereEmitter(xpos[i],ypos[i],zpos[i], 60.0, radius[i], &twidth, &theight, &indices, &emit_shader);
    }
  }
  for(int i=0; i < numOfPE; i++){
    if(reuseParticles)
      pe[i]->setParticleReuse(&indicesInUse, lifeTime);

    pe[i]->emit = false;

    //Set up the ParticleEmitter method to release particles
    //Release particles per time step only if duration is defined and
    //there is a fixed time step.
    if(duration != 0 && !useRealTime){
      pe[i]->releasePerTimeStep = true;
    }
    else{
      pe[i]->releaseOne = false;
      pe[i]->releasePerTimeStep = false;
      pe[i]->releasePerSecond = true;
    }

    if(pe[i]->releasePerTimeStep){
      //set number of particles to emit = (number of particles/ total number of time steps);
      int num = (int)floor((double)(twidth*theight) / (duration/(double)time_step));
      pe[i]->setNumToEmit(num);
    }
  }
  glEnable(texType);
  glGenTextures(8, texid);
  /////////////////////////////
  //Textures used:
  //texid[0] and texid[1] are the double buffered position textures
  positions0 = texid[0];
  positions1 = texid[1];
  //texid[2] can be used to initialize the positions
  //texid[3] is the wind field texture
  windField = texid[3];
  //texid[4] is random values
  randomValues = texid[4];
  prime0 = texid[5];
  prime1 = texid[6];
  lambda = texid[7];

  /////////////////////////////
  setupTextures();

  //
  // set up vertex buffer
  // 
  glGenBuffersARB(1, &vertex_buffer);
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, vertex_buffer);
  glBufferDataARB(GL_ARRAY_BUFFER_ARB, twidth*theight*4*sizeof(GLfloat), 0, GL_STREAM_COPY);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  //Initialize FBO
  initFBO();

  //This shader is used to advect the particles using the windfield
  pc->setupAdvectShader(&time_step, &numInRow, lifeTime);

  //This shader is used to update the prime values
  pc->setupPrimeShader(&numInRow); //Included argument -- Balli(04/12/07)

  //This shader is used to emmit particles
  emit_shader.addShader("Shaders/emitParticle_vp.glsl", GLSLObject::VERTEX_SHADER);
  emit_shader.addShader("Shaders/emitParticle_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  emit_shader.createProgram();

  //Initializes the particle positions in the domain
  //We don't need to do this anymore since we now can initialize data values
  //past 0-1 range in a RGBA32F texture without a shader.
  //pc->initParticlePositions(fbo, texid[2]); 
  CheckErrorsGL("END of init");

  if(useRealTime){
    display_clock = new Timer(true);
    sim->init();
  } 

}
//Error checking variable.
//Delete once emit particles problem solved. 
int d = 0;

void PlumeControl::display(){ 
  if(osgPlume){
    glGetFloatv(GL_MODELVIEW_MATRIX,mvm);
    glGetFloatv(GL_PROJECTION_MATRIX,pm);
    glGetIntegerv(GL_VIEWPORT,vp);

    glViewport(vp[0],vp[1],vp[2],vp[3]);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(-1, 1, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
  }

  //Makes sure the display loop is run once first
  //before emitting particles.  There is an error
  //in starting to emit particles the first time.
  ///////////////////////
  if(d < 1){
    d++;
  }
  else{
    //emit = true;
  //////////////////////
  //Set start time the second display call.
    //Once error is fixed, we can call this the first time.
    //Also can get rid of all the checks on if(!firstTime).
    if(firstTime){
      if(!osgPlume)
	pe[0]->emit = true;
      sim->setStartTime();
      firstTime = false;
    }
  }
 
  glGetIntegerv(GL_DRAW_BUFFER, &draw_buffer);
  
  // bind the framebuffer object so we can render to the 2nd texture
  fbo->Bind();
  
  //update simulation information such as total time elapsed and current time step
  //if set to run for a total number of time steps and that value has been reached,
  //clean up anything needed and close the window
  if(!firstTime){
    if(sim->update(&time_step)){
      if(!osgPlume){
	quitSimulation = true;
      }
    }
  }
  /////////////////////////////////////////////////////////////
  //Reuse particles if their lifespan is up
  /////////////////////////////////////////////////////////////
  if(reuseParticles){
    particleReuse();
  }
  ////////////////////////////////////////////////////////////
  // Emit Particles
  ////////////////////////////////////////////////////////////
  for(int i = 0; i < numOfPE; i++){
    if(pe[i]->emit){    
      if(pe[i]->releasePerTimeStep){
	//Releases particles per time step.
	//This can only be used with a time duration > 0
	//and a fixed time step
	pe[i]->setPosTexID(positions0, positions1);
	totalNumPar += (double)pe[i]->EmitParticle(fbo,odd);
      }
      else if(pe[i]->releaseOne){
	//Release one particle every time key is pressed
	pe[i]->setNumToEmit(1);
	pe[i]->setPosTexID(positions0, positions1);
	totalNumPar += (double)pe[i]->EmitParticle(fbo,odd);
 
	stream->addNewStream(pe[i]);
	pe[i]->emit = false;
      }
      else if(pe[i]->releasePerSecond){
	//Release particle using the defined particles per second
	if(pe[i]->timeToEmit(time_step))
	  {
	    pe[i]->setPosTexID(positions0, positions1);
	    totalNumPar += (double)pe[i]->EmitParticle(fbo, odd); 
	  }
      }   
    }
  }
  ////////////////////////////////////////////////////////////
  // Collection Boxes
  ////////////////////////////////////////////////////////////
  if(!endCBox){
    
    if(sim->totalTime >= endCBoxTime){
      endCBox = true;
      if(endCBoxTime != 0)
	output_CollectionBox = true;
    }       
  }

  if((sim->totalTime >= startCBoxTime) && !endCBox && !firstTime){
    
    glReadPixels(0, 0, twidth, theight, GL_RGBA, GL_FLOAT, pos_buffer); 
    for(int i = 3; i <= (theight*twidth*4); i+=4){
      //If particle has been emitted
      if(pos_buffer[i] == -1){
	
	//Get the x,y,z position of the particle
	float x = pos_buffer[i-3];
	float y = pos_buffer[i-2];
	float z = pos_buffer[i-1];

	//Check to see if particle is inside a collection box
	for(int j = 0; j < num_cBoxes; j++){
	  //if a particle is in a box the concentration value is updated
	  //cBoxes[j]->calculateConc(x,y,z,time_step,totalNumPar);
	  cBoxes[j]->calcSimpleConc(x,y,z);
	}	
      }
    } 
    if(sim->totalTime >= avgTime){
      avgTime += averagingTime;
      output_CollectionBox = true;
    }

  }
  ////////////////////////////////////////////////////////////
  // Update Prime Values
  ////////////////////////////////////////////////////////////
  if(!firstTime)
    pc->updatePrime(fbo, odd,positions0, positions1, prime0, prime1, 
		    windField, randomValues, lambda,time_step);
 
  ////////////////////////////////////////////////////////////
  // Update Particle Positions 
  ////////////////////////////////////////////////////////////
  if(!firstTime)
    pc->advect(fbo, odd, randomValues, windField, positions0, positions1, 
	       prime0,prime1,time_step);

  ////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////
  // Get Position for Streams
  ////////////////////////////////////////////////////////////
  if(!firstTime && stream->doUpdate()){
  //if(!firstTime && releaseOne){  
    stream->updateStreamPos();
  }
  ///////////////////////////////////////////////////////////

  CheckErrorsGL("END : after 1st pass");
  
  //Switches the frame buffer and binding texture
  odd = !odd;

  // In some circumstances, we may want to dump the contents of
  // the FBO to a file.
  if (dump_contents)
     {
       pc->dumpContents();
       dump_contents = false;
     }
  if(output_CollectionBox)
     {
       for(int j = 0; j < num_cBoxes; j++){
	 cBoxes[j]->outputConc(output_file,sim->totalTime,sim->curr_timeStep);
       }
       output_CollectionBox = false;
     }

  // We only need to do PASS 2 (copy to VBO) and PASS 3 (visualize) if
  // we actually want to render to the screen.  Rendering to the
  // screen will make the simulation run more slowly. This feature is
  // mainly included to allow some idea of how much faster the
  // simulation can run if left to run on the GPU.
  if (show_particle_visuals)
    {
      
      // //////////////////////////////////////////////////////////////
      // PASS 2 - copy the contents of the 2nd texture (the new positions)
      // into the vertex buffer
      // //////////////////////////////////////////////////////////////
      
      glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB, vertex_buffer);
      glReadPixels(0, 0, twidth, theight, GL_RGBA, GL_FLOAT, 0);
      glBindBuffer(GL_PIXEL_PACK_BUFFER_ARB, 0);
      CheckErrorsGL("after glReadPixels");
      
      // Disable the framebuffer object
      FramebufferObject::Disable();
      glDrawBuffer(draw_buffer); // send it to the original buffer
      CheckErrorsGL("END : after 2nd pass");

      // //////////////////////////////////////////////////////////////
      // PASS 3 - draw the vertices; This represents the visualization
      // of the PLUME particle field.
      // //////////////////////////////////////////////////////////////

      // clear the color and depth buffer before drawing the scene, and
      // set the viewport to the window dimensions
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      if(osgPlume){
	glViewport(vp[0], vp[1], vp[2], vp[3]);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glMultMatrixf(pm);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();		
	glMultMatrixf(mvm);
      }
      else{
	glViewport(0, 0, glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60.0, glutGet(GLUT_WINDOW_WIDTH)/float(glutGet(GLUT_WINDOW_HEIGHT)), 1.0, 250.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();		
      }

      //plume->displayVisual(vertex_buffer);
      
      dc->drawVisuals(vertex_buffer, windField, numInRow, twidth, theight);
      glDisable(texType);
      stream->draw();
      if(!osgPlume){
	for(int i=0; i < numOfPE; i++){
	  pe[i]->Draw();
	}
      }
      
      // If we've chose to display the 3D particle domain, we need to
      // set the projection and modelview matrices back to what is
      // needed for the particle advection step
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      gluOrtho2D(-1, 1, -1, 1);
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
     
      glDisable(texType);
      CheckErrorsGL("END : visualization");
      
      // Finally, swap the front and back buffers to display the
      // particle field to the monitor
      if(!osgPlume)
	glutSwapBuffers();
    }

  if(quitSimulation){
    std::cout << "Simulation ended after " << sim->simDuration << " seconds."<< std::endl;
    std::cout << "Total number of particles used: " << totalNumPar << std::endl;
    glutDestroyWindow(winid);
  }

}

void PlumeControl::initFBO(void){
  // Get the Framebuffer Object ready
  fbo = new FramebufferObject();
  fbo->Bind();
      
  rb = new Renderbuffer();
  rb->Set(GL_DEPTH_COMPONENT24, twidth, theight);
  fbo->AttachRenderBuffer(GL_DEPTH_ATTACHMENT_EXT, rb->GetId() );

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
	CheckErrorsGL("BEGIN : Creating textures");

	int sz = 4;
	GLfloat *data = new GLfloat[ twidth * theight * sz];
  
	for (int j=0; j<theight; j++)
		for (int i=0; i<twidth; i++)
		{
			int idx = j*twidth*sz + i*sz;
	
			//
			// Generate random positions for the particles within the
			// domain.  Currently, the domain is positive.
			//
			// With floating point textures, we have to create the inital
			// values between 0 and 1 and then use an initial shader to
			// transform the normalized coordinates to the correct domain.
      
			data[idx] = Random::uniform();
			data[idx+1] = Random::uniform();
			data[idx+2] = Random::uniform();
			data[idx+3] = Random::uniform();
		}
	pc->createTexture(texid[2], int_format_init, twidth, theight, data);

	CheckErrorsGL("\tcreated texid[2]...");
		
	// Creates wind field data texture
	pc->initWindTex(windField, lambda, &numInRow, testcase);
	CheckErrorsGL("\tcreated texid[3], the wind field texture...");

	for (int j=0; j<theight; j++)
		for (int i=0; i<twidth; i++)
		{
			int idx = j*twidth*sz + i*sz;
			data[idx] = -10.0;
			data[idx+1] = -10.0;
			data[idx+2] = 10.0;
			data[idx+3] = lifeTime+1;
		}
  
	// create the base texture with inital vertex positions
	pc->createTexture(texid[0], int_format, twidth, theight, data);
	CheckErrorsGL("\tcreated texid[0], the position texture...");

	// create a second texture to double buffer the vertex positions
	pc->createTexture(texid[1], int_format, twidth, theight, data);
	CheckErrorsGL("\tcreated texid[1], the position texture (double buffer)...");


	std::vector<float> random_values;
	//These two textures are to store the prime values(previous and updated values)
	//We will need to initialize some data into prime0
	for (int j=0; j<theight; j++)
	  for (int i=0; i<twidth; i++)
	    {
	      int idx = j*twidth*sz + i*sz;

	      data[idx] = Random::normal();
	      data[idx+1] = Random::normal();
	      data[idx+2] = Random::normal();
	      data[idx+3] = 0.0;

	      // record the random values to compute a mean, std, and var
	      random_values.push_back( data[idx] );
	      random_values.push_back( data[idx + 1] );
	      random_values.push_back( data[idx + 2] );

	      // normalize
	      float mag = sqrt(data[idx]*data[idx] + data[idx+1]*data[idx+1] 
			       + data[idx+2]*data[idx+2]);
	      //calculating u',v' and w'and storing them in data
	      //Quicplume coordinates
	      data[idx] = sigU*(data[idx]/mag);   
	      data[idx+1] = sigV*(data[idx+1]/mag);
	      data[idx+2] = sigW*(data[idx+2]/mag);
	      //Need to switch to gpuPlume coordinates
	      //data[idx] = sigV*(data[idx]/mag);   
	      //data[idx+1] = sigW*(data[idx+1]/mag);
	      //data[idx+2] = sigU*(data[idx+2]/mag);
	    }

	// Sum random values to determine if they have mean of zero and variance of 1
	double mean = 0.0;
	// std::cout << "randvalsprime = [" << std::endl;
	for (unsigned int i=0; i<random_values.size(); i++)
	  {
	    // std::cout << random_values[i] << std::endl;
	    mean += random_values[i];
	  }
	// std::cout << "];" << std::endl;
	mean = mean / (double)random_values.size();

	double stddev = 0.0, variance = 0.0, tmp_sum = 0.0;
	for (unsigned int i=0; i<random_values.size(); i++)
	  {
	    tmp_sum += ((random_values[i] - mean) * (random_values[i] - mean));
	  }
	variance = tmp_sum / (double)random_values.size();
	stddev = sqrt( variance );
	std::cout << "Prime Textures Random Values: Mean = " << mean << ", Standard Deviation = " << stddev << ", Variance = " << variance << std::endl;

	pc->createTexture(prime0, int_format, twidth, theight, data);
	CheckErrorsGL("\tcreated texid[5], the initial prime value texture...");

	pc->createTexture(prime1, int_format, twidth, theight, data);

	//
	// create random texture for use with particle simulation and turbulence
	//
	random_values.clear();
	for (int j=0; j<theight; j++)
	  for (int i=0; i<twidth; i++)
	    {
		int idx = j*twidth*sz + i*sz;
	
		//
		// Generate random values should of normal distribution with zero mean and standard deviation of one.
		// Need to pull classes from sim_fast that handle this... 
		// For now, generate random values between -1 and 1.... shader subtracts 1.0
		//
		data[idx] = Random::normal();
		data[idx+1] = Random::normal();
		data[idx+2] = Random::normal();
		data[idx+3] = 0.0;
		
		random_values.push_back( data[idx] );
		random_values.push_back( data[idx + 1] );
		random_values.push_back( data[idx + 2] );
	    }

	// Sum random values to determine if they have mean of zero and variance of 1
	mean = 0.0;
	for (unsigned int i=0; i<random_values.size(); i++)
	  {
	    mean += random_values[i];
	  }
	mean = mean / (double)random_values.size();

	stddev = 0.0;
	variance = 0.0; 
	tmp_sum = 0.0;
	for (unsigned int i=0; i<random_values.size(); i++)
	  {
	    tmp_sum += ((random_values[i] - mean) * (random_values[i] - mean));
	  }
	variance = tmp_sum / (double)random_values.size();
	stddev = sqrt( variance );
	std::cout << "texid[4] Random Values: Mean = " << mean << ", Standard Deviation = " << stddev << ", Variance = " << variance << std::endl;

	pc->createTexture(texid[4], int_format, twidth, theight, data);
	CheckErrorsGL("\tcreated texid[4], the random number texture...");

  delete [] data;

  CheckErrorsGL("END : Creating textures");
}

void PlumeControl::particleReuse(){
    if(frameCount == 0)
      if(useRealTime)
	reuse_time[0] = display_clock->tic();//reuse_clock->tic();

    frameCount++;
    if(frameCount == 10){
      

      if(useRealTime)
	reuse_time[1] = display_clock->tic();//reuse_clock->tic();
      //Iterate list indicesInUse; add time difference; 
      //if over lifetime; remove from list in indicesInUse
      //add that index into list indices
      iter = indicesInUse.begin();
      bool exit = false;
      while(iter != indicesInUse.end() && !exit){
	pIndex &reIndex = *iter;
	if(useRealTime)
	  reIndex.time += display_clock->deltas(reuse_time[0],reuse_time[1]);
	else
	  reIndex.time += frameCount*time_step;
	     
	if(reIndex.time >= lifeTime){
	  totalNumPar -= 1.0;
	  indicesInUse.erase(iter);
	  indices.push_front(reIndex.id);
	  std::cout << reIndex.time << " " << reIndex.id << std::endl;
	  exit = true;
	}
	iter++;
      }
      frameCount = 0;
    }


}
