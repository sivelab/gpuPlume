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

/*PlumeControl::PlumeControl(Util* u){

  util = u;

  //from util
  twidth = util->twidth;
  theight = util->theight;
  nx = util->nx;
  ny = util->ny;
  nz = util->nz;
  advectChoice = util->advectChoice;
  time_step = util->time_step;

  //Sets up the type of simulation to run
  sim = new Simulation(util->useRealTime,util->duration,&time_step);
  
  texType = GL_TEXTURE_RECTANGLE_ARB;
  int_format = GL_RGBA32F_ARB;
  int_format_init = GL_RGBA;

  totalNumPar = 0.0;  
  pos_buffer = new GLfloat[ twidth * theight * 4 ];

  //CollectionBox Settings
  avgTime = util->averagingTime + util->startCBoxTime;
  //If we want to output concentrations only at end set this
  avgTime = util->endCBoxTime;
  cBoxes[0] = new CollectionBox(util->numBox_x,util->numBox_y,util->numBox_z,util->bounds,util->averagingTime);
  
  num_cBoxes = 1;
  
  firstTime = true;
  endCBox = false;
  output_CollectionBox = false;
  odd = true; 
  dump_contents = false;
  createImages = false;
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
   
    }*/

void PlumeControl::init(bool OSG){
  /*
  osgPlume = OSG;

  pc = new ParticleControl(texType, twidth,theight,nx,ny,nz,util->u,util->v,util->w);
  pc->setUstarAndSigmas(util->ustar);

  dc = new DisplayControl(nx,ny,nz, texType);  
  dc->initVars(util->numBuild,util->xfo,util->yfo,util->zfo,util->ht,util->wti,util->lti);
  buildParam = new double[6];
  if(util->numBuild == 0){
    dc->draw_buildings = false;
    buildParam[0] = 0;
    buildParam[1] = 0;
    buildParam[2] = 0;
    buildParam[3] = 0;
    buildParam[4] = 0;
    buildParam[5] = 0;
  }
  else{
    dc->draw_buildings = true;
    buildParam[0] = util->xfo[0];
    buildParam[1] = util->yfo[0];
    buildParam[2] = util->zfo[0];
    buildParam[3] = util->ht[0];
    buildParam[4] = util->wti[0];
    buildParam[5] = util->lti[0];
  }

  if(osgPlume){
    vp = new GLint[4];
    mvm = new GLfloat[16];
    pm = new GLfloat[16];
    dc->osgPlume = true;
  }
  
  setupEmitters();
  
  glEnable(texType);
  glGenTextures(10, texid);
  /////////////////////////////
  //Textures used:
  positions0 = texid[0];
  positions1 = texid[1];
  windField = texid[3];
  randomValues = texid[4];
  prime0 = texid[5];
  prime1 = texid[6];
  lambda = texid[7];
  tau_dz = texid[8];
  duvw_dz = texid[9];
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

  //Setup Shaders for advection
  if(advectChoice == 1){
    //This shader performs both the update prime and advection by
    //writing to two buffers.
    pc->setupPrime_and_AdvectShader(numInRow,lifeTime);
  }
  else if(advectChoice == 2){
    pc->setupNonGaussianShader(numInRow,lifeTime);
  }
  else{
  //This shader is used to advect the particles using the windfield
    pc->setupAdvectShader(numInRow, lifeTime);

  //This shader is used to update the prime values
    pc->setupPrimeShader(numInRow); //Included argument -- Balli(04/12/07)
  }

  //This shader is used to emmit particles
  emit_shader.addShader("Shaders/emitParticle_vp.glsl", GLSLObject::VERTEX_SHADER);
  emit_shader.addShader("Shaders/emitParticle_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  emit_shader.createProgram();

  //Initializes the particle positions in the domain
  //We don't need to do this anymore since we now can initialize data values
  //past 0-1 range in a RGBA32F texture without a shader.
  //pc->initParticlePositions(fbo, texid[2]); 
  CheckErrorsGL("END of init");

  if(util->useRealTime){
    display_clock = new Timer(true);
    sim->init();
    } 
  */
}
bool imagesDone = false;
int PlumeControl::display(){ 
  /*if(osgPlume){
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
  //Store simulation start time and turn on one particle emitter
  if(firstTime){
    if(!osgPlume)
      pe[0]->emit = true;
    sim->setStartTime(&time_step);
    firstTime = false;
  }

  glGetIntegerv(GL_DRAW_BUFFER, &draw_buffer);
  glEnable(texType);
  // bind the framebuffer object so we can render to the 2nd texture
  fbo->Bind();
  
  /////////////////////////////////////////////////////////////
  //Reuse particles if their lifespan is up
  /////////////////////////////////////////////////////////////
  if(reuseParticles){
    particleReuse();
  }
  ////////////////////////////////////////////////////////////
   // Emit Particles
  ////////////////////////////////////////////////////////////
  for(int i = 0; i < util->numOfPE; i++){
    if(pe[i]->emit){    
      if(pe[i]->releasePerTimeStep){
	//Releases particles per time step.
	//This can only be used with a time duration > 0
	//and a fixed time step
	pe[i]->setPosTexID(positions0, positions1);
	pe[i]->setVertices();
	totalNumPar += (double)pe[i]->EmitParticle(fbo,odd);
      }
      else if(pe[i]->releaseOne){
	//Release one particle every time key is pressed
	pe[i]->setNumToEmit(1);
	pe[i]->setVertices();
	pe[i]->setPosTexID(positions0, positions1);
	totalNumPar += (double)pe[i]->EmitParticle(fbo,odd);
 
	stream->addNewStream(pe[i]);
	pe[i]->emit = false;
      }
      else if(pe[i]->releasePerSecond){
	
	//Release particle using the defined particles per second
	if(pe[i]->timeToEmit(time_step)){	  
	    pe[i]->setVertices();
	    pe[i]->setPosTexID(positions0, positions1);
	    totalNumPar += (double)pe[i]->EmitParticle(fbo, odd); 
	}
      }
      else if(pe[i]->instantRelease){
	pe[i]->setNumToEmit(twidth*theight);
	pe[i]->setVertices();
	pe[i]->setPosTexID(positions0,positions1);
	totalNumPar += (double)pe[i]->EmitParticle(fbo, odd);
	pe[i]->emit = false;
      }
    }
  }
  ////////////////////////////////////////////////////////////
  // Collection Boxes
  ////////////////////////////////////////////////////////////
  if(!endCBox){
    
    if(sim->totalTime >= util->endCBoxTime){
      sim->curr_timeStep += 1.0;
      endCBox = true;
      if(util->endCBoxTime != 0)
	output_CollectionBox = true;
    }       
  }

  if((sim->totalTime >= util->startCBoxTime) && !endCBox){

    sim->curr_timeStep +=1.0;

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
      avgTime += util->averagingTime;
      output_CollectionBox = true;
    }

  }
  
  ////////////////////////////////////////////////////////////
  // Update Prime Values and Particle Positions
  ////////////////////////////////////////////////////////////
   
  if(advectChoice == 1){
    pc->updatePrimeAndAdvect(odd,windField,positions0,positions1,
  			   prime0,prime1,randomValues,lambda,time_step);
  }
  //NOTE: The framebuffer object is no longer passed as an argument.
  //I don't think it's needed.
  else if(advectChoice == 2){
    pc->nonGaussianAdvect(odd,windField,positions0,positions1,prime0,prime1,
			  randomValues,lambda,tau_dz,duvw_dz,time_step);
  }
  else{
    pc->updatePrime(odd,positions0, positions1, prime0, prime1, 
  		    windField, randomValues, lambda,time_step);
 
    pc->advect(odd, windField, positions0, positions1, 
  	       prime0,prime1,time_step);  
  }  
  
  if(sim->totalTime >= 10.0 && !imagesDone){
    createImages = true;
    imagesDone = true;
  }
  if(createImages){
    pc->createPosImages(odd);
    pc->createPrimeImages(odd);
    
    createImages = false;
    }
  ////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////
  // Get Position for Streams
  ////////////////////////////////////////////////////////////
  if(stream->doUpdate()){
    stream->updateStreamPos();
  }
  ///////////////////////////////////////////////////////////

  CheckErrorsGL("END : after 1st pass");
  
  // In some circumstances, we may want to dump the contents of
  // the FBO to a file.
  if (dump_contents)
     {
       glGetIntegerv(GL_READ_BUFFER, &readbuffer);

       if(odd)
	 glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
       else
	 glReadBuffer(GL_COLOR_ATTACHMENT1_EXT);

       std::cout << "Previous Positions" << std::endl;
       pc->dumpContents();

       if(odd)
	 glReadBuffer(GL_COLOR_ATTACHMENT1_EXT);
       else
	 glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);
       
       std::cout << "Updated Positions" << std::endl;
       pc->dumpContents();

       dump_contents = false;
     }
  if(output_CollectionBox)
     {
       for(int j = 0; j < num_cBoxes; j++){
		cBoxes[j]->outputConc(util->output_file,sim->totalTime,sim->curr_timeStep);
       }
       output_CollectionBox = false;
     }

  //Switches the frame buffer and binding texture
  odd = !odd;

  // We only need to do PASS 2 (copy to VBO) and PASS 3 (visualize) if
  // we actually want to render to the screen.  Rendering to the
  // screen will make the simulation run more slowly. This feature is
  // mainly included to allow some idea of how much faster the
  // simulation can run if left to run on the GPU.
  if (util->show_particle_visuals)
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
	//glClearColor(0.93,0.93,0.93,1.0);	
      }
      
      dc->drawVisuals(vertex_buffer, windField, numInRow, twidth, theight);
      stream->draw();
         

      dc->drawLayers(windField, numInRow);

      if(!osgPlume){
	for(int i=0; i < util->numOfPE; i++){
	  pe[i]->Draw();
	}

	if(util->show_collectionBox_visuals){
	  cBoxes[0]->sort(dc->eye_pos[0],dc->eye_pos[1],dc->eye_pos[2]);
	  cBoxes[0]->draw(sim->curr_timeStep);
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
    // glutDestroyWindow(winid);
    return 0;
  }
  */
  return 1;
  
}
void PlumeControl::setupEmitters(){
  for(int i=0; i < util->numOfPE; i++){
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
      pe[i] = new SphereEmitter(util->xpos[i],util->ypos[i],util->zpos[i], 
				util->rate[i], util->radius[i], twidth, theight, &indices, &emit_shader,
				&random_values,pc->sig,nxdx,nydy,nzdz);
  }
 
  for(int i=0; i < util->numOfPE; i++){
    if(reuseParticles)
      pe[i]->setParticleReuse(&indicesInUse, lifeTime);

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

    if(pe[i]->releaseType == perTimeStep){
      //set number of particles to emit = (number of particles/ total number of time steps);
      float number = floorf((twidth*theight) / (util->duration/time_step));
      // Num to emit should be averaged across all emitters (for now
      // until we figure out a more controllable scheme)
      number = number / (float)util->numOfPE;

      int num = (int)number;

      pe[i]->setNumToEmit(num);
      //std::cout << (twidth*theight)/(util->duration/time_step) << std::endl;
      std::cout << "Emitting " << num << " particles per time step" <<std::endl;
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
  /*
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
	pc->initWindTex(windField, &numInRow, util->windFieldData);
	CheckErrorsGL("\tcreated texid[3], the wind field texture...");

	if(advectChoice == 2){
	  //pc->initLambdaTex(lambda, numInRow);
	  pc->initLambda_and_TauTex(lambda, tau_dz, duvw_dz, numInRow);
	}
	else
	  pc->initLambdaTex(lambda, numInRow);
	CheckErrorsGL("\tcreated texid[7], the lambda texture...");


	for (int j=0; j<theight; j++)
		for (int i=0; i<twidth; i++)
		{
			int idx = j*twidth*sz + i*sz;
			data[idx] = 100.0;
			data[idx+1] = 100.0;
			data[idx+2] = 100.0;
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

	//idx lookup using the source position to determine
	//the sigma values at that point in space. Only accurate
	//if there's only one source and it's a point emitter.
	if(advectChoice == 2){
	
	  int xs = (int)util->xpos[0];
	  int ys = (int)util->ypos[0];
	  int zs = (int)util->zpos[0];

	  int p2idx = zs*ny*nx + ys*nx + xs;

	  util->sigU = pc->sig[p2idx].u;
	  util->sigV = pc->sig[p2idx].v;
	  util->sigW = pc->sig[p2idx].w;
	}

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
	      //float mag = sqrt(data[idx]*data[idx] + data[idx+1]*data[idx+1] --Commented as it is not required now --Balli(05/01/07)
		  //	       + data[idx+2]*data[idx+2]);
	      //calculating u',v' and w'and storing them in data
	      //Quicplume coordinates

	      //data[idx] = sigU*(data[idx]/mag);   
	      //data[idx+1] = sigV*(data[idx+1]/mag);
	      //data[idx+2] = sigW*(data[idx+2]/mag);

	      data[idx] = util->sigU*(data[idx]);   
	      data[idx+1] = util->sigV*(data[idx+1]);
	      data[idx+2] = util->sigW*(data[idx+2]);
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
  */
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
