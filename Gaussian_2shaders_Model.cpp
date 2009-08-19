#include <math.h>
#include "Gaussian_2shaders_Model.h"
#include "glErrorUtil.h"

#ifdef WIN32
#include <windows.h>
#include <stdio.h>
#include <conio.h>
#include <tchar.h>

#endif

Gaussian_2shaders_Model::Gaussian_2shaders_Model(Util* u){
  util = u;

  //from util
  twidth = util->twidth;
  theight = util->theight;
  nx = util->nx;
  ny = util->ny;
  nz = util->nz;
  time_step = util->time_step;

  nxdx = (int)(nx*(1.0/util->dx));
  nydy = (int)(ny*(1.0/util->dy));
  nzdz = (int)(nz*(1.0/util->dz));

  //Sets up the type of simulation to run
  sim = new Simulation(util->useRealTime,util->duration,&time_step);
  
  texType = GL_TEXTURE_RECTANGLE_ARB;
  int_format = GL_RGBA32F_ARB;
  int_format_init = GL_RGBA;

  totalNumPar = 0.0;  

  //CollectionBox Settings
  cBoxes[0] = new CollectionBox(util);
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
Gaussian_2shaders_Model::~Gaussian_2shaders_Model(){}

void Gaussian_2shaders_Model::init(bool OSG){
  osgPlume = OSG;

  pc = new ParticleControl(texType, twidth,theight,nx,ny,nz,util->dx,util->dy,util->dz, util);
  pc->setUstarAndSigmas(util->ustar);

  dc = new DisplayControl(nx,ny,nz, texType, inPauseMode, util);
  dc->initVars(util->numBuild,util->numSides,util->xfo,util->yfo,util->zfo,util->ht,util->wti,util->lti,util->gamma);
  dc->draw_buildings = false;
  
  if(osgPlume){
    vp = new GLint[4];
    mvm = new GLfloat[16];
    pm = new GLfloat[16];
    dc->osgPlume = true;
  }

  
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
  /////////////////////////////
  setupTextures();

  setupEmitters();
  //
  // set up vertex buffer
  // 
  glGenBuffersARB(1, &vertex_buffer);
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, vertex_buffer);
  glBufferDataARB(GL_ARRAY_BUFFER_ARB, twidth*theight*4*sizeof(GLfloat), 0, GL_STREAM_COPY);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  //Initialize FBO
  initFBO();

  pc->setupPrimeShader();
  pc->setupAdvectShader(lifeTime);

  //This shader is used to emmit particles
  emit_shader.addShader("Shaders/emitParticle_vp.glsl", GLSLObject::VERTEX_SHADER);
  emit_shader.addShader("Shaders/emitParticle_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  emit_shader.createProgram();

  CheckErrorsGL("END of init");

  if(util->useRealTime){
    display_clock = new Timer(true);
    sim->init();
  } 
}

int Gaussian_2shaders_Model::display(){
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
      totalNumPar += (double)pe[i]->EmitParticle(odd,positions0,positions1,time_step,
						 prime0,prime1);
      if(pe[i]->releaseType == onePerKeyPress){
	stream->addNewStream(pe[i]);
      }
    }
  }
  ////////////////////////////////////////////////////////////
  // Collection Boxes
  ////////////////////////////////////////////////////////////

  if(!endCBox){
    output_CollectionBox = cBoxes[0]->findConc(sim,&endCBox,odd);
  }
  
  ////////////////////////////////////////////////////////////
  // Update Prime Values and Particle Positions
  ////////////////////////////////////////////////////////////
   
  pc->updatePrime(odd,positions0, positions1, prime0, prime1, 
		  windField, randomValues, lambda,time_step);
 
  
  pc->advect(odd, windField, positions0, positions1, 
	     prime0,prime1,time_step);  
 
  ////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////
  // Get Position for Streams
  ////////////////////////////////////////////////////////////
  if(stream->doUpdate()){
    stream->updateStreamPos();
  }
  ///////////////////////////////////////////////////////////

  // In some circumstances, we may want to dump the contents of
  // the FBO to a file.
  if (dump_contents)
     {
       pc->printPositions(odd);
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

  CheckErrorsGL("END : after 1st pass");
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
	glClearColor(util->bcolor[0],util->bcolor[1],util->bcolor[2],1.0);
      }

      dc->drawVisuals(vertex_buffer, windField, 0, numInRow, twidth, theight, 0, 0);
      stream->draw();
      dc->drawLayers(windField, numInRow, 1.0);

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
  return 1;

}

void Gaussian_2shaders_Model::setupTextures(){
  CheckErrorsGL("BEGIN : Creating textures");

  int sz = 4;
  GLfloat *data = new GLfloat[ twidth * theight * sz];
		
  // Creates wind field data texture
  pc->initWindTex(windField, &numInRow, util->windFieldData);
  CheckErrorsGL("\tcreated texid[3], the wind field texture...");

  //Creates lambda texture
  pc->initLambdaTex(lambda);
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


  //These two textures are to store the prime values(previous and updated values)
  //We will need to initialize some data into prime0

  int iterations = 0;
  double mean = 1.0;
  double stddev = 0.0, variance = 0.0, tmp_sum = 0.0;

  while( !(-0.01 < mean && mean < 0.01 && 0.90 < variance && variance < 1.01) && iterations < 50){
    iterations++;

    random_values.clear();

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

	  data[idx] = util->sigU*(data[idx]);   
	  data[idx+1] = util->sigV*(data[idx+1]);
	  data[idx+2] = util->sigW*(data[idx+2]);
	}

    // Sum random values to determine if they have mean of zero and variance of 1
    mean = 0.0;
    // std::cout << "randvalsprime = [" << std::endl;
    for (unsigned int i=0; i<random_values.size(); i++)
      {
	// std::cout << random_values[i] << std::endl;
	mean += random_values[i];
      }
    // std::cout << "];" << std::endl;
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

  }
  std::cout << "Prime Textures Random Values: Mean = " << mean << ", Standard Deviation = " << stddev << ", Variance = " << variance << std::endl;
  std::cout << "Number of iterations to get random values: " << iterations << std::endl;

  pc->createTexture(prime0, int_format, twidth, theight, data);
  CheckErrorsGL("\tcreated texid[5], the initial prime value texture...");

  pc->createTexture(prime1, int_format, twidth, theight, data);

  //
  // create random texture for use with particle simulation and turbulence
  //
 
  iterations = 0;
  mean = 1.0;
  variance = 0.0;

  while( !(-0.01 < mean && mean < 0.01 && 0.90 < variance && variance < 1.01) && iterations < 50){
    random_values.clear();

    iterations++;
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

  }

  std::cout << "texid[4] Random Values: Mean = " << mean << ", Standard Deviation = " << stddev << ", Variance = " << variance << std::endl;
  std::cout << "Number of iterations to get random values: " << iterations << std::endl;

  pc->createTexture(texid[4], int_format, twidth, theight, data);
  CheckErrorsGL("\tcreated texid[4], the random number texture...");

  delete [] data;

  CheckErrorsGL("END : Creating textures");
 
}
