#include <math.h>
#include "MultipleBuildingsModel.h"
#include "glErrorUtil.h"


#ifdef WIN32
#include <windows.h>
#include <stdio.h>
#include <conio.h>
#include <tchar.h>

#endif

MultipleBuildingsModel::MultipleBuildingsModel(Util* u){
  util = u;
  //pwidth = util->pwidth;
  //pheight = util->pheight;
  //pathNum = 0;

  //from util
  twidth = util->twidth;
  theight = util->theight;
  nx = util->nx;
  ny = util->ny;
  nz = util->nz;
  time_step = util->time_step;

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
  print_MeanVel = false;
  createImages = false;
  quitSimulation = false;
	
  //stream = new StreamLine(twidth,theight,nx,ny,nz);
  
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
MultipleBuildingsModel::~MultipleBuildingsModel(){}

void MultipleBuildingsModel::init(bool OSG){
  osgPlume = OSG;
  
  pathLines = new PathLine(util->pwidth,util->pheight,texType);

  pc = new ParticleControl(texType, twidth,theight,nx,ny,nz);
  pc->setUstarAndSigmas(util->ustar);
  pc->setBuildingParameters(util->numBuild,util->xfo,util->yfo,util->zfo,util->ht,util->wti,util->lti);
  pc->setQuicFilesPath(util->quicFilesPath);

  dc = new DisplayControl(nx,ny,nz, texType);  
  dc->initVars(util->numBuild,util->xfo,util->yfo,util->zfo,util->ht,util->wti,util->lti);

  if(util->numBuild == 0){
    dc->draw_buildings = false;  
  }
  else{   
    dc->draw_buildings = true;
  }

  if(osgPlume){
    vp = new GLint[4];
    mvm = new GLfloat[16];
    pm = new GLfloat[16];
    dc->osgPlume = true;
  }

  setupEmitters();
  
  glEnable(texType);
  glGenTextures(18, texid);
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

  //Texture for mean velocities
  meanVel0 = texid[10];
  meanVel1 = texid[11];
  //Texture used to hold current velocity
  currVel = texid[12];
  //Texture used for building information
  buildings = texid[13];
  //Texture used for cell type information
  cellType = texid[14];
 
  //Texture to store path lines
  pathTex = texid[15];
  
  //Texture for tau11,tau22,tau33,and tau13
  tau = texid[16];

  /////////////////////////////
  setupTextures(); 
  /////////////////////////////

  //Create isocontours
  contours = new Contour(pc,util->num_contour_regions);

  //
  // set up vertex buffer
  // 
  glGenBuffersARB(2, vbo_buffer);
  vertex_buffer = vbo_buffer[0];
  color_buffer = vbo_buffer[1];

  glBindBufferARB(GL_ARRAY_BUFFER_ARB, vertex_buffer);
  glBufferDataARB(GL_ARRAY_BUFFER_ARB, twidth*theight*4*sizeof(GLfloat), 0, GL_STREAM_COPY);
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, color_buffer);
  glBufferDataARB(GL_ARRAY_BUFFER_ARB, twidth*theight*4*sizeof(GLfloat), 0, GL_STREAM_COPY);

  glBindBuffer(GL_ARRAY_BUFFER, 0);

  //Initialize FBO
  initFBO();

  if(util->windFieldData >= 5)
    pc->setupMultipleBuildingsShader(lifeTime,0);
  else
    pc->setupMultipleBuildingsShader(lifeTime,1);

  pc->setupMeanVel_shader();

  pc->setupCurrVel_shader();

  //This shader is used to emmit particles
  emit_shader.addShader("Shaders/emitParticle_vp.glsl", GLSLObject::VERTEX_SHADER);
  emit_shader.addShader("Shaders/emitParticle_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  emit_shader.createProgram();

  dc->setupTurbulenceShader(pc->tauMax, pc->tauMin, pc->tauLocalMax, pc->tauLocalMin);
  dc->setupWindFieldShader(pc->windMax, pc->windMin);

  CheckErrorsGL("END of init");

  if(util->useRealTime){
    display_clock = new Timer(true);
    sim->init();
  } 

  int maxtextures;
  glGetIntegerv(GL_MAX_TEXTURE_IMAGE_UNITS, (GLint*)&maxtextures);
  std::cout << "Max Textures: " << maxtextures << std::endl;

}

int MultipleBuildingsModel::display(){

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

  glGetIntegerv(GL_DRAW_BUFFER, &draw_buffer);
  glEnable(texType);
  // bind the framebuffer object so we can render to the 2nd texture
  fbo->Bind();

  //update simulation information such as total time elapsed and current time step
  //if set to run for a total number of time steps and that value has been reached,
  //clean up anything needed and close the window
  if(!paused || !inPauseMode){

    if(!firstTime){
      if(sim->update(&time_step)){
	if(!osgPlume)
	  quitSimulation = true;
      }
    }
    //Store simulation start time and turn on one particle emitter
    if(firstTime){
      if(!osgPlume)
	pe[0]->emit = true;
      sim->setStartTime(&time_step);
      firstTime = false;
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
    for(int i = 0; i < util->numOfPE; i++){
      if(pe[i]->emit){    
	totalNumPar += (double)pe[i]->EmitParticle(odd,positions0,positions1,time_step);
	if(pe[i]->releaseType == onePerKeyPress){
	  //stream->addNewStream(pe[i]);
	  pathLines->addNewPath(pe[i]);

	}
      }
    }
    ////////////////////////////////////////////////////////////
    // Update Path Lines
    ////////////////////////////////////////////////////////////
    //FramebufferObject::Disable();
    pathFbo->Bind();
    pathLines->updatePathLines(positions0,positions1,odd);
    
    //Make sure to bind this fbo, because the pathLines use a
    //different fbo
    fbo->Bind();

    ////////////////////////////////////////////////////////////
    // Collection Boxes
    ////////////////////////////////////////////////////////////

    if(!endCBox){
      output_CollectionBox = cBoxes[0]->findConc(sim,&endCBox,odd); 
    }

    ////////////////////////////////////////////////////////////
    // Update Prime Values and Particle Positions
    ////////////////////////////////////////////////////////////
   
    pc->multipleBuildingsAdvect(odd,windField,positions0,positions1,prime0,prime1,
    			 randomValues,lambda,tau_dz,duvw_dz,time_step,buildings, cellType);
   
    ////////////////////////////////////////////////////////////
    // Get Position for Streams
    ////////////////////////////////////////////////////////////
    //if(stream->doUpdate()){
    //stream->updateStreamPos();
    //}
    
    ////////////////////////////////////////////////////////////
    // Update Mean Velocities
    ////////////////////////////////////////////////////////////
    if(util->calculateMeanVel){
      if(maxColorAttachments <= 4){
	FramebufferObject::Disable();
	fbo2->Bind();
      }
      pc->findMeanVel(odd,prime0,prime1,meanVel0,meanVel1,positions0,positions1,windField);

      if(maxColorAttachments <= 4){
	FramebufferObject::Disable();
	fbo->Bind();
      }
    }
    
    ///////////////////////////////////////////////////////////
    // Update Current Velocities
    ///////////////////////////////////////////////////////////
    if(maxColorAttachments <= 4){
      FramebufferObject::Disable();
      fbo2->Bind();
    }
    pc->updateCurrVel(odd,prime0,prime1,windField,positions0,positions1);

    if(maxColorAttachments <= 4){
      FramebufferObject::Disable();
      fbo->Bind();
    }

    ///////////////////////////////////////////////////////////
    // In some circumstances, we may want to dump the contents of
    // the FBO to a file

    if (dump_contents)
    {
      //pc->printPositions(odd);
      pathFbo->Bind();
      pathLines->printPathLineTexture();

      fbo->Bind();
      dump_contents = false;
    }
    if(print_MeanVel)
    {
      if(maxColorAttachments <= 4){
	FramebufferObject::Disable();
	fbo2->Bind();
      }
      pc->printMeanVelocities(odd);
      print_MeanVel = false;
      
      if(maxColorAttachments <= 4){
	FramebufferObject::Disable();
	fbo->Bind();
      }
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
    paused = true;
  }

  CheckErrorsGL("END : after 1st pass");
  // We only need to do PASS 2 (copy to VBO) and PASS 3 (visualize) if
  // we actually want to render to the screen.  Rendering to the
  // screen will make the simulation run more slowly. This feature is
  // mainly included to allow some idea of how much faster the
  // simulation can run if left to run on the GPU.
  if (util->show_particle_visuals)
    {
      glGetIntegerv(GL_READ_BUFFER, &read_buffer);
      // //////////////////////////////////////////////////////////////
      // PASS 2 - copy the contents of the 2nd texture (the new positions)
      // into the vertex buffer
      // //////////////////////////////////////////////////////////////
      if(odd){
	glReadBuffer(GL_COLOR_ATTACHMENT1_EXT);
      }
      else 
	glReadBuffer(GL_COLOR_ATTACHMENT0_EXT);

      glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB, vertex_buffer);
      glReadPixels(0, 0, twidth, theight, GL_RGBA, GL_FLOAT, 0);

      if(maxColorAttachments <= 4){
	FramebufferObject::Disable();
	fbo2->Bind();
      }
      glReadBuffer(currVelBuffer);
      
      glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB, color_buffer);
      glReadPixels(0, 0, twidth, theight, GL_RGBA, GL_FLOAT, 0);

      //update path line vbos
      pathFbo->Bind();
      pathLines->updateVBOS();

      glBindBuffer(GL_PIXEL_PACK_BUFFER_ARB, 0);
      CheckErrorsGL("after glReadPixels");
      glReadBuffer(read_buffer);

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
      
      // dc->drawVisuals(vertex_buffer, windField, color_buffer, numInRow, twidth, theight);
      dc->drawVisuals(vertex_buffer, duvw_dz, color_buffer, numInRow, twidth, theight);
      //stream->draw();
      pathLines->draw();

      glDisable(texType);
      if(dc->tau_visual == draw_contours){
	contours->draw();
	contours->displayContourLayer(pc,tau,numInRow);
      }
      glEnable(texType);

      if(dc->tau_visual == draw_layers)
	dc->drawLayers(windField,tau,numInRow);

      dc->drawScale();

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

void MultipleBuildingsModel::initFBO(void){
  // Get the Framebuffer Object ready
  fbo = new FramebufferObject();
  fbo->Bind();
      
  glGetIntegerv(GL_MAX_COLOR_ATTACHMENTS_EXT, (GLint*)&maxColorAttachments);
  std::cout << "Max color attachments: " << maxColorAttachments << std::endl;
  //rb = new Renderbuffer();
  //rb->Set(GL_DEPTH_COMPONENT24, twidth, theight);
  //fbo->AttachRenderBuffer(GL_DEPTH_ATTACHMENT_EXT, rb->GetId() );

  //Attach textures to framebuffer object
  fbo->AttachTexture(GL_COLOR_ATTACHMENT0_EXT, texType, texid[0]);
  fbo->AttachTexture(GL_COLOR_ATTACHMENT1_EXT, texType, texid[1]);
  fbo->AttachTexture(GL_COLOR_ATTACHMENT2_EXT, texType, prime0);
  fbo->AttachTexture(GL_COLOR_ATTACHMENT3_EXT, texType, prime1);

  if(maxColorAttachments > 4){
    fbo->AttachTexture(GL_COLOR_ATTACHMENT4_EXT, texType, meanVel0);
    fbo->AttachTexture(GL_COLOR_ATTACHMENT5_EXT, texType, meanVel1);
    fbo->AttachTexture(GL_COLOR_ATTACHMENT6_EXT, texType, currVel);

    currVelBuffer = GL_COLOR_ATTACHMENT6_EXT;
    pc->currVelBuffer = GL_COLOR_ATTACHMENT6_EXT;
    pc->meanVelBuffer0 = GL_COLOR_ATTACHMENT4_EXT;
    pc->meanVelBuffer1 = GL_COLOR_ATTACHMENT5_EXT;

  }

  fbo->IsValid();
  FramebufferObject::Disable();

  //If the max number of textures that can be attached to the fbo 
  //is <= 4, then we have to use another fbo.
  if(maxColorAttachments <= 4){
    fbo2 = new FramebufferObject();
    fbo2->Bind();

    fbo2->AttachTexture(GL_COLOR_ATTACHMENT0_EXT, texType, meanVel0); 
    fbo2->AttachTexture(GL_COLOR_ATTACHMENT1_EXT, texType, meanVel1);
    fbo2->AttachTexture(GL_COLOR_ATTACHMENT2_EXT, texType, currVel);
    CheckErrorsGL("FBO init 2");

    currVelBuffer = GL_COLOR_ATTACHMENT2_EXT;
    pc->currVelBuffer = GL_COLOR_ATTACHMENT2_EXT;
    pc->meanVelBuffer0 = GL_COLOR_ATTACHMENT0_EXT;
    pc->meanVelBuffer1 = GL_COLOR_ATTACHMENT1_EXT;

    fbo2->IsValid();
    FramebufferObject::Disable();
  }

  
  pathFbo = new FramebufferObject();
  pathFbo->Bind();
  pathFbo->AttachTexture(GL_COLOR_ATTACHMENT0_EXT, texType, pathTex);
  pathFbo->IsValid();
  FramebufferObject::Disable();

}

void MultipleBuildingsModel::setupTextures(){
  std::ofstream out;
  out.open(texFileName.c_str());

  CheckErrorsGL("BEGIN : Creating textures");

  int sz = 4;
  GLfloat *data = new GLfloat[ twidth * theight * sz];
		
  // Creates wind field data texture
  //The variable numInRow gets set here and should not be changed after being set!
  pc->initWindTex(windField, &numInRow, util->windFieldData);
  CheckErrorsGL("\tcreated texid[3], the wind field texture...");

  //Creates lambda, tau/dz, and duvw/dz textures
  if(util->windFieldData < 5)
    pc->initLambda_and_TauTex(lambda, tau_dz, duvw_dz);
  else if(util->windFieldData == 5)
    pc->initLambda_and_TauTex_fromQUICFILES(windField, lambda, tau_dz, duvw_dz, tau);
  else
    pc->initLambda_and_Taus_withCalculations(windField, lambda, tau_dz, duvw_dz, tau);
  CheckErrorsGL("\tcreated texid[7], the lambda texture...");

  //Print out max and min taus
  /*std::cout << "Max and Min Taus" << std::endl;
  std::cout << pc->tauMax[0] << " " << pc->tauMax[1] << " " << pc->tauMax[2] << " " << 
    pc->tauMax[3] << std::endl;
  
  std::cout << pc->tauMin[0] << " " << pc->tauMin[1] << " " << pc->tauMin[2] << " " << 
  pc->tauMin[3] << std::endl;*/

  pc->addBuildingsInWindField(cellType);
  CheckErrorsGL("\tcreated texid[14], the cell type texture...");

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
	
  int xs = (int)util->xpos[0];
  int ys = (int)util->ypos[0];
  int zs = (int)util->zpos[0];

  int p2idx = zs*ny*nx + ys*nx + xs;

  util->sigU = pc->sig[p2idx].u;
  util->sigV = pc->sig[p2idx].v;
  util->sigW = pc->sig[p2idx].w;

  int iterations = 0;
  double mean = 1.0;
  double stddev = 0.0, variance = 0.0, tmp_sum = 0.0;

  while( !(-0.01 < mean && mean < 0.01 && 0.90 < variance && variance < 1.01) && iterations < 50){
    iterations++;

    for (int j=0; j<theight; j++)
      for (int i=0; i<twidth; i++)
	{
	  int idx = j*twidth*sz + i*sz;

	  data[idx] = Random::normal();
	  data[idx+1] = Random::normal();
	  data[idx+2] = Random::normal();
	  data[idx+3] = 0.0;
	  
	  out << data[idx] << "\n";
	  out << data[idx+1] << "\n";
	  out << data[idx+2] << "\n";

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
  pc->createTexture(prime1, int_format, twidth, theight, data);
  CheckErrorsGL("\tcreated texid[5], the initial prime value texture...");
  pc->createTexture(meanVel0, int_format, twidth, theight, data);
  pc->createTexture(meanVel1, int_format, twidth, theight, data);
  CheckErrorsGL("\tcreated mean velocity textures...");

  for (int j=0; j<theight; j++)
    for (int i=0; i<twidth; i++)
      {
	int idx = j*twidth*sz + i*sz;

	if(i%2 == 0){
	  data[idx] = 1.0;
	  data[idx+1] = 0.0;
	  data[idx+2] = 0.0;
	  data[idx+3] = 1.0;
	}
	else{
	  data[idx] = 1.0;
	  data[idx+1] = 1.0;
	  data[idx+2] = 0.0;
	  data[idx+3] = 1.0;
	}
      }

  pc->createTexture(currVel, int_format, twidth, theight, data);

  //
  // create random texture for use with particle simulation and turbulence
  //
  random_values.clear();

  iterations = 0;
  mean = 1.0;
  variance = 0.0;

  while( !(-0.01 < mean && mean < 0.01 && 0.90 < variance && variance < 1.01) && iterations < 50){

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
		

	  out << data[idx] << "\n";
	  out << data[idx+1] << "\n";
	  out << data[idx+2] << "\n";

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
  
  //Building Texture
  
  GLfloat *bdata = new GLfloat[util->numBuild*2*sz];
  
  for(int j=0; j < util->numBuild; j++){
    int idx = j*2*sz;
          
    
    bdata[idx] = util->xfo[j];
    bdata[idx+1] = util->yfo[j];
    bdata[idx+2] = util->zfo[j];
    bdata[idx+3] = 0.0;
    
    bdata[idx+4] = util->ht[j];
    bdata[idx+5] = util->wti[j];
    bdata[idx+6] = util->lti[j];
    bdata[idx+7] = 0.0;

    /*//Second plane
    bdata[idx+4] = util->xfo[j]+util->lti[j];
    bdata[idx+5] = 0.0;
    bdata[idx+6] = 0.0;
    bdata[idx+7] = 0.0;
    //Third plane
    bdata[idx+8] = util->xfo[j];
    bdata[idx+9] = util->yfo[j]+(util->wti[j]/2.0);
    bdata[idx+10] = 0.0;
    bdata[idx+11] = 0.0;
    //Fourth plane
    bdata[idx+12] = util->xfo[j];
    bdata[idx+13] = util->yfo[j]-(util->wti[j]/2.0);
    bdata[idx+14] = 0.0;
    bdata[idx+15] = 0.0;
    //Fifth plane
    bdata[idx+16] = util->xfo[j];
    bdata[idx+17] = 0.0;
    bdata[idx+18] = util->zfo[j]+util->ht[j];
    bdata[idx+19] = 0.0;*/
    
  }
  pc->createTexture(texid[13], int_format, 2, util->numBuild, bdata);

  
  CheckErrorsGL("\tcreated texid[13], the building param texture...");

  delete [] bdata;

  GLfloat *pdata = new GLfloat[ util->pwidth * util->pheight * 4];
  for (int j=0; j<util->pheight; j++)
    for (int i=0; i<util->pwidth; i++)
      {
	int idx = j*util->pwidth*4 + i*4;
	pdata[idx] = 0.0;
	pdata[idx+1] = 0.0;
	pdata[idx+2] = 0.0;
	pdata[idx+3] = 1.0;
      }
  
  glBindTexture(texType, pathTex);
  
  pc->createTexture(pathTex,int_format,util->pwidth,util->pheight,pdata);
  
  delete [] pdata;


  CheckErrorsGL("END : Creating textures");
 
}

