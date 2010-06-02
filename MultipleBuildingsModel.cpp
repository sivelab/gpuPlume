#include <math.h>
#include "MultipleBuildingsModel.h"
#include "glErrorUtil.h"
#include <fstream>
#include<cstring>
#include<sstream>
#include<stdlib.h>
#include<stdio.h>

#include "PNGImage.h"

#ifdef WIN32
#include <fstream>
#include <windows.h>
#include <stdio.h>
#include <conio.h>
#include <tchar.h>

#endif
GLfloat* int_buffer;
FILE * fp; 
std::ofstream outputFile;
std::vector<float> str;

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

  nxdx = (int)(nx*(1.0/util->dx));
  nydy = (int)(ny*(1.0/util->dy));
  nzdz = (int)(nz*(1.0/util->dz));
  
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
  drawIsoSurface = false;
  color_by_advect_terms = false;
  //stream = new StreamLine(twidth,theight,nx,ny,nz);
  
  //Set whether to reuse particles or not
  //If reuseParticles is set to false: fourth coordinate of particle is -1 if emitted, else 0
  //If reuseParticles is set to true: fourth coordinate is <= lifetime if emiited, else lifetime+1
  reuseParticles = false;
  //Set this flag for continuous flow of particles.
  continuousParticleFlow = util->reuse_particles;
  std::cout << "Particle Reuse is " << ((continuousParticleFlow == true) ? "ON" : "OFF") << std::endl;
  frameCount = 0;
  if(reuseParticles)
    lifeTime = 1.0;
  else lifeTime = -1.0;

  for(int i = twidth*theight-1; i >= 0; i--)
    indices.push_back(i);

  oneTime = 0;

  mbaTimer = new Timer(true);

  std::string str1=util->quicFilesPath+"QU_deposition.dat";
  // std::cout << "Attempting to open deposition file for writing: \"" << str1 << "\"." << std::endl;

  //outputFile.open(str1.c_str(),std::ios::out);
  fp=fopen(str1.c_str(),"wb");

  // Set the default sun angle's
  sun_azimuth = util->sun_azimuth;
  sun_altitude = util->sun_altitude;
  
}

MultipleBuildingsModel::~MultipleBuildingsModel()
{
  fclose(fp);
  //outputFile.close();
  
  // Preform cleanup for the shadow map.
  shadowFBO->Unattach(GL_DEPTH_ATTACHMENT_EXT);
  glDeleteTextures(1, &(dc->shadowMap));
  delete shadowFBO;
  delete cellInShadowShader;
  delete dc->sunProjectionMatrix;
  delete dc->sunModelviewMatrix;
}

void MultipleBuildingsModel::init(bool OSG){
  osgPlume = OSG;
  
  pathLines = new PathLine(util->pwidth,util->pheight,texType);

  pc = new ParticleControl(texType, twidth,theight,nx,ny,nz,util->dx,util->dy,util->dz, util);
  pc->setUstarAndSigmas(util->ustar);
  pc->setBuildingParameters(util->numBuild,util->numSides,util->xfo,util->yfo,util->zfo,util->ht,util->wti,util->lti,util->gamma);
  pc->setQuicFilesPath(util->quicFilesPath);

  inPauseMode = util->pauseMode;
  dc = new DisplayControl(nx, ny, nz, texType, inPauseMode, util);
  dc->initVars(util->numBuild,util->numSides,util->xfo,util->yfo,util->zfo,util->ht,util->wti,util->lti,util->gamma);
  
  //??Don't remember why I wanted to do this
  //Send copy of particle emitter to displayControl
  //dc->setEmitter(pe[0]);

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
  //Texture used to hold current particle direction
  currDirection = texid[12];
  //Texture used for building information
  buildings = texid[13];
  //Texture used for cell type information
  cellType = texid[14];
 
  //Texture to store path lines
  pathTex = texid[15];
  
  //Texture for tau11,tau22,tau33,and tau13
  tau = texid[16];

  advect_terms = texid[17];
  /////////////////////////////
  setupTextures(); 
  /////////////////////////////

  

  setupEmitters();
    
  //Create isocontours
  contours = new Contour(pc,util->num_contour_regions);

  //Create visual plane class to visualize tau layers in the 3D domain
  planeVisual = new VisualPlane(pc, pc->tauMax, pc->tauMin,
				pc->tauLocalMax, pc->tauLocalMin);

 
  glDisable(texType);
  

  //glTexImage2D(GL_PROXY_TEXTURE_RECTANGLE_ARB,0,int_format,8192,8192,0,GL_RGBA,GL_FLOAT,NULL);
  //int wid = 0;
  //glGetTexLevelParameteriv(GL_PROXY_TEXTURE_RECTANGLE_ARB,0,GL_TEXTURE_WIDTH,&wid);
  //glGetIntegerv(GL_MAX_RECTANGLE_TEXTURE_SIZE_ARB, &wid);
  
  //std::cout << "max texture size = " << wid << std::endl;
  
  //Create isosurface

  isoSurface = new IsoSurface(pc);


  //
  // set up vertex buffer
  // 
  glGenBuffersARB(2, vbo_buffer);
  vertex_buffer = vbo_buffer[0];
  color_buffer = vbo_buffer[1];
  
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, vertex_buffer);
  glBufferDataARB(GL_ARRAY_BUFFER_ARB, twidth*theight*4*sizeof(GLfloat), 0, GL_STREAM_COPY);

  if (util->updateParticleColors)
    {
      glBindBufferARB(GL_ARRAY_BUFFER_ARB, color_buffer);
      glBufferDataARB(GL_ARRAY_BUFFER_ARB, twidth*theight*4*sizeof(GLfloat), 0, GL_STREAM_COPY);
    }
  
  glBindBuffer(GL_ARRAY_BUFFER, 0);
 

  //Initialize FBO
  initFBO();

  if(util->windFieldData >= 5)
    pc->setupMultipleBuildingsShader(lifeTime,0);
  else
    pc->setupMultipleBuildingsShader(lifeTime,1);

  pc->setupMeanVel_shader();

  pc->setupParticleColor_shader();

  //This shader is used to emmit particles
  emit_shader.addShader("Shaders/emitParticle_vp.glsl", GLSLObject::VERTEX_SHADER);
  emit_shader.addShader("Shaders/emitParticle_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  emit_shader.createProgram();

  CheckErrorsGL("END of init");

  // display_clock = new Timer(true);
  display_clock = new Timer();
  if(util->useRealTime)
    {
      sim->init();
    } 

#if 0
  int maxtextures;
  glGetIntegerv(GL_MAX_TEXTURE_IMAGE_UNITS, (GLint*)&maxtextures);
  std::cout << "Max Textures: " << maxtextures << std::endl;
#endif
  
  // Setup the shadowMap including the shadow map texutre and FBO.
  shadowMapSetup();

  // Set up the wind field lookup shader.
  pc->setupWindFieldLookupShader();

  // Initialize the treadport.
  dc->initTreadport();
}

int imgCounter = 0;
bool img_notDone = true;

int MultipleBuildingsModel::display(){
  
  // Timer_t displayStart = mbaTimer->tic();    

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
  
  if(isoSurface->once){
    //render the 3D density function texture
    isoSurface->render3DTexture(isoFbo);
    isoSurface->once = false;
  }


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
    
    // Store simulation start time and turn on one particle emitter
    if(firstTime){
      if(!osgPlume) {
	      for (int pe_i = 0; pe_i < util->numOfPE; pe_i++)
	      pe[pe_i]->emit = true;
	    }
      sim->setStartTime(&time_step);
      firstTime = false;
    }

    /////////////////////////////////////////////////////////////
    //Reuse particles if their lifespan is up
    /////////////////////////////////////////////////////////////
    if(reuseParticles)
      {
        particleReuse();
        CheckErrorsGL("MBA: display() - reuse particles");
      } 
    
    ////////////////////////////////////////////////////////////
    // Emit Particles
    ////////////////////////////////////////////////////////////
    for(int i = 0; i < util->numOfPE; i++){
      if(pe[i]->emit){    
	totalNumPar += (double)pe[i]->EmitParticle(odd,positions0,positions1,time_step,
						   prime0,prime1);

	CheckErrorsGL("MBA: display() - emit particle");

	if(pe[i]->releaseType == onePerKeyPress){
	  //stream->addNewStream(pe[i]);
	  pathLines->addNewPath(pe[i]);

	}
      }
    }
   
    ////////////////////////////////////////////////////////////
    // Update Path Lines
    ////////////////////////////////////////////////////////////
    FramebufferObject::Disable();
    pathFbo->Bind();
    
    pathLines->updatePathLines(positions0,positions1,odd);
    CheckErrorsGL("END : pathLines");
    //Make sure to bind this fbo, because the pathLines use a
    //different fbo
    FramebufferObject::Disable();
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
   
    /*else{*/
    pc->multipleBuildingsAdvect(odd,windField,positions0,positions1,prime0,prime1,
				randomValues,lambda,tau_dz,duvw_dz,time_step,buildings, 
				cellType, color_by_advect_terms);
    //}

      CheckErrorsGL("MBA: display() - advection");
     
      if(color_by_advect_terms){
	int_buffer = new GLfloat[ twidth * theight * 4 ];
	glReadBuffer(GL_COLOR_ATTACHMENT7_EXT);
	glReadPixels(0, 0, twidth, theight, GL_RGBA, GL_FLOAT, int_buffer); 
	for(int i = 3; i <= (theight*twidth*4); i+=4){
	  //If particle has been reflected
	  if(int_buffer[i] == 2.0){
	
	    //Get the x,y,z position of the particle
	    float x = int_buffer[i-3];
	    float y = int_buffer[i-2];
	    float z = int_buffer[i-1];
	    //fwrite(x,sizeof(x),1,fp);
	    fprintf(fp,"%f %f %f\n",x,y,z);
	    //outputFile<<float(x)<<"\t"<<float(y)<<"\t"<<float(z)<<"\n";

	    /*char* ch=new char[10];
	      char* ch1=new char[10];
	      char* ch2=new char[10];
	      gcvt(x,4,ch);
	      gcvt(y,4,ch1);
	      gcvt(z,4,ch2);
	      outputFile.write(ch,sizeof(ch));
	      outputFile.write(ch,sizeof(ch)); 
	      outputFile.write(ch,sizeof(ch));*/
	  }
	}    
	delete [] int_buffer;
      }
 
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

      CheckErrorsGL("MBA: display() - calc mean vel");

      if(maxColorAttachments <= 4){
	FramebufferObject::Disable();
	fbo->Bind();
      }
    }
    
    ///////////////////////////////////////////////////////////
    // Update Particle Colors
    ///////////////////////////////////////////////////////////
    if (util->updateParticleColors)
      {
	if(maxColorAttachments <= 4){
	  FramebufferObject::Disable();
	  fbo2->Bind();
	}
	pc->updateParticleColors(odd,prime0,prime1,windField,positions0,positions1);

	CheckErrorsGL("MBA: display() - update particle colors");

	if(maxColorAttachments <= 4){
	  FramebufferObject::Disable();
	  fbo->Bind();
	}
      }

    ///////////////////////////////////////////////////////////
    // In some circumstances, we may want to dump the contents of
    // the FBO to a file

    if (dump_contents)
    {
      pc->printPositions(odd);
      //pathFbo->Bind();
      //pathLines->printPathLineTexture();

      //fbo->Bind();
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
				// outputs concentration in grams per cubic meter
				cBoxes[j]->outputConcStd(util->output_file, util->output_id, 
				util->averagingTime,
				util->volume,
				util->time_step, 
				(util->twidth)*(util->theight));// standard Concentration Calc. - Balli
      }
      output_CollectionBox = false;

      // Pete
      // need to add not as function but tranformed data to the end of the cbox.m file...

      // output matlab code to remove the collections boxes in the
      // emitter and buildings

      std::ostringstream cStringId;
      cStringId << util->output_id << "_conc";

      std::ofstream mfunc_output;
      mfunc_output.open(util->output_file.c_str(), std::ios::app);

  mfunc_output << "x=" << cStringId.str() << "(:,1);" << std::endl;
  mfunc_output << "y=" << cStringId.str() << "(:,2);" << std::endl;
  mfunc_output << "z=" << cStringId.str() << "(:,3);" << std::endl;
  mfunc_output << "con=" << cStringId.str() << "(:,4);" << std::endl;

  mfunc_output << "xUni=unique(x);" << std::endl;
  mfunc_output << "yUni=unique(y);" << std::endl;
  mfunc_output << "zUni=unique(z);" << std::endl;
  mfunc_output << "xLen = length(xUni);" << std::endl;
  mfunc_output << "yLen = length(yUni);" << std::endl;
  mfunc_output << "zLen = length(zUni);" << std::endl;
  mfunc_output << "[X Y] = meshgrid(xUni, yUni);" << std::endl;


  // only care about first three layers for now
  mfunc_output << "for i=1:3" << std::endl;
  mfunc_output << "fid = figure;" << std::endl;
  mfunc_output << "concSlice=con(z==zUni(i));" << std::endl;
  mfunc_output << "C=reshape(concSlice,xLen,yLen);" << std::endl;
  mfunc_output << "pcolor(C')" << std::endl;
  mfunc_output << "colorbar;" << std::endl;
  mfunc_output << "title(['Concentration - Z = ',num2str(zUni(i)),'m'])" << std::endl;
  mfunc_output << "set(gcf,'color','w')" << std::endl;

  mfunc_output << "filename = sprintf('" << cStringId.str() << "_%02d.png', i);" << std::endl;
  mfunc_output << "print('-dpng', filename);" << std::endl;

  mfunc_output << "end" << std::endl;

      mfunc_output.close();
    }

    // Switches the frame buffer and binding texture
    odd = !odd;
    paused = true;
  }

  CheckErrorsGL("MBA: display() - completed all advection, prior to display");

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

      /*if(color_by_advect_terms){
	glReadBuffer(GL_COLOR_ATTACHMENT7_EXT);
	glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB, color_buffer);
	glReadPixels(0, 0, twidth, theight, GL_RGBA, GL_FLOAT, 0);
	
      }
      else{*/
	if (util->updateParticleColors)
	  {
	    if(maxColorAttachments <= 4){
	      FramebufferObject::Disable();
	      fbo2->Bind();
	    }
	    glReadBuffer(particleColorBuffer);
      
	    glBindBufferARB(GL_PIXEL_PACK_BUFFER_ARB, color_buffer);
	    glReadPixels(0, 0, twidth, theight, GL_RGBA, GL_FLOAT, 0);
	  }
      //}
      
      //update path line vbos
      pathFbo->Bind();
      pathLines->updateVBOS();

      glBindBuffer(GL_PIXEL_PACK_BUFFER_ARB, 0);
      CheckErrorsGL("after glReadPixels");
      glReadBuffer(read_buffer);

      // Disable the framebuffer object
      FramebufferObject::Disable();
      if (util->updateParticleColors == false)
	glColor3f(1.0, 0.0, 0.0);  // general color if we're not supplying it from the FBO textures.

      glDrawBuffer(draw_buffer); // send it to the original buffer

      /////////////////////////////////////////////////////
      //Render geometry shader outputs to the vertex buffer
      /////////////////////////////////////////////////////
      CheckErrorsGL("before isosurface");

      if(oneTime < 1){
	isoSurface->createIsoSurface();
	oneTime++;
      }

      CheckErrorsGL("END : after 2nd pass");

      // //////////////////////////////////////////////////////////////
      // PASS 3 - draw the vertices; This represents the visualization
      // of the PLUME particle field.
      // //////////////////////////////////////////////////////////////
      
      // Grab the shadow map. Note, this should really only be done
      // every time the light source moves and does not need to be 
      // done every frame (large waste).
      if(reCalcShadows) {
				generateShadowMap();
				
				for(int i = 0; i < util->nz; i++) {
	  			genGridShadow(i, 0);
	  			genGridShadow(i, 1);
				}
				
				reCalcShadows = false;
      }

      if(util->onlyCalcShadows) {
				writeShadowMapToFile();
				exit(0);
      }

      // GLfloat windDir[3];
      GLfloat pos[3];
      
      pos[0] = dc->eye_pos[0];
      pos[1] = dc->eye_pos[1];
      pos[2] = dc->eye_pos[2];
      
      // pc->lookupWindField(pos, windField, dc->windDir[0], dc->windDir[1], dc->windDir[2]);

      // std::cout << windDir[0] << " " << windDir[1] << " " << windDir[2] << std::endl;
      // int x;
      // std::cin >> x;

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
      
      // Synchronized the data with the treadport system, if we are not
      // running in treadport mode, then this function does nothing.
      dc->syncTreadportData();
      
      // Synchronize the data over the network, if we don't have 
      // networking enabled, this does nothing.
      dc->syncDataOverNetwork();
      
      // Now that the values have been swapped over the network, we can
      // set them here.
      inPauseMode = dc->getInPauseMode();

      // Initialize the view frustum.
      dc->initializeView();
			
      // Draw the visuals (done in DisplayControl).
      dc->drawVisuals(vertex_buffer, duvw_dz, color_buffer, numInRow, twidth, theight, texid[0], prime0);
			
      // Uninitialize view (this is done because of treadport compatibility).
      dc->deinitializeView();

      CheckErrorsGL("MBA : called drawVisuals");

#if 0
      //stream->draw();
      if (pathLines)
	{
	  pathLines->draw();
	  CheckErrorsGL("MBA : after drawing pathlines");
	}
#endif

      glDisable(texType);
      if(dc->tau_visual == draw_contours){
				//contours->draw();
				contours->displayContourLayer(pc,tau,numInRow);
      }
      glEnable(texType);

      if(dc->tau_visual == draw_layers){
				if(planeVisual->visual_field > 0)
	  			if(planeVisual->rotationPlane)
	    			planeVisual->drawRotationalPlane();
	  			else
	    			planeVisual->drawAxisAlignedPlane();
				else
	  				dc->drawLayers(windField,numInRow, util->calculatedMaxVel);
				CheckErrorsGL("MBA : after drawing layers");
      }
      
      //Draw isosurface
      if(drawIsoSurface){
			isoSurface->draw();
      }

#if 0
      // Draw the wind field vectors...
      glBindBufferARB(GL_ARRAY_BUFFER_ARB, pc->windFieldVector_vbo);
      glColor3f(1.0, 0.0, 0.0);
      glVertexPointer(3, GL_FLOAT, 0, 0);
      glEnableClientState(GL_VERTEX_ARRAY);
      // glDrawArrays(GL_LINES, 0, pc->windFieldVector_w*pc->windFieldVector_h);
      glDrawArrays(GL_POINTS, 0, 10);
      glDisableClientState(GL_VERTEX_ARRAY);
      glBindBuffer(GL_ARRAY_BUFFER, 0);
#endif

      //      planeVisual->drawScale();

      if(!osgPlume){
	for(int i=0; i < util->numOfPE; i++){
	  pe[i]->Draw();      
	  CheckErrorsGL("MBA : after drawing particle emitters");
	}

	if(util->show_collectionBox_visuals){
	  cBoxes[0]->sort(dc->eye_pos[0],dc->eye_pos[1],dc->eye_pos[2]);
	  cBoxes[0]->draw(sim->curr_timeStep);
	}
      }


      // ////////////////////////////////////////////
#if 0
      if (img_notDone && imgCounter > 2000)
	{
	  // int iWidth = glutGet(GLUT_WINDOW_WIDTH);
	  // int iHeight = glutGet(GLUT_WINDOW_HEIGHT);

	  int iWidth = 400;
	  int iHeight = 200;
	  GLfloat *imgBuffer = new GLfloat[ iWidth * iHeight * 3 ];
	  glReadBuffer(GL_BACK);
	  glReadPixels(0, 0, iWidth, iHeight, GL_RGB, GL_FLOAT, imgBuffer); 
	  
	  cs5721::PNGImage pngimage;
	  pngimage.writeFileData("plume.png", iWidth, iHeight, imgBuffer);
	  std::cout << "Wrote PNG file" << std::endl;

	  delete [] imgBuffer;
	  img_notDone = false;
	}
      else
	imgCounter++;
#endif


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
      if(util->offscreenRender == false && (!osgPlume))
	glutSwapBuffers();
    }
  else{
    FramebufferObject::Disable();
    //glDisable(texType);
    glDrawBuffer(draw_buffer); // send it to the original buffer
    CheckErrorsGL("END : after not showing visuals");

    if(util->offscreenRender == false)
      glutSwapBuffers();
  }

  // Timer_t displayEnd = mbaTimer->tic();    
  //  std::cout << "MBA Display Time: " << mbaTimer->deltau(displayStart, displayEnd) << " us." << std::endl;  

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
  // std::cout << "Max color attachments: " << maxColorAttachments << std::endl;
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
    fbo->AttachTexture(GL_COLOR_ATTACHMENT6_EXT, texType, currDirection);

    particleColorBuffer = GL_COLOR_ATTACHMENT6_EXT;
    pc->particleColorBuffer = GL_COLOR_ATTACHMENT6_EXT;
    pc->meanVelBuffer0 = GL_COLOR_ATTACHMENT4_EXT;
    pc->meanVelBuffer1 = GL_COLOR_ATTACHMENT5_EXT;

    fbo->AttachTexture(GL_COLOR_ATTACHMENT7_EXT, texType, advect_terms);

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
    fbo2->AttachTexture(GL_COLOR_ATTACHMENT2_EXT, texType, currDirection);
    CheckErrorsGL("FBO init 2");

    particleColorBuffer = GL_COLOR_ATTACHMENT2_EXT;
    pc->particleColorBuffer = GL_COLOR_ATTACHMENT2_EXT;
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

  isoFbo = new FramebufferObject();
  isoFbo->Bind();
  isoFbo->AttachTexture(GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_3D, isoSurface->tex3d[0]);
  isoFbo->IsValid();
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

  //These two textures are to store the prime values(previous and updated values)
  //We will need to initialize some data into prime0

  //idx lookup using the source position to determine
  //the sigma values at that point in space. Only accurate
  //if there's only one source and it's a point emitter.
	
  int xs = (int)util->xpos[0];
  int ys = (int)util->ypos[0];
  int zs = (int)util->zpos[0];

  int p2idx = zs*nydy*nxdx + ys*nxdx + xs;

  util->sigU = pc->sig[p2idx].u;
  util->sigV = pc->sig[p2idx].v;
  util->sigW = pc->sig[p2idx].w;

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
	  //data[idx] = 100.0;   
	  //data[idx+1] = 100.0;
	  //data[idx+2] = 100.0;
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
  // std::cout << "Number of iterations to get random values: " << iterations << std::endl;

  pc->createTexture(prime0, int_format, twidth, theight, data);
  pc->createTexture(prime1, int_format, twidth, theight, data);

  //Create advect_terms texture
  pc->createTexture(advect_terms, int_format, twidth, theight, data);

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

  pc->createTexture(currDirection, int_format, twidth, theight, data);

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

  std::cout << "Creating Random Values Texture: Mean = " << mean << ", Standard Deviation = " << stddev << ", Variance = " << variance << std::endl;
  // std::cout << "Number of iterations to get random values: " << iterations << std::endl;

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
    bdata[idx+3] = util->numSides[j];
    
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

void MultipleBuildingsModel::shadowMapSetup()
{
  
  reCalcShadows = true;
  
  // Initialize the sun matricies
  dc->sunModelviewMatrix = new GLfloat[16];
  dc->sunProjectionMatrix = new GLfloat[16];

  //
  // Initialize the shadow map texture.
  //
  glGenTextures(1, &(dc->shadowMap));
  glBindTexture(GL_TEXTURE_2D, dc->shadowMap);
  
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  
  // This is to allow the use of the shadow2DProj function in the shader.
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_R_TO_TEXTURE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
  glTexParameteri(GL_TEXTURE_2D, GL_DEPTH_TEXTURE_MODE, GL_INTENSITY);

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  
  glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, 2048, 2048, 0, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, NULL);
  
  //
  // Initialize the framebuffer and assign the shadow map texture to it.
  //
  shadowFBO = new FramebufferObject();
  shadowFBO->AttachTexture(GL_DEPTH_ATTACHMENT_EXT, GL_TEXTURE_2D, dc->shadowMap);

  glDrawBuffer(GL_NONE);
  glReadBuffer(GL_NONE);

  shadowFBO->IsValid();

  FramebufferObject::Disable();

  // Create the shader and load the programs.
  cellInShadowShader = new GLSLObject();
  cellInShadowShader->addShader("Shaders/cellInShadow_vp.glsl", GLSLObject::VERTEX_SHADER);
  cellInShadowShader->addShader("Shaders/cellInShadow_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  cellInShadowShader->createProgram();
  
}

void MultipleBuildingsModel::rotatePoint(float (& pos)[3], float axis[3], float angle) {
  
  float c = cos(angle);
  float s = sin(angle);
  float rotMat[4][4];
  
  // Setup the rotation matrix, this matrix is based off of the rotation matrix used in glRotatef.
  rotMat[0][0] = axis[0] * axis[0] + (1 - axis[0] * axis[0]) * c;           rotMat[0][1] = axis[0] * axis[1] * (1 - c) - axis[2] * s; rotMat[0][2] = axis[0] * axis[2] * (1 - c) + axis[1] * s; rotMat[0][3] = 0;
  rotMat[1][0] = axis[1] * axis[0] * (1 - c) + axis[2] * s; rotMat[1][1] = axis[1] * axis[1] + (1 - axis[1] * axis[1]) * c;           rotMat[1][2] = axis[1] * axis[2] * (1 - c) - axis[0] * s; rotMat[1][3] = 0;
  rotMat[2][0] = axis[0] * axis[2] * (1 - c) - axis[1] * s; rotMat[2][1] = axis[1] * axis[2] * (1 - c) + axis[0] * s; rotMat[2][2] = axis[2] * axis[2] + (1 - axis[2] * axis[2]) * c;           rotMat[2][3] = 0;
  rotMat[3][0] = 0;                                         rotMat[3][1] = 0;                                         rotMat[3][2] = 0;                                         rotMat[3][3] = 1;
  
  float x, y, z;

  // Multiply the rotation matrix with the position vector.
  x = rotMat[0][0] * pos[0] + rotMat[0][1] * pos[1] + rotMat[0][2] * pos[2] + rotMat[0][3];
  y = rotMat[1][0] * pos[0] + rotMat[1][1] * pos[1] + rotMat[1][2] * pos[2] + rotMat[1][3];
  z = rotMat[2][0] * pos[0] + rotMat[2][1] * pos[1] + rotMat[2][2] * pos[2] + rotMat[2][3];
  
  pos[0] = x;
  pos[1] = y;
  pos[2] = z;
  
}

void MultipleBuildingsModel::generateShadowMap()
{
  
  // Set the sun position, so it can be drawn in the world.
  dc->sun_pos[0] = 0;
  dc->sun_pos[1] = 0;
  dc->sun_pos[2] = sqrt((util->nx * util->dx) * (util->nx * util->dx) 
                        + (util->ny * util->dy) * (util->ny * util->dy)
                        + (util->nz * util->dz) * (util->nz * util->dz)
                       );
  
  float axis[3];
  axis[0] = 0.0;
  axis[1] = 1.0;
  axis[2] = 0.0;
  rotatePoint(dc->sun_pos, axis, -(90 - sun_altitude) * M_PI / 180);
  
  axis[0] = 0.0;
  axis[1] = 0.0;
  axis[2] = 1.0;
  rotatePoint(dc->sun_pos, axis, -sun_azimuth * M_PI / 180);
  
  std::cout << "Altitude: " << sun_altitude << " Azimuth: " << sun_azimuth << std::endl;
  
  // Prep the scale and bias matrix.
  dc->sunScaleAndBiasMatrix[0] = 0.5;
  dc->sunScaleAndBiasMatrix[1] = 0.0;
  dc->sunScaleAndBiasMatrix[2] = 0.0;
  dc->sunScaleAndBiasMatrix[3] = 0.0;
  dc->sunScaleAndBiasMatrix[4] = 0.0;
  dc->sunScaleAndBiasMatrix[5] = 0.5;
  dc->sunScaleAndBiasMatrix[6] = 0.0;
  dc->sunScaleAndBiasMatrix[7] = 0.0;
  dc->sunScaleAndBiasMatrix[8] = 0.0;
  dc->sunScaleAndBiasMatrix[9] = 0.0;
  dc->sunScaleAndBiasMatrix[10] = 0.5;
  dc->sunScaleAndBiasMatrix[11] = 0.0;
  dc->sunScaleAndBiasMatrix[12] = 0.5;
  dc->sunScaleAndBiasMatrix[13] = 0.5;
  dc->sunScaleAndBiasMatrix[14] = 0.5;
  dc->sunScaleAndBiasMatrix[15] = 1.0;
  
  //
  // Draw the scene from the light's (Sun's) perspective.
  //
  
  shadowFBO->Bind();

  // Calculate the scene's bounding radius.
  float sceneBoundingRadius = 100.0;
  
  // Calculate the distance between the sceen and the light
  float lightToSceneDistance = sqrt(float(dc->sun_pos[0] * dc->sun_pos[0] + dc->sun_pos[1] * dc->sun_pos[1] + dc->sun_pos[2] * dc->sun_pos[2]));

  // Calculate the near plane for the light
  float lightNearPlane = lightToSceneDistance - sceneBoundingRadius;
  
  float lightFieldOfView = 180 * M_PI * ( 2.0f * atan(sceneBoundingRadius / lightToSceneDistance ));
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glViewport(0, 0, 2048, 2048);
  // gluPerspective(lightFieldOfView, 1.0f, lightNearPlane, lightNearPlane + (2.0f * sceneBoundingRadius));
  glOrtho(-1024, 1024, -1024, 1024, 1.0f, lightToSceneDistance * 2);
  glGetFloatv(GL_PROJECTION_MATRIX, dc->sunProjectionMatrix);
	
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(
	    dc->sun_pos[0], dc->sun_pos[1], dc->sun_pos[2],
	    0,   0,   0,
	    0,   0,   1
	   );
  glGetFloatv(GL_MODELVIEW_MATRIX, dc->sunModelviewMatrix);
  
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glClearDepth(1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  dc->drawFeatures();

  FramebufferObject::Disable();
  
  // Now that we are done capturing the image, we can unattach the texture
  // from the frame buffer.
  shadowFBO->Bind();
  shadowFBO->Unattach(GL_DEPTH_ATTACHMENT_EXT);
  FramebufferObject::Disable();

}

void MultipleBuildingsModel::genGridShadow(int i, int cellPoints) {
  
  // Create and initialize the array of cell positions (the position being the
  // center of each cell).
  GLfloat positions[util->nx][util->ny][4];
  for(int x = 0; x < util->nx; x++) {
    for(int y = 0; y < util->ny; y++) {
			positions[x][y][0] = x * util->dx + 0.5 * util->dx;
			positions[x][y][1] = y * util->dy + 0.5 * util->dy;
			positions[x][y][2] = (float)i*util->dz + 0.5*util->dz;
			positions[x][y][3] = 1.0f;
    }
  }

	// Create and initialize a texture to store the positions on the graphics
	// card.
  GLuint posTex;
  glGenTextures(1, &posTex);
  glBindTexture(GL_TEXTURE_RECTANGLE_ARB, posTex);
  glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, GL_RGBA32F_ARB, util->nx, util->ny, 0, GL_RGBA, GL_FLOAT, &positions);

  // Set up the output texture
  // NOTE NEED LOGIC TO USE THE LARGER OF NX NY!!! IF WE DON'T HAVE SQUARE DOMAIN.
  GLuint cellsInShadow;
  glGenTextures(1, &cellsInShadow);
  glBindTexture(GL_TEXTURE_2D, cellsInShadow);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, util->nx, util->ny, 0, GL_RGBA, GL_FLOAT, NULL);
  
  // Bind the output texture to the fbo.
  shadowFBO->AttachTexture(GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, cellsInShadow);
  shadowFBO->IsValid();
  FramebufferObject::Disable();
  
  // Setup that allows the shader to be able to read the data from our texture.
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(-1, 1, -1, 1);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glViewport(0, 0, util->nx, util->ny);
  
  // Pass the shadowMap texture to the shader along with
  // the other matricies and parameters.
  cellInShadowShader->activate();
  
  // This is just a temp variable which we can reuse.
  GLint tmpID;

  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D, dc->shadowMap);
  
  tmpID = cellInShadowShader->createUniform("shadowMap");
  glUniform1i(tmpID, 1);

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_RECTANGLE_ARB, posTex);

  tmpID = cellInShadowShader->createUniform("posTex");
  glUniform1i(tmpID, 0);

  tmpID = cellInShadowShader->createUniform("zPos");
  glUniform1f(tmpID, (float)i*util->dz + 0.5*util->dz);
  
  tmpID = cellInShadowShader->createUniform("dx");
  glUniform1f(tmpID, (float)util->dx);
  tmpID = cellInShadowShader->createUniform("dy");
  glUniform1f(tmpID, (float)util->dy);
  tmpID = cellInShadowShader->createUniform("dz");
  glUniform1f(tmpID, (float)util->dz);
  
  tmpID = cellInShadowShader->createUniform("cellPoints");
  glUniform1i(tmpID, cellPoints);
  
  tmpID = cellInShadowShader->createUniform("sunModelviewMatrix");
  glUniformMatrix4fv(tmpID, 1, GL_FALSE, dc->sunModelviewMatrix);
  
  tmpID = cellInShadowShader->createUniform("sunProjectionMatrix");
  glUniformMatrix4fv(tmpID, 1, GL_FALSE, dc->sunProjectionMatrix);
  
  tmpID = cellInShadowShader->createUniform("sunScaleAndBiasMatrix");
  glUniformMatrix4fv(tmpID, 1, GL_FALSE, (GLfloat*)&(dc->sunScaleAndBiasMatrix));
    
  glBegin(GL_QUADS);
  {
    glTexCoord2f(0, 0);			        glVertex3f(-1, -1, -0.5f);
    glTexCoord2f(util->nx, 0);			glVertex3f( 1, -1, -0.5f);
    glTexCoord2f(util->nx, util->ny);	        glVertex3f( 1,  1, -0.5f);
    glTexCoord2f(0, util->ny);			glVertex3f(-1,  1, -0.5f);
  }
  glEnd();
  
  // Disable shader here.
  cellInShadowShader->deactivate();

  // To use frame buffer specificy read buffer as color attachment 0.
  if(cellPoints == 0) {
    glReadPixels(0, 0, util->nx, util->ny, GL_RGBA, GL_FLOAT, &(dc->inShadowData[i*util->ny*util->nx*4]));
  } else {
    glReadPixels(0, 0, util->nx, util->ny, GL_RGBA, GL_FLOAT, &(dc->inShadowData2[i*util->ny*util->nx*4]));
  }
  // Cleanup
  glDeleteTextures(1, &cellsInShadow);
}


void MultipleBuildingsModel::swapPauseMode() {
  // We set this in DisplayControl and not here so
  // that the value is synchronized first before it
  // is actually set here. This way we won't be out
  // of sync with the other screens.
  dc->setInPauseMode(!inPauseMode);
}

void MultipleBuildingsModel::writeShadowMapToFile() {
  
  std::ostringstream fileNameBuffer;
  fileNameBuffer << "shadow-" << util->sun_altitude << "-" << util->sun_azimuth << ".txt";
  std::string fileName = fileNameBuffer.str();

  std::cout << "Writing shadow data to shadow.txt..." << std::flush;
  
  std::ofstream file;
  file.open(fileName.c_str(), std::ios_base::out);
  
  for(int i = 0; i < nz; i++) {
	 for(int j = 0; j < nx; j++) {
	  for(int k = 0; k < ny; k++) {
		 int index = i*nx*ny*4 + j*ny*4 + k*4;
		 if(dc->inShadowData[index] == 0) {
		  file << k * util->dy << " " << j * util->dx <<  " " << i * util->dz << " " << "B " << dc->inShadowData[index] << " ";
		 }
		 if(dc->inShadowData[index + 1] == 0) {
		  file << k * util->dy << " " << j * util->dx <<  " " << i * util->dz << " " << "N " << dc->inShadowData[index + 1] << " ";
		 }
		 if(dc->inShadowData[index + 2] == 0) {
		  file << k * util->dy << " " << j * util->dx <<  " " << i * util->dz << " " << "W " << dc->inShadowData[index + 2] << " ";
		 }
		 
		 if(dc->inShadowData2[index] == 0) {
		  file << k * util->dy << " " << j * util->dx <<  " " << i * util->dz << " " << "T " << dc->inShadowData2[index] << " ";
		 }
		 if(dc->inShadowData2[index + 1] == 0) {
		  file << k * util->dy << " " << j * util->dx <<  " " << i * util->dz << " " << "S " << dc->inShadowData2[index + 1] << " ";
		 }
		 if(dc->inShadowData2[index + 2] == 0) {
		  file << k * util->dy << " " << j * util->dx <<  " " << i * util->dz << " " << "E " << dc->inShadowData2[index + 2] << " ";
		 }
	  }
	 }
  }
  
  file.close();

  std::cout << "   done." << std::endl;
}
