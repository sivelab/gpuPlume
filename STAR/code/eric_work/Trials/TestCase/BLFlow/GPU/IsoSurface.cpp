#include <cstdlib>
#include "IsoSurface.h"
#include "glErrorUtil.h"

IsoSurface::IsoSurface(ParticleControl *pc){
  //Mesh size
  //more rigid to smooth
  // 1    -    32? (It probably won't work with 32 though.  4,8,16 work)
  mesh = 4;
  //Contour value 
  contourValue = 3.0;

  nx = pc->nx;
  ny = pc->ny;
  nz = pc->nz;
  
  nxdx = pc->nxdx;
  nydy = pc->nydy;
  nzdz = pc->nzdz;

  int_format = GL_RGBA32F_ARB;
  texType2 = GL_TEXTURE_3D;

  int num_vertices;
  if(GL_EXT_geometry_shader4)
    {
      std::cout << "Ready for geom shader!" << std::endl;
      //int numV = 0;
      //glGetIntegerv(GL_MAX_GEOMETRY_OUTPUT_VERTICES_EXT, &numV);
      
      glGetIntegerv(GL_MAX_GEOMETRY_OUTPUT_VERTICES_EXT, &num_vertices);
      std::cout << num_vertices << " number of output vertices." << std::endl;
    }
  geomShader.addShader("Shaders/isoSurface_vp.glsl", GLSLObject::VERTEX_SHADER);
  geomShader.addShader("Shaders/isoSurface_gp.glsl", GLSLObject::GEOMETRY_SHADER);
  geomShader.setInputandOutput(GL_POINTS,GL_TRIANGLE_STRIP,num_vertices);
  geomShader.createProgram();
  u_tau3D = geomShader.createUniform("tau");
  u_case = geomShader.createUniform("case_to_numpoly");
  u_edge = geomShader.createUniform("edge_connect_list");
  u_dx = geomShader.createUniform("dx");
  u_dy = geomShader.createUniform("dy");
  u_dz = geomShader.createUniform("dz");
  u_mesh = geomShader.createUniform("mesh");

  geomShader.activate();
  
  float scale = 1.0/(float)mesh;
  
  glUniform1fARB(u_dx,scale/(float)(nxdx));
  glUniform1fARB(u_dy,scale/(float)(nydy));
  glUniform1fARB(u_dz,scale/(float)(nzdz));
  glUniform1fARB(u_mesh,scale);

  geomShader.deactivate();

  int* attr = new int[1];
  attr[0] = GL_POSITION;
   
  geomShader.setVaryingOutput(1,attr,GL_INTERLEAVED_ATTRIBS_NV);
 
  render3DShader.addShader("Shaders/render3D_vp.glsl", GLSLObject::VERTEX_SHADER);
  render3DShader.addShader("Shaders/render3D_fp.glsl", GLSLObject::FRAGMENT_SHADER);
  render3DShader.createProgram();
  u_slice = render3DShader.createUniform("slice");
  uniform_tau = render3DShader.createUniform("tau");
  uniform_cValue = render3DShader.createUniform("contourValue");

  /////////////////////////////////////////////////
  //Create Textures
  /////////////////////////////////////////////////
  glEnable(texType2);

  glGenTextures(2,tex3d);

  //Texture to render to
  GLfloat* data = new GLfloat[(nxdx+1)*(nydy+1)*(nzdz+1)*4];
  for(int k=0;k<=nzdz;k++){
    for(int i=0; i <= nydy; i++){
      for(int j=0; j <= nxdx; j++){
	//int idx = k*ny*nx + i*nx + j;

	int tidx = k*nydy*nxdx*4 + i*nxdx*4 + j*4;
	
	data[tidx] = 0.0;
	data[tidx+1] = 0.0;
	data[tidx+2] = 0.0;
	data[tidx+3] = 1.0;
      }
    }
  }
  


  glBindTexture(texType2, tex3d[0]);
  //glTexParameteri(texType2, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  //glTexParameteri(texType2, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(texType2, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(texType2, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(texType2, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(texType2, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(texType2, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
  glTexImage3D(texType2, 0, int_format, nxdx+1, nydy+1, nzdz+1,0,GL_RGBA, GL_FLOAT,data); 
  
  delete [] data;

  GLfloat* data2 = new GLfloat[nxdx*nydy*nzdz*4];

  //Texture to read from
  for(int k=0;k<nzdz;k++){
    for(int i=0; i < nydy; i++){
      for(int j=0; j < nxdx; j++){
	int idx = k*nydy*nxdx + i*nxdx + j;

	int tidx = k*nydy*nxdx*4 + i*nxdx*4 + j*4;
	
	data2[tidx] = pc->tau[idx].t11;
	data2[tidx+1] = pc->tau[idx].t22;
	data2[tidx+2] = pc->tau[idx].t33;
	data2[tidx+3] = pc->tau[idx].t13;
      }
    }
  }
  glBindTexture(texType2, tex3d[1]);
  glTexParameteri(texType2, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(texType2, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  //glTexParameteri(texType2, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  //glTexParameteri(texType2, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(texType2, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(texType2, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(texType2, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
  glTexImage3D(texType2, 0, int_format, nxdx, nydy, nzdz,0,GL_RGBA, GL_FLOAT,data2); 

  delete[] data2;

  glBindTexture(texType2,0);
  glDisable(texType2);

   
  //Read in lookup tables and store in textures
  readInTables();
  
  glGenQueries(1,&query);
  numPrimitives = 0;

  //set up vertex buffer
  
  glGenBuffersARB(10, iso_buffer);
  buffer_num = 0;

  glBindBuffer(GL_ARRAY_BUFFER_ARB, iso_buffer[0]);
  //??The size of this buffer might need to be larger??
  glBufferData(GL_ARRAY_BUFFER_ARB, mesh*(nzdz)*(nxdx)*(nydy)*15*4*sizeof(GLfloat),0, GL_STREAM_COPY);
  glBindBufferBaseNV(GL_TRANSFORM_FEEDBACK_BUFFER_NV,0,iso_buffer[0]);
  
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, iso_buffer[1]);
  //??The size of this buffer might need to be larger??
  glBufferDataARB(GL_ARRAY_BUFFER_ARB, mesh*(nz)*(nx)*(ny)*15*4*sizeof(GLfloat),0, GL_STREAM_COPY);
  //glBindBufferBaseNV(GL_TRANSFORM_FEEDBACK_BUFFER_NV,0,iso_buffer[1]);
  
  glBindBuffer(GL_ARRAY_BUFFER_ARB, 0);

  solid = true;
  once = true;

}
void IsoSurface::readInTables(){
  glEnable(GL_DEPTH_TEST);

  std::ifstream fp;
  std::string path;
  path = "case_to_numpolys";
  std::string path2 = "edge_connect_list";

  fp.open(path.c_str());
  if(!fp.is_open()){
    std::cerr << "Couldn't open file: case_to_numpolys";
    exit(1);
  }

  GLfloat* data = new GLfloat[256*4];
  std::string comma;

  for(int i = 0; i < 256; i++){
    int idx = i*4;
    fp >> data[idx];
    fp >> data[idx+1];
    fp >> data[idx+2];
    fp >> data[idx+3];

    //A comma follows each set 
    fp >> comma;
  }
  
  //for(int i = 0; i < 256; i++){
  //std::cout << data[i*4] << " " << 0 << " " << 0 << " " << 0 << " , ";
  //}

  fp.close();

  fp.open(path2.c_str());
  if(!fp.is_open()){
    std::cerr << "Couldn't open file: edge_connect_list";
    exit(1);
  }

  //read in edge list
  GLfloat* data2 = new GLfloat[256*4*5];
  for(int i=0; i < 256; i++){
    for(int j=0; j < 5; j++){
      int idx = i*5*4 + j*4;
      
      fp >> data2[idx];
      fp >> data2[idx+1];
      fp >> data2[idx+2];
      fp >> data2[idx+3];

      //A comma follows each set 
      fp >> comma;

    }
  }
  /*for(int i=0; i < 256; i++){
    for(int j=0; j < 5; j++){
      int idx = i*5*4 + j*4;
      std::cout << data2[idx] << " " << data2[idx+1] << " " << data2[idx+2] << " " << data2[idx+3] << " , ";
    }
    }*/


  fp.close();

  //Create case to numpoly texture
  glEnable(GL_TEXTURE_1D);
  glGenTextures(1,case_to_numpoly);
  
  glBindTexture(GL_TEXTURE_1D,case_to_numpoly[0]);

  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexImage1D(GL_TEXTURE_1D, 0, int_format, 256, 0,GL_RGBA, GL_FLOAT,data); 

  glBindTexture(GL_TEXTURE_1D,0);

  glDisable(GL_TEXTURE_1D);

  //Create edge connect list texture
  glEnable(GL_TEXTURE_2D);
  glGenTextures(1,edge_connect_list);
  
  glBindTexture(GL_TEXTURE_2D,edge_connect_list[0]);

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexImage2D(GL_TEXTURE_2D, 0, int_format, 5,256,0,GL_RGBA, GL_FLOAT,data2); 

  glDisable(GL_TEXTURE_2D);

  delete [] data;
  delete [] data2;

}

void IsoSurface::render3DTexture(FramebufferObject* fbo){
  glEnable(texType2);
  
  fbo->Bind();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);

  glViewport(0, 0, nxdx+1, nydy+1);

  ///////////////////////////////////////////////////
  //Render to 3D texture
  ////////////////////////////////////////////////////
  //Have to render one slice at a time
  render3DShader.activate();
  //glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
  glUniform1fARB(uniform_cValue,contourValue);

  for(int z = 0; z <= nzdz; z++){
    fbo->AttachTexture(GL_COLOR_ATTACHMENT0_EXT,texType2,tex3d[0],0,z);
    float slice = ((float)(z) / (float)(nz));

    glUniform1iARB(u_slice,z);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(texType2,tex3d[1]);
    glUniform1iARB(uniform_tau,0);
    
    glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);

    glBegin(GL_QUADS);
    glTexCoord3f(0.0f,0.0f,slice); glVertex2f(-1.0f, -1.0f);
    glTexCoord3f(1.0f,0.0f,slice); glVertex2f(1.0f, -1.0f);
    glTexCoord3f(1.0f,1.0f,slice); glVertex2f(1.0f, 1.0f);
    glTexCoord3f(0.0f,1.0f,slice); glVertex2f(-1.0f, 1.0f);
    glEnd();

    fbo->Unattach(GL_COLOR_ATTACHMENT0_EXT);

  }
  render3DShader.deactivate();
  
  FramebufferObject::Disable();
  glDisable(texType2);

}
void IsoSurface::createIsoSurface(){

  

  glEnable(GL_RASTERIZER_DISCARD_NV);

  geomShader.activate();
   
  glDisable(texType2);

  glEnable(GL_TEXTURE_1D);
  glActiveTexture(GL_TEXTURE2);
  glBindTexture(GL_TEXTURE_1D,case_to_numpoly[0]);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glUniform1iARB(u_case,2);
  glDisable(GL_TEXTURE_1D);

  glEnable(GL_TEXTURE_2D);
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D,edge_connect_list[0]);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glUniform1iARB(u_edge,1);
  glDisable(GL_TEXTURE_2D);

  glEnable(texType2);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(texType2,tex3d[0]);
  glTexParameteri(texType2, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(texType2, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glUniform1iARB(u_tau3D,0);

  glDisable(texType2);
  
  //Start recording triangle made into a buffer
  glBeginTransformFeedbackNV(GL_TRIANGLES);

  //start query to determine number of triangles made
  glBeginQuery(GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN_NV, query);

  glBindBuffer(GL_ARRAY_BUFFER_ARB, iso_buffer[buffer_num]);
    
  glPointSize(1.0);
  glBegin(GL_POINTS);

  //Visit each voxel in the domain 
  for(int k=0; k < (nzdz*mesh); k++){
    for(int i=0; i < nydy*mesh; i++){
      for(int j=0; j < nxdx*mesh; j++){
	
	glVertex4f((float)j/(float)mesh,(float)i/(float)mesh,
		   (float)k/(float)mesh,1.0);
	
      }
    }
  }
  
  glEnd();
  
  glBindBuffer(GL_ARRAY_BUFFER_ARB,0);

  //Query
  glEndQuery(GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN_NV);

  //End transform feedback and deactivate shader
  glEndTransformFeedbackNV();
  geomShader.deactivate();
   
  
  glDisable(GL_RASTERIZER_DISCARD_NV);
    
  numPrimitives = 0;
  glGetQueryObjectiv(query,GL_QUERY_RESULT,&numPrimitives);
  std::cout << numPrimitives << std::endl;


}

void IsoSurface::draw(){
  glDisable(texType2);
  
  glEnableClientState(GL_VERTEX_ARRAY);
  glBindBufferARB(GL_ARRAY_BUFFER, iso_buffer[buffer_num]);
  glVertexPointer(4,GL_FLOAT,0,0);
  
  if(!solid)
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
  
  //glDrawArrays(GL_TRIANGLES,0,15*(nx)*(ny)*(nz));

  glPushMatrix();
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  //isoShader.activate();
  glColor4f(1.0,0.0,0.0,0.8);
  glDrawArrays(GL_TRIANGLES,0,numPrimitives*3); 
  //isoShader.deactivate();

  glDisable(GL_BLEND);
  glDisable(GL_COLOR_MATERIAL);
  glPopMatrix();

  glBindBuffer(GL_ARRAY_BUFFER,0);
  glDisableClientState(GL_VERTEX_ARRAY);
  ///////////////////////////////////////
  
  glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
}

//These methods don't work yet
void IsoSurface::increaseMesh(){
  if(mesh < 16){
    mesh = mesh*2;    
    buffer_num++;
    glBindBuffer(GL_ARRAY_BUFFER_ARB, iso_buffer[buffer_num]);
    glBindBufferBaseNV(GL_TRANSFORM_FEEDBACK_BUFFER_NV,0,iso_buffer[buffer_num]);
    glBindBuffer(GL_ARRAY_BUFFER_ARB, 0);
  }

 
  
}
void IsoSurface::decreaseMesh(){
  if(mesh > 1){
    mesh = mesh/2;
    buffer_num--;
  }

  
}
