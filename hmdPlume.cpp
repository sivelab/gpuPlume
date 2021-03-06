/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
#include "plumeControl.h"
#include "nonGaussianModel.h"
#include "GaussianModel.h"
#include "Gaussian_2shaders_Model.h"
#include "ReflectionModel.h"
#include "emitterTransform.h"

#include <assert.h>

#include <osg/LightSource>
#include <osg/Light>
#include <osg/MatrixTransform>
#include <osg/Material>
#include <osg/ShapeDrawable>
#include <osgDB/ReadFile>

#include "vissimViewer.h"
#include "vissimViewerEventHandler.h"
#include "MatrixEuler.h"

#include "Settings.h"

bool firstTime;
PlumeControl* plume;
Util* util;

osg::TessellationHints* hints;
emitterTransform* peTransform[10];

int peNum = 0;

class DrawableDrawCallback : public osg::Drawable::DrawCallback
{
  

  virtual void drawImplementation(osg::State& state,const osg::Drawable* drawable) const
  {	  	  
      glPushAttrib(GL_ALL_ATTRIB_BITS);

      glMatrixMode(GL_PROJECTION);
      glPushMatrix();
      glMatrixMode(GL_MODELVIEW);
      glPushMatrix();

      if(firstTime){
         
	glewInit();
	glEnable(GL_DEPTH_TEST);     
	plume->init(true);
	firstTime = false;
        srand48( 4 );
      }		 

      glDisable(GL_LIGHTING);
      int i = plume->display();     
    
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      glPopMatrix();
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
      glPopMatrix();

      glPopAttrib();
  }  

};
class CBoxDrawCallback : public osg::Drawable::DrawCallback
{
  
  virtual void drawImplementation(osg::State& state,const osg::Drawable* drawable) const
  {	  	  
      GLfloat* mvm = new GLfloat[16];

      glPushAttrib(GL_ALL_ATTRIB_BITS);
      
      //Get  Eye Position
      glGetFloatv(GL_MODELVIEW_MATRIX,mvm);

      //std::cout << mvm[12] << " " << mvm[13] << " " << mvm[14] << std::endl;
      float x = mvm[12];//(mvm[12]-plume->nx/2);
      float y = mvm[13];//(mvm[13]-plume->ny/2);
      float z = mvm[14];

      plume->cBoxes[0]->sort(x,y,z);
      plume->cBoxes[0]->draw(plume->sim->curr_timeStep);

      
      glPopAttrib();
  }  

};
osg::Node* emitter(){
  osg::Geode* source = new osg::Geode;
  osg::ShapeDrawable* p;

  if(util->radius[peNum] == 0){
    p = new osg::ShapeDrawable(
	  new osg::Box(osg::Vec3(util->xpos[peNum],util->ypos[peNum],util->zpos[peNum]),0.25),hints);
  }
  else{
    p = new osg::ShapeDrawable(
	  new osg::Sphere(osg::Vec3(util->xpos[peNum],util->ypos[peNum],
                                    util->zpos[peNum]),util->radius[peNum]),hints);	
  
  }
  p->setColor(osg::Vec4(0.0, 0.0, 1.0, 1.0));
  
  source->addDrawable(p);

  osg::MatrixTransform* transform = new osg::MatrixTransform();
  transform->setUpdateCallback(peTransform[peNum]);
  transform->addChild(source);

  peNum++;
  
  return transform;
 
}

// This struct will be included in the packets sent across the
// network.  Do not use pointers in this struct.
typedef struct
{
    float x,y,z;
    //bool releaseOne;
    //bool releasePerSecond;
    particleReleaseType releaseType;
    bool emit;
    bool clear;
    int curr;
    bool keystroke;
    float time_step;
    //int ID;  // store stuff in here that you want sent to the slave...
} sharedDataStruct;


// If MASTER:
//    - the data pointed to by sharedData will be copied into the
//      packets that are broadcast.
// If SLAVE:
//    - the data sent from the master will automatically be at the
//      location of this pointer.
sharedDataStruct* sharedData;


/* An event handler that is used all of the time regardless of Viewer
   and CameraManipulator.  */
class MyEventHandler : public osgGA::GUIEventHandler
{
  protected:
    vissimViewer* _viewer;
    osg::MatrixTransform* _modelNode;
    osg::Group* _rootNode;

    //current particle emitter selected
    int curr;
    bool keystroke;
    bool clear;
    

  public:
    MyEventHandler(vissimViewer* viewer,
                   osg::Group* rootNode,
                   osg::MatrixTransform* modelNode,
                   osg::MatrixTransform* hmdNode)
    {
	assert(viewer);
	_viewer = viewer;

	_modelNode = modelNode;
        _rootNode = rootNode;

        curr = 0;
        keystroke = false;
        clear = false;
   
        // This checks to see if the mode of operation is the SENDER (aka master)
        if(_viewer->m_settings->getNetMode() != RECEIVE)
        {
            // do things only the master would need to do upon initialization...
            // and then copy into (or set) in the shared data structure
            //sharedData->ID = 2;
            
            sharedData->emit = false;
            //sharedData->releaseOne = false;
            //sharedData->releasePerSecond = true;
            sharedData->releaseType = perSecond;
            sharedData->curr = curr;
            sharedData->clear = false;

        }
        else
        {
            sharedData->time_step = 0.0012;
            sharedData->keystroke = false;
            sharedData->curr = 0;
            sharedData->emit = false;
            sharedData->x = util->xpos[0];
            sharedData->y = util->ypos[0];
            sharedData->z = util->zpos[0];
            sharedData->clear = false;
            /* Necessary for the slave not to crash ??? (PeteW) */
            // sharedData->ID = 100000;
            // strcpy(sharedData->targetName, TargetPlacement::getRandomTargetName().c_str());
        }
    }

    /* Handle events -- return true if a keystroke is handled, false
       otherwise.  This function will be called continuously so you can
       use it to make animations, etc. */
    virtual bool handle(const osgGA::GUIEventAdapter& ea,
			osgGA::GUIActionAdapter&)
    {
	//keystroke = false;
        // read out the data from the shared data structure... note
        // that both the master and slave could do things here...
        // This acts like your keyboard handler...

	/* --- handle KEYSTROKES --- */
	if (ea.getEventType() == osgGA::GUIEventAdapter::KEYDOWN)
	{
            if(ea.getKey() == 'p')
            {
                for(int i=0; i<2; i++)
                {
                    osg::Matrix mat;
                    if(i==0)
                    {
                        printf("=== Tracker coordinate frame ===\n");
                        mat = _viewer->vvEventHandler->getMatrix();
                    }
                    else
                    {
                        printf("=== Model coordinate frame ===\n");
                        printf("Tracker frame with model matrix applied (from settings system)\n");
                        mat = _viewer->vvEventHandler->getMatrix() * osg::Matrix::inverse(_viewer->m_settings->getCurrentModelMatrix()) ;
                    }

                    osg::Vec3 pos = mat.getTrans();
                    printf("x/y/z: %0.4f %0.4f %0.4f\n", pos.x(), pos.y(), pos.z());

                    osg::Vec3 orient = MatrixEuler::matrix2euler(mat);
                    printf("yaw/pitch/roll: %0.4f %0.4f %0.4f\n", orient.x(), orient.y(), orient.z());

                    printf("matrix = [\n");
                    for(int i=0; i<4; i++)
                    {
                        for(int j=0; j<4; j++)
                            printf("%0.4f ", mat(j,i));
                        printf("\n");
                    }
                    printf("]\n");
                }
            }

            if(ea.getKey() == ' ' || ea.getKey() == 't')
            {
                // space is hit... increment our ID... (as an example!)
                //sharedData->ID++;
                //return true;
                keystroke = true;
            }
            
            if( ea.getKey() == 'e')
            {
                plume->pe[curr]->emit = !plume->pe[curr]->emit;
                keystroke = true;
                //return true;
            }
            else if( ea.getKey() == 'o')
            {
                //plume->pe[curr]->releasePerSecond = false;
                //plume->pe[curr]->releasePerTimeStep = false;
                //plume->pe[curr]->releaseOne = true;
                plume->pe[curr]->releaseType = onePerKeyPress;
                keystroke = true;
                //return true;
            }
            else if(ea.getKey() == 'i')
            {
                //plume->pe[curr]->releaseOne = false;
                //plume->pe[curr]->releasePerTimeStep = false;
                //plume->pe[curr]->releasePerSecond = true;
                plume->pe[curr]->releaseType = perSecond;
                keystroke = true;
                //return true;
            }
            else if (ea.getKey() == 'W')
            {
                plume->pe[curr]->setPosition(plume->pe[curr]->xpos, plume->pe[curr]->ypos-1.0, plume->pe[curr]->zpos);
                peTransform[curr]->setPosition(plume->pe[curr]->xpos, plume->pe[curr]->ypos, plume->pe[curr]->zpos);
                keystroke = true;
                //return true;
            }
            else if (ea.getKey() == 'S')
            {
                plume->pe[curr]->setPosition(plume->pe[curr]->xpos, plume->pe[curr]->ypos+1.0, plume->pe[curr]->zpos);
                peTransform[curr]->setPosition(plume->pe[curr]->xpos, plume->pe[curr]->ypos, plume->pe[curr]->zpos);
                keystroke = true;
                //return true;
            }
            else if (ea.getKey() == 'A')
            {
                plume->pe[curr]->setPosition(plume->pe[curr]->xpos+1.0, plume->pe[curr]->ypos, plume->pe[curr]->zpos);
                peTransform[curr]->setPosition(plume->pe[curr]->xpos, plume->pe[curr]->ypos, plume->pe[curr]->zpos);
                keystroke = true;
                //return true;
            }
            else if (ea.getKey() == 'D')
            {
                plume->pe[curr]->setPosition(plume->pe[curr]->xpos-1.0, plume->pe[curr]->ypos, plume->pe[curr]->zpos);
                peTransform[curr]->setPosition(plume->pe[curr]->xpos, plume->pe[curr]->ypos, plume->pe[curr]->zpos);
                keystroke = true;
                //return true;
            }
            else if (ea.getKey() == 'Q')
            {
                plume->pe[curr]->setPosition(plume->pe[curr]->xpos, plume->pe[curr]->ypos, plume->pe[curr]->zpos+1.0);
                peTransform[curr]->setPosition(plume->pe[curr]->xpos, plume->pe[curr]->ypos, plume->pe[curr]->zpos);
                keystroke = true;
                //return true;
            }
            else if (ea.getKey() == 'E')
            {
                plume->pe[curr]->setPosition(plume->pe[curr]->xpos, plume->pe[curr]->ypos, plume->pe[curr]->zpos-1.0);
                peTransform[curr]->setPosition(plume->pe[curr]->xpos, plume->pe[curr]->ypos, plume->pe[curr]->zpos);
                keystroke = true;
                //return true;
            }
            else if (ea.getKey() == '>')
            {
                curr++;
                if(curr == util->numOfPE)
                    curr = 0;

                keystroke = true;
                //return true;
            }
            else if(ea.getKey() == 'c')
            {
                //plume->stream->clear();
                clear = true;
                keystroke = true;
                //return true;

            }
            
	}

        /* --- handle SHARED DATA --- */
	if(_viewer->m_settings->getNetMode() != RECEIVE)   // MASTER
	{
                      
            if(keystroke){
                
                sharedData->curr = curr;
                sharedData->releaseType = plume->pe[curr]->releaseType;
                //sharedData->releaseOne = plume->pe[curr]->releaseOne;
                //sharedData->releasePerSecond = plume->pe[curr]->releasePerSecond;
                if(clear){
                    plume->stream->clear();
                    sharedData->clear = true;
                    clear = false;
                } 
                else sharedData->clear = false;
                
            }
            if(!firstTime){
                sharedData->emit = plume->pe[curr]->emit;
                plume->pe[curr]->getPosition(&sharedData->x,&sharedData->y,&sharedData->z);
            }

            sharedData->time_step = plume->time_step;
            sharedData->keystroke = keystroke;

            //if (sharedData->keystroke == true)
            //  std::cout << "Master: keystroke = " << sharedData->keystroke << std::endl;
	}
        else
            { 
                plume->time_step = sharedData->time_step;
               
                // std::cout << "Slave: " << plume->time_step << std::endl;
                
                keystroke = sharedData->keystroke;
                if (keystroke == true)
                    std::cout << "Slave: keystroke = " << keystroke << std::endl;

                //if(sharedData->keystroke && !firstTime){
                if(!firstTime){
                    plume->pe[sharedData->curr]->emit = sharedData->emit;
                    //plume->pe[sharedData->curr]->releaseOne = sharedData->releaseOne;
                    //plume->pe[sharedData->curr]->releasePerSecond = sharedData->releasePerSecond;
                    plume->pe[sharedData->curr]->releaseType = sharedData->releaseType;
                    plume->pe[sharedData->curr]->setPosition(sharedData->x,sharedData->y,sharedData->z);
                    peTransform[sharedData->curr]->setPosition(sharedData->x,sharedData->y,sharedData->z);

                    if(sharedData->clear){
                        plume->stream->clear();
                        sharedData->clear = false;
                    }

                }
            }
        
        if(keystroke){
            keystroke = false;
            //sharedData->keystroke = false;
            return true;
        }

        return keystroke;
    }

};
osg::Node* floor()
{
    osg::Geode* geode = new osg::Geode();

    osg::Geometry* ground = new osg::Geometry();
    osg::Vec3Array* vertices = new osg::Vec3Array(4);
    
    (*vertices)[0].set(0.0, 0.0, 0.0);
    (*vertices)[1].set(0.0, 60.0, 0.0);
    (*vertices)[2].set(60.0, 60.0, 0.0);
    (*vertices)[3].set(60.0, 0.0, 0.0);
    
    ground->setVertexArray(vertices);

    osg::Vec2Array* tex = new osg::Vec2Array(4);
    (*tex)[0].set(0.0, 0.0);
    (*tex)[1].set(4.0, 0.0);
    (*tex)[2].set(4.0, 4.0);
    (*tex)[3].set(0.0, 4.0);

    ground->setTexCoordArray(0,tex);

    osg::Vec3Array* normals = new osg::Vec3Array;
    normals->push_back(osg::Vec3(0.0f, 0.0f, 1.0f));
    ground->setNormalArray(normals);
    ground->setNormalBinding(osg::Geometry::BIND_OVERALL);

    ground->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::QUADS,0,4));
    
    osg::StateSet* stateset = new osg::StateSet();

    osg::Image* image = osgDB::readImageFile( "concrete.jpeg");
    if(image == NULL){
	std::cout << "Couldn't load image" << std::endl;
    }
    osg::Texture2D* texture = new osg::Texture2D;
    //texture->setDataVariance(osg::Object::DYNAMIC); // protect from being optimized away as static state.
    texture->setFilter(osg::Texture2D::MIN_FILTER,osg::Texture2D::LINEAR);
    texture->setFilter(osg::Texture2D::MAG_FILTER,osg::Texture2D::LINEAR);
    texture->setWrap(osg::Texture2D::WRAP_S,osg::Texture2D::REPEAT);
    texture->setWrap(osg::Texture2D::WRAP_T,osg::Texture2D::REPEAT);
    texture->setImage(image);
    stateset->setTextureAttributeAndModes(0,texture,osg::StateAttribute::ON);

    geode->setStateSet( stateset );

    geode->addDrawable(ground);

    return geode;

}

void setupModel(vissimViewer *viewer, osg::MatrixTransform* modelRoot)
{
    assert(modelRoot);

    // Do your stuff here... creating a node called model that can be
    // added to the root at the end
    osg::Group* model = new osg::Group;

    osg::ShapeDrawable* shape = new osg::ShapeDrawable;
    osg::Geode* geode = new osg::Geode;
    shape->setDrawCallback(new DrawableDrawCallback());
    shape->setUseDisplayList(false);
    shape->setSupportsDisplayList(false);
    
    const osg::BoundingBox* bb = new osg::BoundingBox(0.0, 0.0,0.0, plume->nx,plume->ny,plume->nz);
    shape->setInitialBound(*bb);

    geode->addDrawable(shape);
    
    //Collection box
    
    osg::ShapeDrawable* boxShape = new osg::ShapeDrawable;
    osg::Geode* cboxes = new osg::Geode;
    boxShape->setDrawCallback(new CBoxDrawCallback());
    boxShape->setUseDisplayList(false);
    boxShape->setSupportsDisplayList(false);
    
    //const osg::BoundingBox* cbb = new osg::BoundingBox(0,0,0,plume->bounds[3]-plume->bounds[0],
    //                                                  plume->bounds[4]-plume->bounds[1],plume->bounds[5]-plume->bounds[2]);
    const osg::BoundingBox* cbb = new osg::BoundingBox(util->bounds[0],util->bounds[1],util->bounds[2],util->bounds[3],
                                                      util->bounds[4],util->bounds[5]);
    boxShape->setInitialBound(*cbb);

    cboxes->addDrawable(boxShape);
    ////////////////////////////
   
    osg::StateSet* modelStateSet = new osg::StateSet();
    model->setStateSet(modelStateSet);
    
    model->addChild(geode);

    osg::Group* PEGroup = new osg::Group;
    osg::StateSet* peState = new osg::StateSet();
    PEGroup->setStateSet(peState);

    for(int i=0; i < util->numOfPE; i++)
      PEGroup->addChild(emitter());  

    osg::Group* fGroup = new osg::Group;
    fGroup->addChild(floor());

    model->addChild(fGroup);
    model->addChild(PEGroup);
    
    //Add Collection box visual last
    osg::Group* boxGroup = new osg::Group;
    osg::StateSet* boxState = new osg::StateSet();
    boxGroup->setStateSet(boxState);
    //model->addChild(cboxes);
    boxGroup->addChild(cboxes);
    model->addChild(boxGroup);

    modelRoot->addChild(model);
   
}

void setupWand(osg::MatrixTransform* wand)
{
    assert(wand);

    // To get wand location/orientation:  wand->getMatrix()

    // To make wand visible, add geometry under wand node:

    osg::Geode* geode = new osg::Geode;

    // draw a small sphere
    geode->addDrawable(new osg::ShapeDrawable(new osg::Sphere(osg::Vec3(0.0f,0.0f,0.0f), 0.05)));

    osg::StateSet *stateset = geode->getOrCreateStateSet();
    osg::Material* mat = new osg::Material();
    mat->setDiffuse( osg::Material::FRONT, osg::Vec4( 1.0, 0.0, 0.0, 1.0 ) );
    stateset->setAttribute( mat );
    wand->addChild(geode);

}

int main( int argc, char **argv )
{

     //This is my plume initialization
    util = new Util();
    util->readInput("Settings/input.txt");
    switch(util->advectChoice){
    case 0:
        //util->windFieldData = 4;
        plume = new Gaussian_2shaders_Model(util);
        break;
    case 1:
        //util->windFieldData = 4;
        plume = new GaussianModel(util);
        break;
    case 2:
        util->windFieldData = 5;
        plume = new NonGaussianModel(util);
        break;
    case 3:
        plume = new ReflectionModel(util);
        break;
    default:
        std::cout << "Error in advection Choice in Settings file!" << std::endl;
    }

    //plume = new PlumeControl();
    firstTime = true;

    

    for(int i=0; i < util->numOfPE; i++){
      peTransform[i] = new emitterTransform(util->xpos[i],util->ypos[i],util->zpos[i]);   
    } 
   
    // use an ArgumentParser object to manage the program arguments.
    osg::ArgumentParser arguments(&argc,argv);

    // set up the usage document, in case we need to print out how to use this program.
    arguments.getApplicationUsage()->setApplicationName(arguments.getApplicationName());
    arguments.getApplicationUsage()->setDescription(arguments.getApplicationName()+" is an example program that loads and visualizes 3 models.");
    arguments.getApplicationUsage()->setCommandLineUsage(arguments.getApplicationName()+" [options] filename ...");
    arguments.getApplicationUsage()->addCommandLineOption("-f <filename>","File to view.");

    // The vissim viewer should initially be set to 0.  We'll create the
    // exact viewer based on the arguments, depending on whether we're
    // running in standalone, master, or slave mode.
    vissimViewer *viewer = 0;

    viewer = new vissimViewer( arguments );


    // We must have a viewer by this point in the code.  If not, fail.
    if(viewer == NULL)
    {
        osg::notify(osg::FATAL) << "Could not create a viewer." << std::endl;
        exit(1);
    }


    // Initialize the viewer components: the multicast network sync'ing,
    // the treadport controller manipulator, the HMD head tracker
    // manipulator, and anything else that needs to be setup prior to
    // starting the simulation loop.
    viewer->initialize();

    // allocate space for the shared data.
    sharedData = (sharedDataStruct*) malloc(sizeof(sharedDataStruct));

    // if this process is SLAVE, this tells the viewer to put shared
    // data received across the network in the area pointed to by
    // sharedData.  If this process is MASTER, this tells the viewer to
    // send the data pointed to by sharedData across the network
    // whenever a packet is sent.
    viewer->setSharedData(sharedData, sizeof(sharedDataStruct));


    // Get the ROOT node for our scene.
    osg::Group* rootNode = viewer->getRootNode();

    osg::MatrixTransform* rootTransform = new osg::MatrixTransform();
    rootTransform->setMatrix(osg::Matrix::translate(osg::Vec3(-plume->nx/2,-plume->ny/2,0.0)));
    

    // Create a MatrixTransform node to put the geometry under.
    osg::MatrixTransform* modelNode = new osg::MatrixTransform();

    // Don't mess with this yet - Pete
    // get a wand node.  Include any geometry associated with the wand
    // under this node.  The matrix in wandNode will be automatically
    // updated to show the wand moving in the scene.  If you just want
    // the location of wand (no geometry), do not attach this node to
    // the root node.  Then, use wandNode->getMatrix() to get location
    // information.
    osg::MatrixTransform* hmdNode = viewer->getStationNode(viewer->m_settings->getIsenseStationHMD());
    osg::MatrixTransform* wandNode = viewer->getStationNode(viewer->m_settings->getIsenseStationWand());


    // Make some assertions
    assert(rootNode);
    assert(modelNode);
    assert(wandNode);
    assert(hmdNode);

    // Add our nodes to the root node
    rootTransform->addChild(modelNode);

    rootNode->addChild(rootTransform);
    rootNode->addChild(wandNode);  
    // (we don't actually want to add the hmdNode to the rootNode, we
    // just want to have access to it so we can know where the HMD is.)
    
    hints = new osg::TessellationHints;
    hints->setDetailRatio(1.0f);

    

    // Attach our geometry to the nodes ... this is where you can stick your geometry.... look at
    // the setup function...
    setupModel(viewer, modelNode);

    // Register our event handler with the viewer.
    viewer->getEventHandlerList().push_front(new MyEventHandler(viewer, rootNode, modelNode, hmdNode));
 
    viewer->optimize();

    while( !viewer->done() && !viewer->master_killed() )
    {
        //
        // Wait for all cull and draw threads to complete.
        //
        viewer->sync();

        //
        // Update the scene by traversing it with the the update visitor
        // which will call all node update callbacks and animations.
        // This call also takes care of the synchronization between the
        // various rendering systems-- see the Slave and Master update
        // functions.
        //
        viewer->update();

        //
        // Fire off the cull and draw traversals of the scene.
        //
        if (!viewer->master_killed())
            viewer->frame();
    }

    //
    // Upon completion/master_killed, wait for all cull and draw threads
    // to complete before exit.
    //
    viewer->sync();

    //
    // Perform any cleanup in the post update loop cleanup function.
    //
    viewer->cleanup();
    osg::notify( osg::WARN ) << "Renderer exited normally" << std::endl;

    delete viewer;

    return 0;
}

