/* -*-c++-*- OpenSceneGraph - Copyright (C) 1998-2006 Robert Osfield 
 *
 * This application is open source and may be redistributed and/or modified   
 * freely and without restriction, both in commericial and non commericial applications,
 * as long as this copyright notice is maintained.
 * 
 * This application is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
*/
#include <iostream>

#include "plumeControl.h"

#include <osgProducer/Viewer>

#include <osg/Transform>
#include <osg/Billboard>
#include <osg/Geode>
#include <osg/ShapeDrawable>
#include <osg/Group>
#include <osg/Notify>

#include <osgDB/Registry>
#include <osgDB/ReadFile>

#include <osgGA/TrackballManipulator>
#include <osgGA/FlightManipulator>
#include <osgGA/DriveManipulator>

#include <osgUtil/Optimizer>
bool firstTime;
PlumeControl* plume;

class DrawableDrawCallback : public osg::Drawable::DrawCallback
{


        virtual void drawImplementation(osg::State& state,const osg::Drawable* drawable) const
        {	  
	  	  
	  glPushAttrib(GL_ALL_ATTRIB_BITS);

	  glMatrixMode(GL_PROJECTION);
	  glPushMatrix();
	  //glLoadIdentity();
	  glMatrixMode(GL_MODELVIEW);
	  glPushMatrix();
	  //glLoadIdentity();

	  if(firstTime){
	    
	    glewInit();
	    glEnable(GL_DEPTH_TEST);
	    plume->init();	    	    
	    firstTime = false;
	  }		 

	  glDisable(GL_LIGHTING);
	  plume->display();
	  /*glBegin(GL_QUADS);
	  {
	    glColor3f(1.0, 1.0, 0.0);
	    glNormal3f(0.0, 1.0, 0.0);
	    glVertex3f(-100.0, 0.0, 100.0);
	    glVertex3f(100.0, 0.0, 100.0);
	    glVertex3f(100.0, 0.0, -100.0);
	    glVertex3f(-100.0, 0.0, -100.0);

	  }
	  glEnd();*/	 
	  glMatrixMode(GL_PROJECTION);
	  glLoadIdentity();
	  glPopMatrix();
	  glMatrixMode(GL_MODELVIEW);
	  glLoadIdentity();
	  glPopMatrix();

	  glPopAttrib();	  	  
        }

};

int main( int argc, char **argv )
{
    plume = new PlumeControl(128,128,4);
    firstTime = true;
    // use an ArgumentParser object to manage the program arguments.
    osg::ArgumentParser arguments(&argc,argv);
   
    // set up the usage document, in case we need to print out how to use this program.
    arguments.getApplicationUsage()->setDescription(arguments.getApplicationName()+" is the example which demonstrates use of the range of different types of callbacks supported in the OpenSceneGraph.");
    arguments.getApplicationUsage()->setCommandLineUsage(arguments.getApplicationName()+" [options] filename ...");
    arguments.getApplicationUsage()->addCommandLineOption("-h or --help","Display this information");
   
    // initialize the viewer.
    osgProducer::Viewer viewer(arguments);

    // set up the value with sensible default event handlers.
    viewer.setUpViewer(osgProducer::Viewer::STANDARD_SETTINGS);

    // get details on keyboard and mouse bindings used by the viewer.
    viewer.getUsage(*arguments.getApplicationUsage());

    // if user request help write it out to cout.
    if (arguments.read("-h") || arguments.read("--help"))
    {
        arguments.getApplicationUsage()->write(std::cout);
        return 1;
    }

    // any option left unread are converted into errors to write out later.
    arguments.reportRemainingOptionsAsUnrecognized();
    
    // report any errors if they have occured when parsing the program aguments.
    if (arguments.errors())
    {
        arguments.writeErrorMessages(std::cout);
        return 1;
    }
    
    // load the nodes from the commandline arguments.
    //osg::Node* rootnode = osgDB::readNodeFiles(arguments);
    osg::Group* rootnode = new osg::Group;
    osg::ShapeDrawable* shape = new osg::ShapeDrawable;
    osg::Geode* geode = new osg::Geode;
 

    //setVertexProgramSupported(true);
    shape->setDrawCallback(new DrawableDrawCallback());
    shape->setUseDisplayList(false);
    geode->addDrawable(shape);

    rootnode->addChild(geode);

    // set the scene to render
    viewer.setSceneData(rootnode);
   

    // create the windows and run the threads.
    viewer.realize();
    

    while( !viewer.done() )
    {
        // wait for all cull and draw threads to complete.
        viewer.sync();

        // update the scene by traversing it with the the update visitor which will
        // call all node update callbacks and animations.
        viewer.update();
         
        // fire off the cull and draw traversals of the scene.
        viewer.frame();
        
    }
    
    // wait for all cull and draw threads to complete.
    viewer.sync();

    // run a clean up frame to delete all OpenGL objects.
    viewer.cleanup_frame();

    // wait for all the clean up frame to complete.
    viewer.sync();

    return 0;
}

