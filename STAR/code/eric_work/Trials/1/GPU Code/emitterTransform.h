#include <iostream>
#include <osg/Transform>
#include <osg/MatrixTransform>

class emitterTransform : public osg::NodeCallback
{

    public:
		
        emitterTransform(float x,float y,float z){	 
	  xpos = x;
	  ypos = y;
	  zpos = z;
                    			
        }
	void setPosition(float x,float y, float z){
	    dx = x-xpos;
	    dy = y-ypos;
	    dz = z-zpos;
	}

        virtual void operator() (osg::Node* node, osg::NodeVisitor* nv)
        {
	  osg::MatrixTransform* transform = dynamic_cast<osg::MatrixTransform*>(node);  
	  if (nv && transform && nv->getFrameStamp())
            {	     
		transform->setMatrix(osg::Matrix::translate(osg::Vec3(dx,dy,dz)));	    
	    }
            
            // must continue subgraph traversal.
            traverse(node,nv);            
            
        }
        
    protected:
    
        float dx,dy,dz,xpos,ypos,zpos;
	

};
