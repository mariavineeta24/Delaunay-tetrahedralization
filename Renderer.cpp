//----------------------------------------------------------------------------------------------------------------------
/// @file Renderer.cpp
/// @brief Class that draws the GL elements
//----------------------------------------------------------------------------------------------------------------------

#include "Renderer.h"

#include <iostream>
#include <ngl/Colour.h>
#include <ngl/Mat4.h>
#include <ngl/ShaderLib.h>
#include <ngl/Transformation.h>

//----------------------------------------------------------------------------------------------------------------------
/// @brief the increment for x/y translation with mouse movement
//----------------------------------------------------------------------------------------------------------------------
const static  float INCREMENT=0.01;
//----------------------------------------------------------------------------------------------------------------------
/// @brief the increment for the wheel zoom
//----------------------------------------------------------------------------------------------------------------------
const static float ZOOM=0.1;

//----------------------------------------------------------------------------------------------------------------------
Renderer::Renderer(
                   QWidget *_parent
                  )
                    : QGLWidget( new CreateCoreGLContext(QGLFormat::defaultFormat()), _parent )
{

  // set this widget to have the initial keyboard focus
  setFocus();
  // re-size the widget to that of the parent (in this case the GLFrame passed in on construction)
  this->resize(_parent->size());
  // Now set the initial GLWindow attributes to default values
  // Rotate is false
  m_rotate=false;
  // mouse rotation values set to 0
  m_spinXFace=0;
  m_spinYFace=0;
  m_origX=0;
  m_origY=0;
  m_pointType = 0;
  m_mesh = new MeshSampler();
  m_isSurfacePoints = false;
  m_isHitPoints = false;
  m_dCompute = false;
  m_vCompute = false;
  m_update = true;
  m_dWireframe = false;
  m_showPoints = true;
}

//----------------------------------------------------------------------------------------------------------------------
Renderer::~Renderer()
{
  ngl::NGLInit *Init = ngl::NGLInit::instance();
  std::cout<<"Shutting down NGL, removing VAO's and Shaders\n";
  Init->NGLQuit();
  delete m_mesh;
}

//----------------------------------------------------------------------------------------------------------------------
// This virtual function is called once before the first call to paintGL() or resizeGL(),
//and then once whenever the widget has been assigned a new QGLContext.
// This function should set up any required OpenGL context rendering flags, defining VBOs etc.
//----------------------------------------------------------------------------------------------------------------------
void Renderer::initializeGL()
{
  // Grey Background
  glClearColor(0.4f, 0.4f, 0.4f, 1.0f);
  // enable depth testing for drawing
  glEnable(GL_DEPTH_TEST);
  // we need to initialise the NGL lib, under windows and linux we also need to
  // initialise GLEW, under windows this needs to be done in the app as well
  // as the lib hence the WIN32 define
  ngl::NGLInit *Init = ngl::NGLInit::instance();
  // init GLEW for the lib
  Init->initGlew();
  #ifdef WIN32
    glewInit(); // need a local glew init as well as lib one for windows
  #endif

    // Now we will create a basic Camera from the graphics library
    // This is a static camera so it only needs to be set once
    // First create Values for the camera position
    ngl::Vec3 from(0,0,4);
    ngl::Vec3 to(0,0,0);
    ngl::Vec3 up(0,1,0);
    m_cam= new ngl::Camera(from,to,up,ngl::PERSPECTIVE);
    // set the shape using FOV 45 Aspect Ratio based on Width and Height
    // The final two are near and far clipping planes of 0.5 and 10
    m_cam->setShape(45,(float)720.0/576.0,0.05,350,ngl::PERSPECTIVE);

    // now to load the shader and set the values
    // grab an instance of shader manager
    ngl::ShaderLib *shader=ngl::ShaderLib::instance();
    // load a frag and vert shaders
    glEnable(GL_DEPTH_TEST); // for removal of hidden surfaces

    // Color shader
    shader->createShaderProgram("Colour");

    shader->attachShader("ColourVertex",ngl::VERTEX);
    shader->attachShader("ColourFragment",ngl::FRAGMENT);
    shader->loadShaderSource("ColourVertex","shaders/Colour.vs");
    shader->loadShaderSource("ColourFragment","shaders/Colour.fs");

    shader->compileShader("ColourVertex");
    shader->compileShader("ColourFragment");
    shader->attachShaderToProgram("Colour","ColourVertex");
    shader->attachShaderToProgram("Colour","ColourFragment");

    shader->bindAttribute("Colour",0,"inVert");
    shader->linkProgramObject("Colour");

    // load the mesh
    m_mesh->loadMesh();
    // calculate, calls the sampling function
    calculate();
}

//----------------------------------------------------------------------------------------------------------------------
//This virtual function is called whenever the widget has been resized.
// The new size is passed in width and height.
//----------------------------------------------------------------------------------------------------------------------
void Renderer::resizeGL(
                        int _w,
                        int _h
                       )
{
    // set the viewport for openGL
    glViewport(0,0,_w,_h);
    m_cam->setShape(45,(float)_w/_h,0.05,350,ngl::PERSPECTIVE);
}

//----------------------------------------------------------------------------------------------------------------------
//This virtual function is called whenever the widget needs to be painted.
// this is our main drawing routine
//----------------------------------------------------------------------------------------------------------------------
void Renderer::paintGL()
{
    // clear the screen and depth buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    // Rotation based on the mouse position for our global transform
    ngl::Transformation trans;
    ngl::Mat4 rotX;
    ngl::Mat4 rotY;
    // create the rotation matrices
    rotX.rotateX(m_spinXFace);
    rotY.rotateY(m_spinYFace);
    // multiply the rotations
    ngl::Mat4 final=rotY*rotX;
    // add the translations
    final.m_m[3][0] = m_modelPos.m_x;
    final.m_m[3][1] = m_modelPos.m_y;
    final.m_m[3][2] = m_modelPos.m_z;
    // set this in the TX stack
    trans.setMatrix(final);
    m_transformStack.setGlobal(trans);

    if(m_update)
    {
        m_mesh->drawMesh(m_transformStack,m_cam);
        m_mesh->drawBBox(m_transformStack,m_cam);
        if(m_showPoints)
        {
            if(m_pointType != 0)
            {
                m_mesh->drawSurfaceMeshPoints(m_transformStack,m_cam);
            }
            else
            {
                m_mesh->drawVolumePoints(m_transformStack,m_cam);
            }
        }
        if(m_isSurfacePoints)
        {
            m_mesh->drawSurfacePointsBBox(m_transformStack,m_cam);
        }
        if(m_isHitPoints)
        {
            m_mesh->drawHitPoints(m_transformStack,m_cam);
        }

        if(m_dWireframe)
        {
            glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
        }
        else
        {
            glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
        }
        if(m_dCompute)
        {
            m_mesh->drawTetrahedron(m_transformStack,m_cam);
        }
        if(m_vCompute)
        {
            m_mesh->drawVoronoi(m_transformStack,m_cam);
            m_mesh->drawVoronoiVertices(m_transformStack,m_cam);
        }

    }

}

//----------------------------------------------------------------------------------------------------------------------
void Renderer::mouseMoveEvent (
                               QMouseEvent * _event
                              )
{
  // note the method buttons() is the button state when event was called
  // this is different from button() which is used to check which button was
  // pressed when the mousePress/Release event is generated
  if(m_rotate && _event->buttons() == Qt::LeftButton)
  {
    m_spinYFace = ( m_spinYFace + (_event->x() - m_origX) ) % 360 ;
    m_spinXFace = ( m_spinXFace + (_event->y() - m_origY) ) % 360 ;
    m_origX = _event->x();
    m_origY = _event->y();
    updateGL();

  }
  // right mouse translate code
  else if(m_translate && _event->buttons() == Qt::RightButton)
  {
    int diffX = (int)(_event->x() - m_origXPos);
    int diffY = (int)(_event->y() - m_origYPos);
    m_origXPos=_event->x();
    m_origYPos=_event->y();
    m_modelPos.m_x += INCREMENT * diffX;
    m_modelPos.m_y -= INCREMENT * diffY;
    updateGL();
  }
}

//----------------------------------------------------------------------------------------------------------------------
void Renderer::mousePressEvent (
                                QMouseEvent * _event
                               )
{
  // this method is called when the mouse button is pressed in this case we
  // store the value where the maouse was clicked (x,y) and set the Rotate flag to true
  if(_event->button() == Qt::LeftButton)
  {
    m_origX = _event->x();
    m_origY = _event->y();
    m_rotate =true;
  }
  // right mouse translate mode
  else if(_event->button() == Qt::RightButton)
  {
    m_origXPos = _event->x();
    m_origYPos = _event->y();
    m_translate=true;
  }

}

//----------------------------------------------------------------------------------------------------------------------
void Renderer::mouseReleaseEvent (
                                  QMouseEvent * _event
                                 )
{
  // this event is called when the mouse button is released
  // we then set Rotate to false
  if (_event->button() == Qt::LeftButton)
  {
    m_rotate=false;
  }
  // right mouse translate mode
  if (_event->button() == Qt::RightButton)
  {
    m_translate=false;
  }
}

//----------------------------------------------------------------------------------------------------------------------
void Renderer::wheelEvent(QWheelEvent *_event)
{

    // check the diff of the wheel position (0 means no change)
    if(_event->delta() > 0)
    {
        m_modelPos.m_z+=ZOOM;
    }
    else if(_event->delta() < 0 )
    {
        m_modelPos.m_z-=ZOOM;
    }
    updateGL();
}

//----------------------------------------------------------------------------------------------------------------------
void Renderer::timerEvent(
                           QTimerEvent *_event
                         )
{
  Q_UNUSED(_event);
}

//----------------------------------------------------------------------------------------------------------------------
void Renderer::setObjFilename(
                               std::string _filename
                             )
{
    m_mesh->setObjFilename(_filename);
    m_mesh->loadMesh();
    calculate();
    m_dCompute = false;
    updateGL();
}

//----------------------------------------------------------------------------------------------------------------------
void Renderer::isPointsActivated( int _isCheckBoxEnabled )
{
    if(_isCheckBoxEnabled == 0)
    {
        m_showPoints = false;
    }
    else
    {
        m_showPoints = true;
    }
    updateGL();
}

//----------------------------------------------------------------------------------------------------------------------
void Renderer::isSurfacePointsActivated( int _isCheckBoxEnabled )
{
    if(_isCheckBoxEnabled == 0)
    {
        m_isSurfacePoints = false;
    }
    else
    {
        m_isSurfacePoints = true;
    }
    updateGL();
}

//----------------------------------------------------------------------------------------------------------------------
void Renderer::isHitPointsActivated( int _isCheckBoxEnabled  )
{
    if(_isCheckBoxEnabled == 0)
    {
        m_isHitPoints = false;
    }
    else
    {
        m_isHitPoints = true;
    }
    updateGL();
}

//----------------------------------------------------------------------------------------------------------------------
void Renderer::updateMethod(int _method)
{
    m_method = _method;
    m_mesh->setMethod(_method);
    m_dCompute = false;
    calculate();
    updateGL();
}

//----------------------------------------------------------------------------------------------------------------------
void Renderer::isDelaunayComputeClicked()
{
    m_dCompute = true;
    computeDelaunay();
    updateGL();
}

//----------------------------------------------------------------------------------------------------------------------
void Renderer::isVoronoiComputeClicked()
{
    m_vCompute = true;
    updateGL();
}

//----------------------------------------------------------------------------------------------------------------------
void Renderer::pointLocation( int _location )
{
    m_pointType = _location;
    m_mesh->setPointLocationType(_location);
    m_dCompute = false;
    calculate();
    updateGL();
}

//----------------------------------------------------------------------------------------------------------------------
void Renderer::isDelaunayWireframeClicked( int _isCheckBoxEnabled )
{
    if(_isCheckBoxEnabled == 0)
    {
        m_dWireframe = false;
    }
    else
    {
        m_dWireframe = true;
    }
    updateGL();
}

//----------------------------------------------------------------------------------------------------------------------
void Renderer::updateDensity(int _density)
{
    m_mesh->setDensity(_density);
    calculate();
    m_dCompute = false;
}

//----------------------------------------------------------------------------------------------------------------------
void Renderer::computeDelaunay()
{
    m_mesh->delaunay();
}

//----------------------------------------------------------------------------------------------------------------------
void Renderer::calculate()
{
   m_mesh->SampleMesh();
}

//----------------------------------------------------------------------------------------------------------------------
void Renderer::clearDelaunay()
{
    m_dCompute = false;
    updateGL();
}

//----------------------------------------------------------------------------------------------------------------------
void Renderer::clearVoronoi()
{
    m_vCompute = false;
    updateGL();
}

//----------------------------------------------------------------------------------------------------------------------
void Renderer::clear()
{
    m_update = false;
    updateGL();
}
//----------------------------------------------------------------------------------------------------------------------

