#ifndef RENDERER_H
#define RENDERER_H

//----------------------------------------------------------------------------------------------------------------------
/// @file Renderer.h
/// @brief a basic Qt class for Rendering
/// @author Maria Vineeta Bagya Seelan
/// @version 1.0
/// @date 16/08/13
/// Revision History :
/// Initial Version 16/08/13
/// @class Renderer
/// @brief our main glwindow widget for the application and all drawing elements are
/// included in this file
//----------------------------------------------------------------------------------------------------------------------

#include <QEvent>
#include <QResizeEvent>
#include "ngl/Light.h"
#include "ngl/NGLInit.h"
#include "ngl/Obj.h"
#include "ngl/Camera.h"
#include "ngl/TransformStack.h"
#include "MeshSampler.h"

class Renderer : public QGLWidget
{
Q_OBJECT
public :
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Constructor for GLWindow
    /// @param [in] _parent the parent window to create the GL context in
    //----------------------------------------------------------------------------------------------------------------------
    Renderer(
              QWidget *_parent
            );
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief dtor
    //----------------------------------------------------------------------------------------------------------------------
    ~Renderer();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief sets the obj file received from UI
    /// @param [in] _filename stores the filename
    //----------------------------------------------------------------------------------------------------------------------
    void setObjFilename( std::string _filename );

public slots :
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief updates the density of points
    /// @parm[in] _density the value to update
    //----------------------------------------------------------------------------------------------------------------------
    void updateDensity( int _density );
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief slot to check if Compute button under Delaunay is clicked
    //----------------------------------------------------------------------------------------------------------------------
    void isDelaunayComputeClicked();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief slot to check if Compute button under Voronoi is clicked
    //----------------------------------------------------------------------------------------------------------------------
    void isVoronoiComputeClicked();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief slot to check if the volume points checkbox is activated
    /// @parm[in] _isCheckBoxEnabled the value to update
    //----------------------------------------------------------------------------------------------------------------------
    void isPointsActivated(int _isCheckBoxEnabled );
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief slot to check if the Surface points checkbox is activated
    /// @parm[in] _isCheckBoxEnabled the value to update
    //----------------------------------------------------------------------------------------------------------------------
    void isSurfacePointsActivated(int _isCheckBoxEnabled );
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief slot to check if the Intersection points checkbox is activated
    /// @parm[in] _isCheckBoxEnabled the value to update
    //----------------------------------------------------------------------------------------------------------------------
    void isHitPointsActivated(int _isCheckBoxEnabled );
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief slot to check if Clear button is clicked
    //----------------------------------------------------------------------------------------------------------------------
    void clear();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief slot to check if Clear button under Delaunay is clicked
    //----------------------------------------------------------------------------------------------------------------------
    void clearDelaunay();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief slot to check if Clear button under Voronoi is clicked
    //----------------------------------------------------------------------------------------------------------------------
    void clearVoronoi();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief updates the method the samples have to be generated
    /// @parm[in] _meshSampling the value to update
    //----------------------------------------------------------------------------------------------------------------------
    void updateMethod( int _method );
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief slot to check if the Wireframe checkbox under Delaunay is activated
    /// @parm[in] _isCheckBoxEnabled the value to update
    //----------------------------------------------------------------------------------------------------------------------
    void isDelaunayWireframeClicked( int _isCheckBoxEnabled );
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief updates the point location the samples have to be generated
    /// @parm[in] _location the value set
    //----------------------------------------------------------------------------------------------------------------------
    void pointLocation( int _location );

private :
    /// @brief used to store the x rotation mouse value
    int m_spinXFace;
    /// @brief used to store the y rotation mouse value
    int m_spinYFace;
    /// @brief flag to indicate if the mouse button is pressed when dragging
    bool m_rotate;
    /// @brief flag to indicate if the Right mouse button is pressed when dragging
    bool m_translate;
    /// @brief the previous x mouse value
    int m_origX;
    /// @brief the previous y mouse value
    int m_origY;
    /// @brief the previous x mouse value for Position changes
    int m_origXPos;
    /// @brief the previous y mouse value for Position changes
    int m_origYPos;
    /// @brief the model position for mouse movement
    ngl::Vec4 m_modelPos;
    /// @brief Our Camera
    ngl::Camera *m_cam;
    /// @brief our transform stack for drawing the teapots
    ngl::TransformStack m_transformStack;
    /// @brief flag to show bounding box
    bool m_showBBox;
    /// @brief flag to show bounding sphere
    bool m_showBSphere;
    /// @brief the mesh with all the data in it
    ngl::VertexArrayObject *m_vaoMesh;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Mesh class that creates the Mesh object
    //----------------------------------------------------------------------------------------------------------------------
    MeshSampler *m_mesh;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief store the value when show/hide points checkbox is activated
    //----------------------------------------------------------------------------------------------------------------------
    bool m_showPoints;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief store the value when surface points checkbox is activated
    //----------------------------------------------------------------------------------------------------------------------
    bool m_isSurfacePoints;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief store the value when intersection points checkbox is activated
    //----------------------------------------------------------------------------------------------------------------------
    bool m_isHitPoints;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief store the value when draw rays checkbox is activated
    //----------------------------------------------------------------------------------------------------------------------
    bool m_isDrawRays;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief store the value when Compute button under Delaunay is clicked
    //----------------------------------------------------------------------------------------------------------------------
    int m_dCompute;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief store the value when Compute button under Voronoi is clicked
    //----------------------------------------------------------------------------------------------------------------------
    int m_vCompute;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief stores the value when Clear button is clicked
    //----------------------------------------------------------------------------------------------------------------------
    bool m_update;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief stores the sample points creation method
    //----------------------------------------------------------------------------------------------------------------------
    int m_method;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief stores the value when Wireframe check box under Delaunay is activated
    //----------------------------------------------------------------------------------------------------------------------
    int m_dWireframe;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief stores the value when Point location is selected
    //----------------------------------------------------------------------------------------------------------------------
    int m_pointType;

protected:
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief  The following methods must be implimented in the sub class
    /// this is called when the window is created
    //----------------------------------------------------------------------------------------------------------------------
    void initializeGL();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this is called whenever the window is re-sized
    /// @param[in] _w the width of the resized window
    /// @param[in] _h the height of the resized window
    //----------------------------------------------------------------------------------------------------------------------
    void resizeGL(
                    const int _w,
                    const int _h
                  );
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this is the main gl drawing routine which is called whenever the window needs to
    /// be re-drawn
    //----------------------------------------------------------------------------------------------------------------------
    void paintGL();

private :

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this method is called every time a mouse is moved
    /// @param _event the Qt Event structure
    //----------------------------------------------------------------------------------------------------------------------
    void mouseMoveEvent (
                          QMouseEvent * _event
                        );
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this method is called everytime the mouse button is pressed
    /// inherited from QObject and overridden here.
    /// @param _event the Qt Event structure
    //----------------------------------------------------------------------------------------------------------------------
    void mousePressEvent (
                            QMouseEvent *_event
                         );

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this method is called everytime the mouse button is released
    /// inherited from QObject and overridden here.
    /// @param _event the Qt Event structure
    //----------------------------------------------------------------------------------------------------------------------
    void mouseReleaseEvent (
                            QMouseEvent *_event
                            );

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this method is called everytime the mouse wheel is moved
    /// inherited from QObject and overridden here.
    /// @param _event the Qt Event structure
    //----------------------------------------------------------------------------------------------------------------------
    void wheelEvent(
                      QWheelEvent *_event
                   );

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this method is called everytime a timer is triggered
    /// inherited from QObject and overridden here.
    /// @param _event the Qt Event structure
    //----------------------------------------------------------------------------------------------------------------------
    void timerEvent(QTimerEvent *_event);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief a function that computes the Delaunay Triangulations
    //----------------------------------------------------------------------------------------------------------------------
    void computeDelaunay();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief a function that calculates the sample points for the mesh
    //----------------------------------------------------------------------------------------------------------------------
    void calculate();
};

#endif // RENDERER_H

//----------------------------------------------------------------------------------------------------------------------
