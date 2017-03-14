#ifndef PREDICATES_H
#define PREDICATES_H

#include "ngl/Vec3.h"
#include "Tetrahedron.h"

//----------------------------------------------------------------------------------------------------------------------
/// @file Predicates.h
/// @author Maria Vineeta Bagya Seelan
/// @version 1.0
/// @date 16/08/13
/// Revision History :
/// Initial Version 16/08/13
/// @class Predicates
/// @brief This is the class that defines the two basic geometric tests required
/// @brief for constructing Delaunay Tetraherons
//----------------------------------------------------------------------------------------------------------------------

class Predicates
{
public:
    Predicates();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this method is called when a point p is checked if it is above,
    /// below or on a plane enclosed by points a,b,c
    /// @param _a Point A of the plane
    /// @param _b Point B of the plane
    /// @param _c Point C of the plane
    /// @param _p Point to be checked
    //----------------------------------------------------------------------------------------------------------------------
    float orient3d(ngl::Vec3 _a, ngl::Vec3 _b, ngl::Vec3 _c, ngl::Vec3 _p);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this method is called when a point p is checked if it is within,
    /// outside or on a sphere enclosed by points a,b,c,d
    /// @param _a Point A of the sphere
    /// @param _b Point B of the sphere
    /// @param _c Point C of the sphere
    /// @param _d Point D of the sphere
    /// @param _p Point to be checked
    //----------------------------------------------------------------------------------------------------------------------
    float insphere(ngl::Vec3 _a, ngl::Vec3 _b, ngl::Vec3 _c, ngl::Vec3 _d, ngl::Vec3 _p);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this method is called when a point p is checked if it is within,
    /// outside or on a sphere enclosed by a Tetrahedron
    /// @param _t defines a Tetrahedron with four vertices
    /// @param _point Point to be checked
    //----------------------------------------------------------------------------------------------------------------------
    float insphere3d(Tetrahedron *_t, ngl::Vec3 _point);
    //----------------------------------------------------------------------------------------------------------------------

private :
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief a function that initializes the variables used in Predicates class
    //----------------------------------------------------------------------------------------------------------------------
    void exactinit();
    //----------------------------------------------------------------------------------------------------------------------

};

#endif // PREDICATES_H
