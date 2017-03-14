#ifndef DELAUNAY_H
#define DELAUNAY_H

//----------------------------------------------------------------------------------------------------------------------
/// @file Delaunay.h
/// @author Maria Vineeta Bagya Seelan
/// @version 1.0
/// @date 28/03/13
/// Revision History :
/// Initial Version 28/03/13
/// @class Delaunay
/// @brief the Delaunay class that holds all members and functions needed to compute Delaunay Tetrahedralization
//----------------------------------------------------------------------------------------------------------------------

#include "Tetrahedron.h"
#include "Point3.h"
#include "Predicates.h"
#include "Voronoi.h"
#include <stack>

//----------------------------------------------------------------------------------------------------------------------
/// @brief structure that stores the data required for flipping
//----------------------------------------------------------------------------------------------------------------------
struct flipData
{
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the point to be inserted
    //----------------------------------------------------------------------------------------------------------------------
    int m_ptPos;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the tetrahedron that has the point p
    //----------------------------------------------------------------------------------------------------------------------
    Tetrahedron* m_t;
};

class Delaunay
{
public:
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Constructor for Delaunay class
    /// @param [in] _tetrahedron the tetrahedron that has the point p
    //----------------------------------------------------------------------------------------------------------------------
    Delaunay( Tetrahedron* _tetrahedron );
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Destructor for Delaunay class
    //----------------------------------------------------------------------------------------------------------------------
    ~Delaunay();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief main function that does the delaunay computations
    /// @param [in] _points all the points needed to construct DT
    /// @param [out] std::vector<Tetrahedron*> the tetrahedrons used in constructing DT
    //----------------------------------------------------------------------------------------------------------------------
    std::vector<Tetrahedron*> compute( std::vector<ngl::Vec3> _points );
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function that generates Voronoi diagram
    /// @param [in] _point single point around which Voronoi cells are constructed
    /// @param [out] Voronoi the voronoi cell after construction
    //----------------------------------------------------------------------------------------------------------------------
    Voronoi *generateVoronoi(ngl::Vec3 _point);

private :
    std::stack<flipData> m_flipStack;
    std::vector<Tetrahedron*> m_tetrahedra;
    Tetrahedron* m_tetrahedron;
    Tetrahedron* m_third;
    Predicates *m_predicates;
    int m_a;
    int m_b;
    int m_cid;
    int m_tetCount;

private :
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function that performs flip14
    /// @param [in] _tetra the tetrahedra to be flipped
    /// @param [in] _point the point to be inserted into DT
    /// @param [out] returns the new tetrahedra created
    //----------------------------------------------------------------------------------------------------------------------
    Tetrahedron* flip14(Tetrahedron* _tetra, ngl::Vec3 _point);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function that performs flip23
    /// @param [in] _t the tetrahedra that contains point p
    /// @param [in] _pid position of the point in _t
    /// @param [in] _did position of the apex in tetrahedron adjacent to _t
    /// @param [out] returns the new tetrahedra created
    //----------------------------------------------------------------------------------------------------------------------
    Tetrahedron* flip23(Tetrahedron* _t, int _pid, int _did);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function that performs flip14
    /// @param [in] _t1 the tetrahedra to be flipped
    /// @param [in] _t2 the tetrahedra to be flipped
    /// @param [in] _t3 the tetrahedra to be flipped
    /// @param [in] _pid position of the point in _t1
    /// @param [in] _did position of the apex in tetrahedron adjacent to _t2
    /// @param [out] returns the new tetrahedra created
    //----------------------------------------------------------------------------------------------------------------------
    Tetrahedron* flip32(Tetrahedron* _t1, Tetrahedron* _t2, Tetrahedron* _t3, int _pid, int _did);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function that performs flip14
    /// @param [in] _t1 the tetrahedra to be flipped
    /// @param [in] _t2 the tetrahedra to be flipped
    /// @param [in] _p the point to be inserted
    /// @param [in] _d the point in _t2 opp _p
    /// @param [in] _cid the position of the point shared by _t1 and _t2
    /// @param [out] returns the new tetrahedra created
    //----------------------------------------------------------------------------------------------------------------------
    Tetrahedron* flip44(Tetrahedron* _t1, Tetrahedron* _t2, ngl::Vec3 _p, ngl::Vec3 _d, int _cid);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function that calculates distance between two points A and B
    /// @param [in] _a point A
    /// @param [in] _b point B
    /// @param [out] returns the distance
    //----------------------------------------------------------------------------------------------------------------------
    double distance(ngl::Vec3 _a, ngl::Vec3 _b);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function that calculates the orientation of two points w.r.t plane A,B,C
    /// @param [in] _a point A
    /// @param [in] _b point B
    /// @param [in] _c point C
    /// @param [in] _p1 point to be checked for orientation
    /// @param [in] _p2 point to be checked for orientation
    /// @param [out] returns the valve depending on the orientation of _p1 and _p2
    //----------------------------------------------------------------------------------------------------------------------
    int oppSides(ngl::Vec3 _a, ngl::Vec3 _b, ngl::Vec3 _c, ngl::Vec3 _p1, ngl::Vec3 _p2);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function that creates the data for a flip
    /// @param [in] _pid position of point p in _t
    /// @param [in] _t the tetrahedron containing the point p
    /// @param [out] returns the stored flip structure
    //----------------------------------------------------------------------------------------------------------------------
    flipData createFlip(int _pid, Tetrahedron* _t);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief main function that checks for delaunay criterion
    /// @brief and performs the necessary flip based on that
    /// @param [out] returns the next tetrahedron to be checked
    //----------------------------------------------------------------------------------------------------------------------
    Tetrahedron* checkDelaunay();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function that checks if T and Ta are concave, convex or coplanar
    /// @param [in] _t1 tetrahedron to be checked
    /// @param [in] _t2 tetrahedron to be checked
    /// @param [out] returns the distance
    //----------------------------------------------------------------------------------------------------------------------
    int checkcase(Tetrahedron* _t1, Tetrahedron *_t2, ngl::Vec3 _point, ngl::Vec3 _d);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function that checks if T and Ta are concave, convex or coplanar
    /// @param [in] _oldt tetrahedron to be checked
    /// @param [in] _newt tetrahedron to be replaced with
    /// @param [in] _findt tetrahedron to be searched for
    //----------------------------------------------------------------------------------------------------------------------
    void updateNeighbour(Tetrahedron *_oldt, Tetrahedron *_newt, Tetrahedron *_findt);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function that returns the position of teh apex in t2
    /// @param [in] _t1 tetrahedron to be checked
    /// @param [in] _t2 tetrahedron to be checked
    /// @param [out] returns the position of the apex
    //----------------------------------------------------------------------------------------------------------------------
    int getApex(Tetrahedron *_t1, Tetrahedron *_t2);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function that checks for the position index of p in T
    /// @param [in] _t tetrahedron to be checked
    /// @param [in] _point the point to be searched
    /// @param [out] return the position index of point p in T
    //----------------------------------------------------------------------------------------------------------------------
    int findPoint(Tetrahedron *_t, ngl::Vec3 _point);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function that finds if a tetrahedron exist between t1 and t2
    /// @param [in] _t1 tetrahedron to be checked
    /// @param [in] _t2 tetrahedron to be checked
    /// @param [in] _c1 position of a vertex shared by t1 and t2 whose neighbour is checked
    /// @param [out] returns the third tetrahedron
    //----------------------------------------------------------------------------------------------------------------------
    Tetrahedron* findthird(Tetrahedron *_t1, Tetrahedron *_t2, int _c1);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function walks through the structure to find the tetrahedron that has the point p
    /// @param [in] _p point to be inserted
    /// @param [in] _t tetrahedron to be checked
    /// @param [in] returns the tetrahedron that contains p
    //----------------------------------------------------------------------------------------------------------------------
    Tetrahedron* walk(ngl::Vec3 _p, Tetrahedron *_t);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function that checks if vertices of t are coplanar
    /// @param [in] _t tetrahedron to be checked
    /// @param [out] returns if t is coplanar or not
    //----------------------------------------------------------------------------------------------------------------------
    bool isCoplanar(Tetrahedron *_t);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function that checks if t1 and t2 are in config with t3 and t4
    /// @param [in] _t1 tetrahedron to be checked
    /// @param [in] _t2 tetrahedron to be checked
    /// @param [in] _p point to be inserted
    /// @param [in] _cid position of a vertex shared by t1 and t2
    /// @param [in] _t3 tetrahedron to be checked
    /// @param [in] _t4 tetrahedron to be checked
    /// @param [out] returns if they are in config44 or not
    //----------------------------------------------------------------------------------------------------------------------
    bool config44(Tetrahedron *_t1, Tetrahedron *_t2, ngl::Vec3 _p, int _cid, Tetrahedron *&_t3, Tetrahedron *&_t4);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function that finds if a tetrahedron exist between t1 and t2
    /// @param [in] _t1 tetrahedron to be checked
    /// @param [in] _t2 tetrahedron to be checked
    /// @param [out] returns the third tetrahedron if exists
    //----------------------------------------------------------------------------------------------------------------------
    Tetrahedron* findthird(Tetrahedron *_t1, Tetrahedron *_t2);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief function that finds if a shared edge exists between t1,t2 and t3
    /// @param [in] _t1 tetrahedron to be checked
    /// @param [in] _t2 tetrahedron to be checked
    /// @param [in] _t3 tetrahedron to be checked
    //----------------------------------------------------------------------------------------------------------------------
    void findSharedEdge(Tetrahedron *_t1, Tetrahedron *_t2, Tetrahedron* _t3);
};

#endif // DELAUNAY_H
