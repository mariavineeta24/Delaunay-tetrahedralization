#ifndef POINT4_H
#define POINT4_H

//----------------------------------------------------------------------------------------------------------------------
/// @file Point4.h
/// @author Maria Vineeta Bagya Seelan
/// @version 1.0
/// @date 16/08/13
/// Revision History :
/// Initial Version 16/08/13
/// @class Point4
/// @brief This is the class that encapsulates a four ngl::Vec3 object
/// @brief (defines a tetrahedron with four vertices defined by OABC)
//----------------------------------------------------------------------------------------------------------------------

#include "ngl/Vec3.h"

class Point4
{
public:
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Constructor for Point4
    //----------------------------------------------------------------------------------------------------------------------
    Point4();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief constructor of Point4 that sets 4 ngl::Vec3 values(OABC)
    /// @param[in]  _o the O component
    /// @param[in]  _a the A component
    /// @param[in]  _b the B component
    /// @param[in]  _c the C component
    //----------------------------------------------------------------------------------------------------------------------
    Point4(
            ngl::Vec3 _o,
            ngl::Vec3 _a,
            ngl::Vec3 _b,
            ngl::Vec3 _c
          )
        : m_o(_o),
          m_a(_a),
          m_b(_b),
          m_c(_c){}
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Destructor for Point4
    //----------------------------------------------------------------------------------------------------------------------
    ~Point4(){}
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief accessor to get the point O
    //----------------------------------------------------------------------------------------------------------------------
    ngl::Vec3 getPointO() const { return m_o; }
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief accessor to get the point A
    //----------------------------------------------------------------------------------------------------------------------
    ngl::Vec3 getPointA() const { return m_a; }
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief accessor to get the point B
    //----------------------------------------------------------------------------------------------------------------------
    ngl::Vec3 getPointB() const { return m_b; }
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief accessor to get the point C
    //----------------------------------------------------------------------------------------------------------------------
    ngl::Vec3 getPointC() const { return m_c; }
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief sets the O component
    //----------------------------------------------------------------------------------------------------------------------
    void setPointO(ngl::Vec3 _o) { m_a = _o; }
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief sets the A component
    //----------------------------------------------------------------------------------------------------------------------
    void setPointA(ngl::Vec3 _a) { m_a = _a; }
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief sets the B component
    //----------------------------------------------------------------------------------------------------------------------
    void setPointB(ngl::Vec3 _b) { m_b = _b; }
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief sets the C component
    //----------------------------------------------------------------------------------------------------------------------
    void setPointC(ngl::Vec3 _c) { m_c = _c; }
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief insertion operator to print out the Point4
    /// @param[in] _output the stream to write to
    /// @param[in] _s the Point4 to write
    //----------------------------------------------------------------------------------------------------------------------
    friend std::ostream& operator<<(std::ostream& _output, const Point4 &_s);

public :
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the O component
    //----------------------------------------------------------------------------------------------------------------------
    ngl::Vec3 m_o;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the A component
    //----------------------------------------------------------------------------------------------------------------------
    ngl::Vec3 m_a;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the B component
    //----------------------------------------------------------------------------------------------------------------------
    ngl::Vec3 m_b;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the C component
    //----------------------------------------------------------------------------------------------------------------------
    ngl::Vec3 m_c;
    //----------------------------------------------------------------------------------------------------------------------

private:

};

#endif // POINT4_H
