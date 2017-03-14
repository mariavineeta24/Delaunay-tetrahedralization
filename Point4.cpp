//----------------------------------------------------------------------------------------------------------------------
/// @file Point4.cpp
/// @brief This is the class that stores four points
/// (defines a tetrahedron)
//----------------------------------------------------------------------------------------------------------------------

#include "Point4.h"

std::ostream& operator<<(std::ostream& _output,const Point4 &_s)
{
    return _output<<_s.m_o<<"\n"<<
                    _s.m_a<<"\n"<<
                    _s.m_b<<"\n"<<
                    _s.m_c<<"\n";
}


