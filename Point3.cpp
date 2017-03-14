//----------------------------------------------------------------------------------------------------------------------
/// @file Point3.cpp
/// @brief This is the class that stores three points
/// (defines a triangle)
//----------------------------------------------------------------------------------------------------------------------

#include "Point3.h"

std::ostream& operator<<(std::ostream& _output,const Point3 &_s)
{
    return _output<<_s.m_a<<"\n"<<
                    _s.m_b<<"\n"<<
                    _s.m_c<<"\n";
}
