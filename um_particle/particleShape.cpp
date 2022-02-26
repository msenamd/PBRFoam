
#include "particleShape.h"

particleShape::particleShape()
{
    currentSize = 0.0;
    air = NULL;
}

particleShape::~particleShape()
{
    //dtor
}

particleShape::particleShape(const particleShape& rhs)
{
    currentSize = rhs.currentSize;
    air = rhs.air;
}

particleShape& particleShape::operator=(const particleShape& rhs)
{
    if (&rhs != this){
        currentSize = rhs.currentSize;
        air = rhs.air;
    } ; // handle self assignment
    //assignment operator
    return *this;
}
