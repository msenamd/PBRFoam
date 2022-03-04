#include "cylinder.h"

cylinder::cylinder() : particleShape()
{
    initialRadius = 0.0;
    length = 0.0;
    currentSize = initialRadius;
}

cylinder::cylinder(double& radius_, double& length_) : particleShape()
{
    initialRadius = radius_;
    length = length_;
    currentSize = initialRadius;
}

cylinder::~cylinder()
{

}

particleShape* cylinder::clone() const
{
    return new cylinder(*this);
}

cylinder::cylinder(const cylinder& rhs) : particleShape(rhs)
{
    initialRadius = rhs.initialRadius;
    length = rhs.length;
}

cylinder& cylinder::operator=(const cylinder& rhs)
{
    if (&rhs != this){
        particleShape::operator=(rhs);
        initialRadius = rhs.initialRadius;
        length = rhs.length;

    } ; // handle self assignment
    //assignment operator
    return *this;
}

double cylinder::get_volume()
{
    return pi * pow(currentSize, 2.0) * length;
}

/**
* Computes particle surface area using the stored value of currentSize (radius or half-thickness).
* @return Single-sided surface area of particle (m^2).
*/
double cylinder::get_surfaceArea()
{
    return 2.0 * pi * currentSize * length;
}

/**
* Computes particle surface area (single sided) given a size (radius or thickness).  Does not use stored currentSize of particle.
* @param size Radius or thickness of particle (m).
* @return Single-sided surface area of particle (m^2).
*/
double cylinder::get_surfaceArea(double size)
{
    return 2.0 * pi * size * length;
}

/**
* Computes particle surface area to volume ratio using the stored value of currentSize (radius).
* @return Surface area to volume ratio of particle (m^-1).
*/
double cylinder::get_surfaceAreaToVolumeRatio()
{
    return 2.0 / currentSize;
}

/**
* Gets a particle radius given a particle volume.  Does not alter values stored in particle.
* @param volume_ Volume of cylinder (m^3).
* @return Radius of cylinder (m).
*/
double cylinder::computeSizeFromVolume(const double volume_)
{
    return pow((volume_ / length / pi) , 0.5);
}

/**
* Computes and sets the cellVolume array from given xFacePostive array.
* @param numCells Number of cells in the arrays.
* @param xFacePositive Array of radial coordinates of cell positive (outward direction) faces (m).
* @param cellVolume Cell volume array filled by this function (m^3).
* @return
*/
void cylinder::set_cellVolumes(int& numCells, std::vector<double>& xFacePositive, std::vector<double>& cellVolume)
{
    cellVolume[0] = pi * pow(xFacePositive[0] , 2.0) * length;
    for (int i = 1; i < numCells; i++)
    {
        cellVolume[i] = pi * (pow(xFacePositive[i] , 2.0) - pow(xFacePositive[i-1] , 2.0)) * length;
    }

//    cellVolume[0] = pi * sizeCell[0] * (xFacePositive[0]) * lengthCylinder;
//    for (int i=1; i < numCells; i++)
//    {
//        cellVolume[i] = pi * sizeCell[i] * (xFacePositive[i] + xFacePositive[i-1]) * lengthCylinder;
//    }
}

/**
* Computes a correction to the computed quantities based on shape.  Rectangular (slab) shapes 
* are doubled to account for the two outer surfaces, other shapes (cylinder, sphere) are not doubled.
* @param
* @return a correcton factor
*/
double cylinder::correctForShape()
{
    return 1.0;
}

/**
* Initialize shape object.
* @return
*/
void cylinder::initialize()
{
    currentSize = initialRadius;
}

/**
* Computes particle convective heat transfer coefficient.
* @param T_g Gas temperature (K).
* @param T_surf Particle surface temperature (K).
* @param u_g Gas velocity (m/s).
* @param particleSize Radius or half-thickness of particle (m).
* @return
*/
double cylinder::get_convectiveHeatTransferCoefficient(double T_g , double T_surf , double u_g , double particleSize)
{
    //    This function calculates the heat transfer coefficient from an empirical
    //     * correlation of a cylinder
    //     * Inputs: gas temperature, surface temp, gas velocity, particle size, particle geometry
    //     *
    //     * cylinder Ref: Churchill and Bernstein


    double Pr, T_film , nu_g, k_g;
    double D_eff, Re_D, Nu_D , h;

    T_film      = 0.5*(T_surf+T_g);         // Film temperature [K]

    Pr = air->get_pr(T_film);
    nu_g = air->get_v(T_film);
    k_g = air->get_k(T_film);

    D_eff       = 2*particleSize;              // Effective diameter of particle [m]
    Re_D        = u_g*D_eff/nu_g;       // Reynolds number at film temp.
    // Nusselt number at film temp. , Ref: Churchill and Bernstein
    Nu_D        = 0.3 + (0.62*pow(Re_D,0.5)*pow(Pr,0.3333))/pow((1+pow(0.4/Pr,2.0/3.0)),0.25)
            * pow((1+pow(Re_D/282000,5.0/8.0)),4.0/5.0) ;

    h    = Nu_D*k_g/D_eff;   // Convective heat transfer coef. [W/m2/K]

    return h;
}

/**
* Computes particle drag coefficient.
* @param T_g Gas temperature (K).
* @param u_g Gas velocity (m/s).
* @param particleSize Radius or half-thickness of particle (m).
* @return
*/
double cylinder::get_dragCoefficient(double T_g, double u_g, double particleSize)
{
    //    This function calculates the drag coefficient from an empirical
    //     * correlation of a cylinder
    //     * Inputs: gas temperature, gas velocity, particle size
    //     *
    //     * cylinder Ref: R. Panton, Incompressible Flow


    double Re, CD;

    Re = u_g * particleSize*2 / air->get_v(T_g);

    if (Re <= 0.001)
    {
        CD = 1000;
    }
    else if (Re <= 1.0)
    {
        CD = 10.0 / pow(Re, 0.8);
    }
    else if (Re > 1.0 && Re <= 1000)
    {
        CD = 10 * (0.6 + 0.4 * pow(Re, 0.8)) / Re;
    }
    else
    {
        CD = 1.0;
    }

    return CD;
}

/**
* Computes particle projected area to surface area ratio.
* @return projected to surface area ratio of particle (-).
*/
double cylinder::get_projectedAreaRatio()
{
    return 1.0/pi;
}

void cylinder::setDimensions(double& radius_, double& length_)
{
    initialRadius = radius_;
    length = length_;
    currentSize = initialRadius;
}
