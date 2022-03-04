#include "sphere.h"

sphere::sphere() : particleShape()
{
    initialRadius = 0.0;
    currentSize = initialRadius;
}

sphere::sphere(double& radius_) : particleShape()
{
    initialRadius = radius_;
    currentSize = initialRadius;
}

sphere::~sphere()
{

}

particleShape* sphere::clone() const
{
    return new sphere(*this);
}

sphere::sphere(const sphere& rhs) : particleShape(rhs)
{
    initialRadius = rhs.initialRadius;
}

sphere& sphere::operator=(const sphere& rhs)
{
    if (&rhs != this){
        particleShape::operator=(rhs);
        initialRadius = rhs.initialRadius;

    } ; // handle self assignment
    //assignment operator
    return *this;
}

double sphere::get_volume()
{
    return 4.0 / 3.0 * pi * pow(currentSize, 3.0);
}

/**
* Computes particle surface area using the stored value of currentSize (radius or half-thickness).
* @return Single-sided surface area of particle (m^2).
*/
double sphere::get_surfaceArea()
{
    return 4.0 * pi * pow(currentSize,2.0);
}

/**
* Computes particle surface area (single sided) given a size (radius or thickness).  Does not use stored currentSize of particle.
* @param size Radius or thickness of particle (m).
* @return Single-sided surface area of particle (m^2).
*/
double sphere::get_surfaceArea(double size)
{
    return 4.0 * pi * pow(size,2.0);
}

/**
* Computes particle surface area to volume ratio using the stored value of currentSize (radius).
* @return Surface area to volume ratio of particle (m^-1).
*/
double sphere::get_surfaceAreaToVolumeRatio()
{
    return 3.0 / currentSize;
}

/**
* Gets a particle radius given a particle volume.  Does not alter values stored in particle.
* @param volume_ Volume of sphere (m^3).
* @return Radius of spherical particle (m).
*/
double sphere::computeSizeFromVolume(const double volume_)
{
    return pow((3.0 * volume_ / 4.0 / pi), (1.0 / 3.0));
}

/**
* Computes and sets the cellVolume array from given xFacePostive array.
* @param numCells Number of cells in the arrays.
* @param xFacePositive Array of radial coordinates of cell positive (outward direction) faces (m).
* @param cellVolume Cell volume array filled by this function (m^3).
* @return
*/
void sphere::set_cellVolumes(int& numCells, std::vector<double>& xFacePositive, std::vector<double>& cellVolume)
{
    cellVolume[0] = 4.0/3.0 * pi * pow(xFacePositive[0] , 3.0);
    for (int i = 1; i < numCells; i++)
    {
        cellVolume[i] = 4.0/3.0 * pi * (pow(xFacePositive[i] , 3.0) - pow(xFacePositive[i-1] , 3.0));
    }

//    cellVolume[0] = 4.0/3.0 * pi * sizeCell[0] * pow(xFacePositive[0] , 2.0);
//    for (int i = 1; i < numCells; i++)
//    {
//        cellVolume[i] = 4.0/3.0 * pi * sizeCell[i] *
//            ( pow(xFacePositive[i] , 2.0) + xFacePositive[i] * xFacePositive[i-1] + pow(xFacePositive[i-1] , 2.0));
//    }
}

/**
* Computes a correction to the computed quantities based on shape.  Rectangular (slab) shapes 
* are doubled to account for the two outer surfaces, other shapes (cylinder, sphere) are not doubled.
* @param
* @return a correction factor
*/
double sphere::correctForShape()
{
    return 1.0;
}

/**
* Initialize shape object.
* @return
*/
void sphere::initialize()
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
double sphere::get_convectiveHeatTransferCoefficient(double T_g , double T_surf , double u_g , double particleSize)
{
    //    This function calculates heat transfer coefficient from an empirical
    //     * correlation of a sphere
    //     * Inputs: gas temperature, surface temp, gas velocity, particle size
    //     *
    //     *
    //     * sphere: Ref: Whitaker


    double Pr, nu_g, k_g;
    double D_eff, Re_D, Nu_D , h;
    double mu_g, mu_s;

    Pr = air->get_pr(T_g);
    nu_g = air->get_v(T_g);
    k_g = air->get_k(T_g);
    mu_g = air->get_mu(T_g);
    mu_s = air->get_mu(T_surf);

    D_eff       = 2*particleSize;                             // Effective diameter of particle [m]
    Re_D        = u_g*D_eff/nu_g;                      // Calculated at gas temperature
    Nu_D        = 2+(0.4*pow(Re_D,0.5)+0.06*pow(Re_D, 2.0/3.0))*pow(Pr,0.4)*pow(mu_g/mu_s,0.25);

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
double sphere::get_dragCoefficient(double T_g, double u_g, double particleSize)
{
    //    This function calculates the drag coefficient from an empirical
    //     * correlation of a sphere
    //     * Inputs: gas temperature, gas velocity, particle size
    //     *
    //     * sphere Ref: R. Panton, Incompressible Flow


    double Re, CD;

    Re = u_g * particleSize*2 / air->get_v(T_g);

    if (Re <= 0.001)
    {
        CD = 100;
    }
    else if (Re <= 1.0)
    {
        CD = 24.0 / Re;
    }
    else if (Re > 1.0 && Re <= 1000)
    {
        CD = 24.0 / Re * (0.85 + 0.15 * pow(Re, 0.687));
    }
    else
    {
        CD = 0.44;
    }

    return CD;
}

/**
* Computes particle projected area to surface area ratio.
* @return projected to surface area ratio of particle (-).
*/
double sphere::get_projectedAreaRatio()
{
    return 1.0 / 4.0;
}

void sphere::setDimensions(double& radius_)
{
    initialRadius = radius_;
    currentSize = initialRadius;
}
