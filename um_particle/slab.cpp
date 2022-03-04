#include "slab.h"

slab::slab() : particleShape()
{
    initialHalfThickness = 0.0;
    length = 0.0;
    width = 0.0;
    currentSize = initialHalfThickness;
}

slab::slab(double& halfThickness_, double& length_, double& width_) : particleShape()
{
    initialHalfThickness = halfThickness_;
    length = length_;
    width = width_;
    currentSize = initialHalfThickness;
}

slab::~slab()
{

}

particleShape* slab::clone() const
{
    return new slab(*this);
}

slab::slab(const slab& rhs) : particleShape(rhs)
{
    initialHalfThickness = rhs.initialHalfThickness;
    length = rhs.length;
    width = rhs.width;
}

slab& slab::operator=(const slab& rhs)
{
    if (&rhs != this){
        particleShape::operator=(rhs);
        initialHalfThickness = rhs.initialHalfThickness;
        length = rhs.length;
        width = rhs.width;

    } ; // handle self assignment
    //assignment operator
    return *this;
}

double slab::get_volume()
{
    return 2.0 * currentSize * length * width;
}

/**
* Computes particle surface area using the stored value of currentSize (radius or half-thickness).
* @return Single-sided surface area of particle (m^2).
*/
double slab::get_surfaceArea()
{
    return length * width;
}

/**
* Computes particle surface area (single sided) given a size (radius or thickness).  Does not use stored currentSize of particle.
* @param size Radius or thickness of particle (m).
* @return Single-sided surface area of particle (m^2).
*/
double slab::get_surfaceArea(double size)
{
    return length * width;
}

/**
* Computes particle surface area to volume ratio using the stored value of currentSize (half-thickness).
* @return Surface area to volume ratio of particle (m^-1).
*/
double slab::get_surfaceAreaToVolumeRatio()
{
    return 1.0 / (2.0*currentSize);
}

/**
* Gets a particle half-thickness given a particle volume.  Does not alter values stored in particle.
* @param volume_ Volume of slab (m^3).
* @return Half-thickness of slab (m).
*/
double slab::computeSizeFromVolume(const double volume_)
{
    return volume_ / (2.0*length*width);
}

/**
* Computes and sets the cellVolume array from given xFacePostive array.
* @param numCells Number of cells in the arrays.
* @param xFacePositive Array of radial coordinates of cell positive (outward direction) faces (m).
* @param cellVolume Cell volume array filled by this function (m^3).
* @return
*/
void slab::set_cellVolumes(int& numCells, std::vector<double>& xFacePositive, std::vector<double>& cellVolume)
{
    cellVolume[0] = xFacePositive[0] * length * width;
    for (int i = 1; i < numCells; i++)
    {
        cellVolume[i] = (xFacePositive[i] - xFacePositive[i-1]) * length * width;
    }

//    for (int i = 0; i < numCells; i++)
//    {
//        cellVolume[i] = sizeCell[i] * areaRectangle;
//    }
}

/**
* Computes a correction to the computed quantities based on shape.  Rectangular (slab) shapes 
* are doubled to account for the two outer surfaces, other shapes (cylinder, sphere) are not doubled.
* @param
* @return a correction factor
*/
double slab::correctForShape()
{
    return 2.0;
}

/**
* Initialize shape object.
* @return
*/
void slab::initialize()
{
    currentSize = initialHalfThickness;
}

/**
* Computes particle convective heat transfer coefficient.
* @param T_g Gas temperature (K).
* @param T_surf Particle surface temperature (K).
* @param u_g Gas velocity (m/s).
* @param particleSize Radius or half-thickness of particle (m).
* @return
*/
double slab::get_convectiveHeatTransferCoefficient(double T_g , double T_surf , double u_g , double particleSize)
{
    //    This function calculates the heat transfer coefficient from an empirical
    //     * correlation of a rectangle
    //     * Inputs: gas temperature, surface temp, gas velocity, particle size
    //     *
    //     * rectangle Ref: Churchill and Bernstein


    double Pr, T_film , nu_g, k_g;
    double D_eff, Re_D, Nu_D , h;

    T_film      = 0.5*(T_surf+T_g);         // Film temperature [K]

    Pr = air->get_pr(T_film);
    nu_g = air->get_v(T_film);
    k_g = air->get_k(T_film);

    D_eff       = (2*particleSize)*2/sqrt(pi); // Effective diameter of particle [m]f, calculated using equal area method
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
double slab::get_dragCoefficient(double T_g, double u_g, double particleSize)
{
    //    This function calculates the drag coefficient from an empirical
    //     * correlation of a rectangle
    //     * Inputs: gas temperature, gas velocity, particle size
    //     *
    //     * rectangle Ref: Sardoy et al., Comb. Flame, 2007


    double AR, CD;

    AR = length / (2 * particleSize);

    CD = 1.98 - 0.8 * (1 - exp(-20.0 / AR));

    return CD;
}

/**
* Computes particle projected area to surface area ratio.
* @return projected to surface area ratio of particle (-).
*/
double slab::get_projectedAreaRatio()
{
    return 1.0;
}


void slab::setDimensions(double& halfThickness_, double& length_, double& width_)
{
    initialHalfThickness = halfThickness_;
    length = length_;
    width = width_;
    currentSize = initialHalfThickness;
}
