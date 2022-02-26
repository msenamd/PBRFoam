#ifndef PARTICLESHAPE_H
#define PARTICLESHAPE_H

#include <vector>
#include <math.h>

#include "air.h"

class particleShape
{
public:
    particleShape();
    virtual ~particleShape();

    virtual particleShape* clone() const = 0;
    particleShape(const particleShape& rhs);
    particleShape& operator=(const particleShape& rhs);

    virtual void initialize() = 0;
    virtual double get_convectiveHeatTransferCoefficient(double T_g , double T_surf , double u_g , double delta) = 0;
    virtual double get_dragCoefficient(double T_g, double u_g, double delta) = 0;
    virtual double get_volume() = 0;
    virtual double get_surfaceArea() = 0;
    virtual double get_surfaceArea(double size) = 0;
    virtual double get_surfaceAreaToVolumeRatio() = 0;
    virtual double get_projectedAreaRatio() = 0;
    virtual double computeSizeFromVolume(const double volume_) = 0;
    virtual void set_cellVolumes(int& numCells, std::vector<double>& xFacePositive, std::vector<double>& cellVolume) = 0;
    virtual double correctForShape() = 0;

    double currentSize;         //instantaneous particle size (radius (cylinder, sphere) or half-thickness (slab) in meters)
    Air* air;                   //For air/gas properties

protected:
    static constexpr double pi = 3.14159265358979323846;
};

#endif // PARTICLESHAPE_H
