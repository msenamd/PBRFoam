#ifndef SLAB_H
#define SLAB_H

#include "particleShape.h"

class slab : public particleShape
{
public:
    slab();
    slab(double& halfThickness_, double& length_, double& width_);
    virtual ~slab();

    particleShape* clone() const;
    slab(const slab& rhs);
    slab& operator=(const slab& rhs);

    double get_convectiveHeatTransferCoefficient(double T_g , double T_surf , double u_g , double particleSize);
    double get_dragCoefficient(double T_g, double u_g, double particleSize);
    double get_volume();
    double get_surfaceArea();
    double get_surfaceArea(double size);
    double get_surfaceAreaToVolumeRatio();
    double get_projectedAreaRatio();
    double computeSizeFromVolume(const double volume_);
    void set_cellVolumes(int& numCells, std::vector<double>& xFacePositive, std::vector<double>& cellVolume);
    double correctForShape();

    void setDimensions(double& halfThickness_, double& length_, double& width_);

    double initialHalfThickness;
    double length;
    double width;
};

#endif // SLAB_H
