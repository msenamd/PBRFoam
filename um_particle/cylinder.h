#ifndef CYLINDER_H
#define CYLINDER_H

#include "particleShape.h"

class cylinder : public particleShape
{
public:
    cylinder();
    cylinder(double& radius_, double& length_);
    virtual ~cylinder();

    particleShape* clone() const;
    cylinder(const cylinder& rhs);
    cylinder& operator=(const cylinder& rhs);

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

    void setDimensions(double& radius_, double& length_);

    double initialRadius;
    double length;
};

#endif // CYLINDER_H
