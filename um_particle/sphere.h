#ifndef SPHERE_H
#define SPHERE_H

#include "particleShape.h"

class sphere : public particleShape
{
public:
    sphere();
    sphere(double& radius_);
    virtual ~sphere();

    particleShape* clone() const;
    sphere(const sphere& rhs);
    sphere& operator=(const sphere& rhs);

    void initialize();
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

    void setDimensions(double& radius_);

    double initialRadius;
};

#endif // SPHERE_H
