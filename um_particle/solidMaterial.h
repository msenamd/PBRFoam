#ifndef SOLIDMATERIAL_H
#define SOLIDMATERIAL_H

#include <string>
#include <iostream>
using namespace std;
#include <cmath>
#include <fstream>
#include <vector>
#include <tuple>
#include <algorithm>


class solidMaterial{

public:
    solidMaterial();
    ~solidMaterial();

    solidMaterial(const solidMaterial& rhs);
    solidMaterial& operator=(const solidMaterial& rhs);

    void set_knownType(const std::string& materialName);
    void set_wetSolid(const double& moistureFraction, const solidMaterial& drySolid);

    void set_bulkDensity(const double& bulkDensity_);
    void set_conductivity(const double& conductivity_, const double& conductivityExponent_);
    void set_radConductivity(const double& radConductivity_);
    void set_specificHeat(const double& specificHeat_, const double& specificHeatExponent_);
    void set_emissivity(const double& emissivity_);
    void set_porosity(const double& porosity_);
    void set_permeability(const double& permeability_);


    double get_bulkDensity(const double& temperature) const;
    double get_conductivity(const double& temperature) const;
    double get_radConductivity(const double& temperature) const;
    double get_specificHeat(const double& temperature) const;
    double get_emissivity(const double& temperature) const;
    double get_porosity(const double& temperature) const;
    double get_permeability(const double& temperature) const;


private:
    double  bulkDensity,                    // Density of porous material (value specific to porosity)
        conductivity, conductivityExponent, // Conductivity of solid, non-porous material (porous effects accounted for elsewhere)
        radConductivity,                    // Radiant conductivity of solid, non-porous material (porous effects accounted for elsewhere)
        specificHeat, specificHeatExponent, // Specific heat of solid, non-porous material (porous effects accounted for elsewhere)
        emissivity,                         // Emissivity of solid, non-porous material (porous effects accounted for elsewhere)
        porosity,
        permeability;
};

#endif //SOLIDMATERIAL_H
