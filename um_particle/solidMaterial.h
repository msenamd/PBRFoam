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

    void set_bulkDensity(double bulkDensity_);
    void set_conductivity(double conductivity_, double conductivityExponent_);
    void set_radConductivity(double radConductivity_);
    void set_specificHeat(double specificHeat_, double specificHeatExponent_);
    void set_emissivity(double emissivity_);
    void set_porosity(double porosity_);
    void set_permeability(double permeability_);


    double get_bulkDensity(double temperature);
    double get_conductivity(double temperature);
    double get_radConductivity(double temperature);
    double get_specificHeat(double temperature);
    double get_emissivity(double temperature);
    double get_porosity(double temperature);
    double get_permeability(double temperature);


private:
    double  bulkDensity,
        conductivity, conductivityExponent,
        radConductivity,
        specificHeat, specificHeatExponent,
        emissivity,
        porosity,
        permeability;
};

#endif //SOLIDMATERIAL_H
