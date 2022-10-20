#ifndef MATERIAL_H
#define MATERIAL_H

#include <string>
#include <iostream>
using namespace std;
#include <fstream>
#include <vector>
#include <tuple>
#include <algorithm>


class Material{

public:
    Material();
    ~Material();
    Material(const Material &rhs);
    Material& operator=(const Material &rhs);

    enum eProperty
    {
        T,
        density,
        specificHeat,
        conductivity,
        emissivity,
        heatOfCombustion
    };

    void set_constantProperties(double density_, double specificHeat_, double conductivity_, double emissivity_, double heatOfCombustion_);
    void set_constant_density(double density_);
    void set_constant_specificHeat(double specificHeat_);
    void set_constant_conductivity(double conductivity_);
    void set_constant_emissivity(double emmisivity_);
    void set_constant_heatOfCombustion(double heatOfCombustion_);
    void set_combustionTemp(double temperature_);
    bool set_knownMaterial(std::string materialName);
    bool set_customName(std::string materialName);
    void push_temperatureRow(double temperature_, double density_, double specificHeat_, double conductivity_, double emissivity_, double heatOfCombustion_);

    double get_density(double temperature);
    double get_specificHeat(double temperature);
    double get_conductivity(double temperature);
    double get_emissivity(double temperature);
    double get_heatOfCombustion(double temperature);
    string& get_name();

    void clearData();

private:
    void memberwiseCopyAssignment(const Material& rhs);
    void checkTemperatureBounds(double temperature);
    double get_property(eProperty property, double temperature);

    string name;
    std::vector<tuple<double, double, double, double, double, double> > propertyTable;
};

#endif //MATERIAL_H
