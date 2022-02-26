#include "solidMaterial.h"

solidMaterial::solidMaterial()
{
    bulkDensity = 0.0;

    conductivity = 0.0;
    conductivityExponent = 0.0;

    radConductivity = 0.0;

    specificHeat = 0.0;
    specificHeatExponent = 0.0;

    emissivity = 0.0;

    porosity = 0.0;

    permeability = 0.0;
}

solidMaterial::~solidMaterial()
{

}

solidMaterial::solidMaterial(const solidMaterial& rhs)
{
    bulkDensity = rhs.bulkDensity;

    conductivity = rhs.conductivity;
    conductivityExponent = rhs.conductivityExponent;

    radConductivity = rhs.radConductivity;

    specificHeat = rhs.specificHeat;
    specificHeatExponent = rhs.specificHeatExponent;

    emissivity = rhs.emissivity;

    porosity = rhs.porosity;

    permeability = rhs.permeability;
}

solidMaterial& solidMaterial::operator=(const solidMaterial& rhs)
{
    if (&rhs != this) {
        bulkDensity = rhs.bulkDensity;

        conductivity = rhs.conductivity;
        conductivityExponent = rhs.conductivityExponent;

        radConductivity = rhs.radConductivity;

        specificHeat = rhs.specificHeat;
        specificHeatExponent = rhs.specificHeatExponent;

        emissivity = rhs.emissivity;

        porosity = rhs.porosity;

        permeability = rhs.permeability;
    }; // handle self assignment
    //assignment operator
    return *this;
}

// Set material properties from external input qnatities

void solidMaterial::set_bulkDensity(double bulkDensity_)
{
    bulkDensity = bulkDensity_;
}

void solidMaterial::set_conductivity(double conductivity_, double conductivityExponent_)
{
    conductivity = conductivity_;
    conductivityExponent = conductivityExponent_;
}

void solidMaterial::set_radConductivity(double radConductivity_)
{
    radConductivity = radConductivity_;
}

void solidMaterial::set_specificHeat(double specificHeat_, double specificHeatExponent_)
{
    specificHeat = specificHeat_;
    specificHeatExponent = specificHeatExponent_;
}

void solidMaterial::set_emissivity(double emissivity_)
{
    emissivity = emissivity_;
}

void solidMaterial::set_porosity(double porosity_)
{
    porosity = porosity_;
}

void solidMaterial::set_permeability(double permeability_)
{
    permeability = permeability_;
}


// Return material properties

double solidMaterial::get_bulkDensity(double temperature)
{
    return bulkDensity;
}

double solidMaterial::get_conductivity(double temperature)
{
    return conductivity*pow(temperature/300.0, conductivityExponent);
}

double solidMaterial::get_radConductivity(double temperature)
{
    return radConductivity;
}

double solidMaterial::get_specificHeat(double temperature)
{
    return specificHeat * pow(temperature/300.0, specificHeatExponent);
}

double solidMaterial::get_emissivity(double temperature)
{
    return emissivity;
}

double solidMaterial::get_porosity(double temperature)
{
    return porosity;
}

double solidMaterial::get_permeability(double temperature)
{
    return permeability;
}