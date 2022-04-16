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

// Set material properties from known type

void solidMaterial::set_knownType(const std::string& materialName)
{
    if(materialName == "white_pine")
    {
        bulkDensity = 361.0;
        conductivity = 0.176;
        conductivityExponent = 0.594;
        radConductivity = 0.0;
        specificHeat = 1664.0;
        specificHeatExponent = 0.660;
        emissivity = 0.759;
        porosity = 0.0744;
        permeability = 1e-10;

    }else if(materialName == "cardboard")
    {
        bulkDensity = 525.9;
        conductivity = 0.22;
        conductivityExponent = 0.0;
        radConductivity = 0.0;
        specificHeat = 2805.0;
        specificHeatExponent = 0.0;
        emissivity = 0.9;
        porosity = 0.0744;
        permeability = 1e-10;

    }else if(materialName == "Char")
    {
        bulkDensity = 73.0;
        conductivity = 0.065;
        conductivityExponent = 0.435;
        radConductivity = 3.3e-3;
        specificHeat = 1219.0;
        specificHeatExponent = 0.283;
        emissivity = 0.957;
        porosity = 0.8128;
        permeability = 1e-10;

    }else if(materialName == "ash")
    {
        bulkDensity = 5.7;
        conductivity = 0.058;
        conductivityExponent = 0.353;
        radConductivity = 6.4e-3;
        specificHeat = 1244.0;
        specificHeatExponent = 0.315;
        emissivity = 0.955;
        porosity = 0.9854;
        permeability = 1e-10;

    }else
    {
        throw exception();
    }
}

void solidMaterial::set_wetSolid(const double& moistureFraction, const solidMaterial& drySolid)
{
    this->operator =(drySolid);

    bulkDensity = drySolid.get_bulkDensity(300.0) * (1.0 + moistureFraction);
    conductivity = drySolid.get_conductivity(300.0) * (1.0 + 2.1*moistureFraction);
    specificHeat = drySolid.get_specificHeat(300.0) * (1.0 + 5.0*moistureFraction);
    porosity = 1.0 - bulkDensity / (drySolid.get_bulkDensity(300.0) / (1 - drySolid.get_porosity(300.0)));
}

// Set material properties from external input quantities

void solidMaterial::set_bulkDensity(const double& bulkDensity_)
{
    bulkDensity = bulkDensity_;
}

void solidMaterial::set_conductivity(const double& conductivity_, const double& conductivityExponent_)
{
    conductivity = conductivity_;
    conductivityExponent = conductivityExponent_;
}

void solidMaterial::set_radConductivity(const double& radConductivity_)
{
    radConductivity = radConductivity_;
}

void solidMaterial::set_specificHeat(const double& specificHeat_, const double& specificHeatExponent_)
{
    specificHeat = specificHeat_;
    specificHeatExponent = specificHeatExponent_;
}

void solidMaterial::set_emissivity(const double& emissivity_)
{
    emissivity = emissivity_;
}

void solidMaterial::set_porosity(const double& porosity_)
{
    porosity = porosity_;
}

void solidMaterial::set_permeability(const double& permeability_)
{
    permeability = permeability_;
}


// Return material properties//////////////////////////////////////

double solidMaterial::get_bulkDensity(const double& temperature) const
{
    return bulkDensity;
}

// Get bulk density of material by adding effects of the input moisture content (should only be used on materials initialized as dry solids)
double solidMaterial::get_bulkDensity(const double& temperature, const double& moistureFraction) const
{
    return bulkDensity * (1.0 + moistureFraction);
}

double solidMaterial::get_conductivity(const double& temperature) const
{
    return conductivity*pow(temperature/300.0, conductivityExponent);
}

// Get conductivity of material by adding effects of the input moisture content (should only be used on materials initialized as dry solids)
double solidMaterial::get_conductivity(const double& temperature, const double& moistureFraction) const
{
    return conductivity*pow(temperature/300.0, conductivityExponent) * (1.0 + 2.1*moistureFraction);
}

double solidMaterial::get_radConductivity(const double& temperature) const
{
    return radConductivity;
}

double solidMaterial::get_specificHeat(const double& temperature) const
{
    return specificHeat * pow(temperature/300.0, specificHeatExponent);
}

// Get specific heat of material by adding effects of the input moisture content (should only be used on materials initialized as dry solids)
double solidMaterial::get_specificHeat(const double& temperature, const double& moistureFraction) const
{
    return specificHeat * pow(temperature/300.0, specificHeatExponent) * (1.0 + 5.0*moistureFraction);
}

double solidMaterial::get_emissivity(const double& temperature) const
{
    return emissivity;
}

double solidMaterial::get_porosity(const double& temperature) const
{
    return porosity;
}

double solidMaterial::get_permeability(const double& temperature) const
{
    return permeability;
}
