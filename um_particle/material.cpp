#include "material.h"

Material::Material()
{

}

Material::~Material()
{

}

/**
* Copy constructor.
* @param rhs material to copy.
* @return
*/
Material::Material(const Material &rhs)
{
    memberwiseCopyAssignment(rhs);
}

/**
* Equality overload.
* @param rhs material to copy.
* @return
*/
Material& Material::operator=(const Material &rhs)
{
    if (this != &rhs)
    {
        memberwiseCopyAssignment(rhs);
    }
    return *this;
}

/**
* Memberwise copy assignment
* @param Makes a deep copy of all member variables
* @return
*/
void Material::memberwiseCopyAssignment(const Material &rhs)
{
    name = rhs.name;
    propertyTable = rhs.propertyTable;
}

/**
* Set material properties based on a built-in type.
* @param materialName Name of built-in material type (wood, woodRothermel, cardboard).
* @return True if set properly, false if not set (because name does not match known names).
*/
bool Material::set_knownMaterial(std::string materialName)
{
    if(materialName == "testSolid")
    {
        name = materialName;
        set_constantProperties(663.0, 2520.0, 0.126, 0.9, 0.0); // Note that heat of combustion set to 0 here because this is used with an explicit combustion reaction rate formulation currently (so heat of combustion not used)
        return true;
    }
    else if(materialName == "moisture")
    {
        name = materialName;
        set_constantProperties(1000.0, 4180.0, 0.6, 0.9, 0.0);  // Note that heat of combustion set to 0.
        return true;
    }
    else if(materialName == "char")
    {
        name = materialName;
        set_constantProperties(132.6, 2520.0, 0.126, 0.9, 0.0); // Note that heat of combustion set to 0 here because this is used with an explicit combustion reaction rate formulation currently (so heat of combustion not used)
        return true;
    }
    else if(materialName == "wood")
    {
        name = materialName;
        set_constantProperties(640.0, 2805.0, 0.15, 0.9, 13000000.0);
        return true;
    }else if(materialName == "woodRothermel")
    {
        name = materialName;
        set_constantProperties(398.0, 2805.0, 0.15, 0.9, 20090000.0);
        return true;
    }else if(materialName == "cardboard")
    {
        name = materialName;
        set_constantProperties(525.9, 2805.0, 0.22, 0.9, 13600000.0); //Babrauskas => hc=13600000.0 ,
                                                                        //Agarwal et al. 2014, Fire Safety Science-Proceedings of the Eleventh International Symposium, 124-137 => hc=16200000.0
        return true;
    }else
    {
        return false;   //unknown name
    }
}

/**
* Set custom material name.  Does not set the properties.
* @param materialName Name of material type (cannot be wood, woodRothermel, cardboard).
* @return True if set properly, false if not set (because name matches a built-in name).
*/
bool Material::set_customName(std::string materialName)
{
    if(materialName == "wood")
    {
        return false;
    }else if(materialName == "woodRothermel")
    {
        return false;
    }else if(materialName == "cardboard")
    {
        return false;
    }else
    {
        name = materialName;
        return true;   //unknown name
    }
}

/**
* Gets material property using linear interpolation in temperature.
* @param property Material property to get (rho, cSubP, k, epsilon).
* @param temperature Temperature to get property at.
* @return Value of the property.
*/
double Material::get_property(eProperty property, double temperature)
{

    int test = 0;
    if(temperature > std::get<T>(propertyTable.back()))
    {
        cout << "Invalid temperature: " << temperature << "\n";
        return -1;
    }
    if(temperature < std::get<T>(propertyTable.front()))
    {
        cout << "Invalid temperature: " << temperature << "\n";
        return -1;
    }
    else
    {
        int i = 0;
        while(temperature > std::get<T>(propertyTable[i]))
            i++;
        if(property == density)
        {
            return (((temperature - std::get<T>(propertyTable[i-1]))/(std::get<T>(propertyTable[i]) - std::get<T>(propertyTable[i-1])))*(std::get<density>(propertyTable[i]) - std::get<density>(propertyTable[i-1]))) + std::get<density>(propertyTable[i-1]);
        }else if(property == specificHeat)
        {
            return (((temperature - std::get<T>(propertyTable[i-1]))/(std::get<T>(propertyTable[i]) - std::get<T>(propertyTable[i-1])))*(std::get<specificHeat>(propertyTable[i]) - std::get<specificHeat>(propertyTable[i-1]))) + std::get<specificHeat>(propertyTable[i-1]);
        }else if(property == conductivity)
        {
            return (((temperature - std::get<T>(propertyTable[i-1]))/(std::get<T>(propertyTable[i]) - std::get<T>(propertyTable[i-1])))*(std::get<conductivity>(propertyTable[i]) - std::get<conductivity>(propertyTable[i-1]))) + std::get<conductivity>(propertyTable[i-1]);
        }else if(property == emissivity)
        {
            return (((temperature - std::get<T>(propertyTable[i-1]))/(std::get<T>(propertyTable[i]) - std::get<T>(propertyTable[i-1])))*(std::get<emissivity>(propertyTable[i]) - std::get<emissivity>(propertyTable[i-1]))) + std::get<emissivity>(propertyTable[i-1]);
        }else if(property == heatOfCombustion)
        {
            return (((temperature - std::get<T>(propertyTable[i-1]))/(std::get<T>(propertyTable[i]) - std::get<T>(propertyTable[i-1])))*(std::get<heatOfCombustion>(propertyTable[i]) - std::get<heatOfCombustion>(propertyTable[i-1]))) + std::get<heatOfCombustion>(propertyTable[i-1]);
        }else
        {
            cout << "Invalid property switch: " << property << "\n";
            return -1;
        }
    }
}

/**
* Gets density using linear interpolation in temperature.
* @param temperature Temperature to get density at.
* @return Value of density (kg/m^3).
*/
double Material::get_density(double temperature)
{
    return get_property(eProperty::density, temperature);
}

/**
* Gets specific heat using linear interpolation in temperature.
* @param temperature Temperature to get specific heat at.
* @return Value of specific heat (J/kg*K).
*/
double Material::get_specificHeat(double temperature)
{
    return get_property(eProperty::specificHeat, temperature);
}

/**
* Gets conductivity using linear interpolation in temperature.
* @param temperature Temperature to get conductivity at.
* @return Value of conductivity (w/m*K).
*/
double Material::get_conductivity(double temperature)
{
    return get_property(eProperty::conductivity, temperature);
}

/**
* Gets emmisivity using linear interpolation in temperature.
* @param temperature Temperature to get emmisivity at.
* @return Value of emmisivity.
*/
double Material::get_emissivity(double temperature)
{
    return get_property(eProperty::emissivity, temperature);
}

/**
* Gets heat of combustion using linear interpolation in temperature.
* @param temperature Temperature to get heat of combustion at.
* @return Value of heat of combustion (J/kg).
*/
double Material::get_heatOfCombustion(double temperature)
{
    return get_property(eProperty::heatOfCombustion, temperature);
}

/**
* Gets name of material.
* @return Value of material name.
*/
string& Material::get_name()
{
    return name;
}

/**
* Sets properties to be constant with temperature.
* @param density_ Density at temperature (kg/m^3).
* @param specificHeat_ Specific heat at temperature (J/kg*K).
* @param conductivity_ Conductivity at temperature (w/m*K).
* @param emissivity_ Emmisivity at temperature.
* @param heatOfCombustion_ Heat of combustion at temperature (J/kg).
* @return
*/
void Material::set_constantProperties(double density_, double specificHeat_, double conductivity_, double emissivity_, double heatOfCombustion_)
{
    propertyTable.clear();
    propertyTable.push_back(std::make_tuple(0.0, density_, specificHeat_, conductivity_, emissivity_, heatOfCombustion_));
    propertyTable.push_back(std::make_tuple(10000.0, density_, specificHeat_, conductivity_, emissivity_, heatOfCombustion_));
    sort(propertyTable.begin(), propertyTable.end());
}

/**
* Sets density to be a constant value with temperature.
* @param density_ Density (kg/m^3).
* @return
*/
void Material::set_constant_density(double density_)
{
    for(int i=0; i<propertyTable.size(); i++)
    {
        std::get<density>(propertyTable[i]) = density_;
    }
}

/**
* Sets specific heat to be a constant value with temperature.
* @param specificHeat_ Specific heat (J/kg*K).
* @return
*/
void Material::set_constant_specificHeat(double specificHeat_)
{
    for(int i=0; i<propertyTable.size(); i++)
    {
        std::get<specificHeat>(propertyTable[i]) = specificHeat_;
    }
}

/**
* Sets conductivity to be a constant value with temperature.
* @param conductivity_ Value of conductivity (w/m*K).
* @return
*/
void Material::set_constant_conductivity(double conductivity_)
{
    for(int i=0; i<propertyTable.size(); i++)
    {
        std::get<conductivity>(propertyTable[i]) = conductivity_;
    }
}

/**
* Sets Emmisivity to be a constant value with temperature.
* @param emmisivity_ Value of emmisivity.
* @return
*/
void Material::set_constant_emissivity(double emmisivity_)
{
    for(int i=0; i<propertyTable.size(); i++)
    {
        std::get<emissivity>(propertyTable[i]) = emmisivity_;
    }
}

/**
* Sets heat of combustion to be a constant value with temperature.
* @param heatOfCombustion_ Value of heat of combustion (J/kg).
* @return
*/
void Material::set_constant_heatOfCombustion(double heatOfCombustion_)
{
    for(int i=0; i<propertyTable.size(); i++)
    {
        std::get<heatOfCombustion>(propertyTable[i]) = heatOfCombustion_;
    }
}

/**
* Push on row of properties onto the property table.
* @param temperature_ Temperature corresponding to the properties (K).
* @param density_ Density at temperature (kg/m^3).
* @param specificHeat_ Specific heat at temperature (J/kg*K).
* @param conductivity_ Conductivity at temperature (w/m*K).
* @param emissivity_ Emmisivity at temperature.
* @param heatOfCombustion_ Heat of combustion at temperature (J/kg).
* @return
*/
void Material::push_temperatureRow(double temperature_, double density_, double specificHeat_, double conductivity_, double emissivity_, double heatOfCombustion_)
{
    propertyTable.push_back(std::make_tuple(temperature_, density_, specificHeat_, conductivity_, emissivity_, heatOfCombustion_));
    sort(propertyTable.begin(), propertyTable.end());
}


/**
* Clears all property data and material name.
* @return
*/
void Material::clearData()
{
    propertyTable.clear();
    name = "";
}


