/*
 * solidReaction.h
 * A class of chemical reactions
 * All properties must me listed here
 * Different reaction mechanisms are to be implemented
 *
 *  Created on: October 7, 2019
 *      Author: mmahmed15
 */
#include "solidReaction.h"

solidReaction::solidReaction()
{
    A = 0.0;
    Ea = 0.0;
    n = 0.0;
    nO2 = 0.0;
    deltaH = 0.0;
    productYield = 0.0;
    O2Yield = 0.0;
}

solidReaction::~solidReaction()
{

}

solidReaction::solidReaction(const solidReaction& rhs)
{
    A = rhs.A;
    Ea = rhs.Ea;
    n = rhs.n;
    nO2 = rhs.nO2;
    deltaH = rhs.deltaH;
    productYield = rhs.productYield;
    O2Yield = rhs.O2Yield;
}

solidReaction& solidReaction::operator=(const solidReaction& rhs)
{
    if (&rhs != this) {
        A = rhs.A;
        Ea = rhs.Ea;
        n = rhs.n;
        nO2 = rhs.nO2;
        deltaH = rhs.deltaH;
        productYield = rhs.productYield;
        O2Yield = rhs.O2Yield;
    }; // handle self assignment
    //assignment operator
    return *this;
}

void solidReaction::set_reaction(const double& A_, const double& Ea_, const double& n_, const double& nO2_,
                                 const double& deltaH_, const double& productYield_, const double& O2Yield_)
{
    A = A_;
    Ea = Ea_;
    n = n_;
    nO2 = nO2_;
    deltaH = deltaH_;
    productYield = productYield_;
    O2Yield = O2Yield_;
}

void solidReaction::set_knownReaction(std::string reactionName, const double& productYield_)
{
    if(reactionName == "drying")  
    {
        // Source: Lautenberger & Fernandez-Pello(2009) Combust. Flame 156 : 1503 - 1513
        A = 4.29e+3; 
        Ea = 43.8e+3;
        n = 0.99;
        nO2 = 0.0;
        deltaH = -2410e+3;
        productYield = productYield_;
        O2Yield = 0.0;
    }
    else if(reactionName == "pinewood_thermalPyrolysis") 
    {
        // Source: Anca-Couce et al. (2012) Combust. Flame 159 : 1708 - 1719
        // Source(DeltaH) : Lautenberger & Fernandez-Pello(2009) Combust. Flame 156 : 1503 - 1513
        A = std::pow(10.0, 6.34);
        Ea = 105e+3;
        n = 0.87;
        nO2 = 0.0;
        deltaH = -533e+3;
        productYield = productYield_;
        O2Yield = 0.0;
    }
    else if(reactionName == "pinewood_oxidativePyrolysis")
    {
        // Source: Anca-Couce et al. (2012) Combust. Flame 159 : 1708 - 1719
        // Source(DeltaH) : Lautenberger & Fernandez-Pello(2009) Combust. Flame 156 : 1503 - 1513
        A = std::pow(10.0, 8.72);
        Ea = 127e+3;
        n = 0.63;
        nO2 = 0.72;
        deltaH = +994e+3;
        productYield = productYield_;
        O2Yield = 0.1 * (1.0 - productYield_);
    }
    else if (reactionName == "pinewood_charOxidation")
    {
        // Source: Anca-Couce et al. (2012) Combust. Flame 159 : 1708 - 1719
        // Source(DeltaH) : Lautenberger & Fernandez-Pello(2009) Combust. Flame 156 : 1503 - 1513
        A = std::pow(10.0, 6.55);
        Ea = 124e+3;
        n = 0.56;
        nO2 = 0.68;
        deltaH = +37700e+3;
        productYield = productYield_;
        O2Yield = 2.0 * (1.0 - productYield_);
    }
    else if (reactionName == "passive")
    {
        A = 0.0;
        Ea = 0.0;
        n = 0.0;
        nO2 = 0.0;
        deltaH = 0.0;
        productYield = 0.0;
        O2Yield = 0.0;
    }
    else
    {
        cout << "Error: Invalid reaction: " << reactionName << endl;
        throw exception();
    }
}


double solidReaction::get_A()
{
    return A;
}

double solidReaction::get_Ta()
{
    return Ea/8.3145;
}

double solidReaction::get_n()
{
    return n;
}

double solidReaction::get_nO2()
{
    return nO2;
}

double solidReaction::get_deltaH()
{
    return deltaH;
}

double solidReaction::get_productYield()
{
    return productYield;
}

double solidReaction::get_O2Yield()
{
    return O2Yield;
}
