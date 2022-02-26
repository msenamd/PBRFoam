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

bool solidReaction::set_reaction(std::string reactionName, double productYield_)
{
    if(reactionName == "drying")  
    {
        A = 4.29e+3;
        Ea = 43.8e+3;
        n = 0.99;
        nO2 = 0.0;
        deltaH = -2410e+3;
        productYield = productYield_;
        O2Yield = 0.0;
        return true;
    }
    else if(reactionName == "thermal pyrolysis") 
    {
        A = 3.29e+9;
        Ea = 135e+3;
        n = 4.78;
        nO2 = 0.0;
        deltaH = -533e+3;
        productYield = productYield_;
        O2Yield = 0.0;
        return true;
    }
    else if(reactionName == "oxidative pyrolysis")
    {
        A = 6.00e+9;
        Ea = 124.2e+3;
        n = 4.99;
        nO2 = 1.16;
        deltaH = +994e+3;
        productYield = productYield_;
        O2Yield = 0.1 * (1.0 - productYield_);
        return true;
    }
    else if (reactionName == "char oxidation")
    {
        A = 9.79e+13;
        Ea = 192.4e+3;
        n = 1.86;
        nO2 = 1.04;
        deltaH = +37700e+3;
        productYield = productYield_;
        O2Yield = 2.0 * (1.0 - productYield_);
        return true;
    }
    else if (reactionName == "passive solidReaction")
    {
        A = 0.0;
        Ea = 0.0;
        n = 0.0;
        nO2 = 0.0;
        deltaH = 0.0;
        productYield = 0.0;
        O2Yield = 0.0;
        return true;
    }
    else
    {
        cout << "Error: Invalid reaction" << endl;
        return false;   //unknown name
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
