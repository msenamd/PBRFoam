#ifndef REACTION_H
#define REACTION_H

#include <string>
#include <iostream>
using namespace std;
#include <fstream>
#include <vector>
#include <tuple>
#include <algorithm>


class solidReaction{

public:
    solidReaction();
    ~solidReaction();

    solidReaction(const solidReaction& rhs);
    solidReaction& operator=(const solidReaction& rhs);

    bool set_reaction(std::string reactionName, double productYield_);

    double get_A();
    double get_Ta();
    double get_n();
    double get_nO2();
    double get_deltaH();
    double get_productYield();
    double get_O2Yield();

private:

    double A;               //Pre-exponential factor [1/s]
    double Ea;              //Activation energy [J/mol]
    double n;               //Solid-phase reactant exponent [-]
    double nO2;             //Gas-phase oxygen exponent [-]
    double deltaH;          //Heat of reaction [J/kg]
    double productYield;    //Mass yield of the product [-]
    double O2Yield;         //Oxygen-to-reactant mass ratio [-]
};

#endif
