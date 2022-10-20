#ifndef PARTICLE_H
#define PARTICLE_H

#include <stdio.h>
#include <stdexcept>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <memory>
#include <vector>
#include <ctime>
#include <cstring>

#include "solidMaterial.h"
#include "air.h"
#include "solidReaction.h"
#include "particleShape.h"
#include "slab.h"
#include "cylinder.h"
#include "sphere.h"

class Particle
{
public:

// Constructors
    Particle();
    virtual Particle* clone() const = 0;
    Particle(const Particle& rhs);
    Particle& operator=(const Particle& rhs);

// Default destructor
    virtual ~Particle();

// Public attributes:
 
    // State of the particle
    enum eState{
        preheat,
        drying,
        pyrolysing,
        flaming,
        charring,
        glowing,
        ashed,
        consumed,
        empty
    };
    eState state;

    // Burning rate model
    enum eBurningRateModel{
        uniformRateResidenceTime,           // Assumes energy is uniformly released over the residence time
        firstOrderArrheniusOneBreakPoint,   // The model from Tinney 1965
        thermalPyrolysisOnly,               // Arhennius reaction for thermal pyrolysis
        pyrolysisCharring                   // Arhenius reactions for thermal and oxydative pyrolysis, and charring
    };
    eBurningRateModel burningRateModel; //model used to compute the energy release rate of the particle as it burns

    // Name of particle (optional, defaults to empty string)
    std::string name;                       

    // Particle shape name as string (can be sphere, cylinder, or slab)
    std::string shapeName;

    // Pointer to particle shape (can be sphere, cylinder, or slab)
    particleShape* shape; 

    // Current global time (s)
    const double* pCurrentTime;

    // Material properties
    solidMaterial* wetSolid;    //For wet solid properties
    solidMaterial* drySolid;    //For dry solid properties
    solidMaterial* Char;        //For char properties (capital because "char" is a keyword)
    solidMaterial* ash;         //For ash properties
    Air* air;                   //For gas properties  

    //reaction mechanism
    solidReaction* R1;    //drying
    solidReaction* R2;    //thermal pyrolysis
    solidReaction* R3;    //oxidative pyrolysis
    solidReaction* R4;    //char oxidation

    // Integration step size [s]
    double localTimeStepSize;           

    // Model output quantities (for public access)
    double particleSize;                // Half thickness or radius of the particle [m]
    double particleVol;                 // Volume of the particle [m3]
    double particleSurfToVolRatio;      // Particle surface to volume ratio [1/m]
    double projectedAreaRatio;          // Particle projected to surface area ratio [-]

    double particleMass;                // Mass of the particle [kg]
    double wetSolidMass;                // Mass of wet solid inside the particle [kg]
    double drySolidMass;                // Mass of dry solid inside the particle [kg]
    double charMass;                    // Mass of char inside the particle [kg]

    double globalMassReleased;          // Integrated mass released for global time-step [kg]
    double globalGasFuelReleased;       // Integrated mass of the gaseous fuel released for global time-step [kg]
    double globalMoistureReleased;      // Integrated mass of water vapor released for global time-step [kg]
    double globalCO2Released;           // Integrated mass of CO2 released for global time-step [kg]
    double globalHeatReleased;          // Integrated heat due to the heterogeneous reactions for global time-step [J]

    double globalMassLossRate;          // Total particle mass loss rate for global time-step [kg/s]
    double globalGasFuelReleaseRate;    // Gaseous fuel release rate for global time-step [kg/s]
    double globalMoistureReleaseRate;   // Moisture release rate for global time-step [kg/s]
    double globalCO2ReleaseRate;        // CO2 release rate for global time-step [kg/s]
    
    double globalHeatReleaseRate;       // Heat consumption/production due to the heterogeneous reactions for global time-step [J/s]

    double globalR1reactionRate;        // Reactoin rate of drying for global time-step [kg/s]
    double globalR2reactionRate;        // Reactoin rate of thermal pyrolysis for global time-step [kg/s]
    double globalR3reactionRate;        // Reactoin rate of oxidative pyrolysis time-step [kg/s]
    double globalR4reactionRate;        // Reactoin rate of char oxidation time-step [kg/s]

    double surfaceTemp;                 // Temperature of the surface [K]
    double coreTemp;                    // Temperature of the core [K]
    double surfaceO2MassFrac;           // Mass fraction of O2 at the surface [-]
    double surfaceHeatFlux;             // Net surface heat flux [W/m2]
    double surfaceHeatFluxConv;         // Convective surface heat flux [W/m2]
    double surfaceHeatFluxRad;          // Radiative surface heat flux [W/m2]
    double surfaceMassFlux;             // Net surface mass flux [kg/s/m2]
    double h_conv;                      // Convective heat transfer coefficient of the particle [W/m2/K]
    double h_mass;                      // Convective mass transfer coefficient [kg/m2/s]
    double dragCoeff;                   // Drag force oefficient of the particle [-]
    double surfaceEmissivity;           // Particle surface emissivity [-]
    double ignitionTime;                // Fuel particle ignition time [s] (TODO)
    double outTime;                     // Fuel particle burned out time [s] (TODO)
    double ignitionTemp;                // Ignition temperature used in some combustion methods [K] (TODO)
    double flamingFraction;             // Fraction of total initial dry mass (including mineralContent) that participates in flaming combustion (TODO)
    double flamingFractionNoMinerals;   // Fraction of dry mass that participates in combustion (if no minerals were present)
    double residenceTime;               // Fuel particle residence time (=duration of the particle burning simulation) [s] (TODO)
    double mineralContent;              // Mass fraction of minerals.


// Public member Functions:

    //set geometry and dimensions
    void setGeometry(std::string geometry_, double radius_, double length, double width);

    //set burning rate model
    void setBurningRateModel(eBurningRateModel burningRateModel_, solidReaction& drying_, solidReaction& thermalPyrolysis_, solidReaction& oxydativePyrolysis_, solidReaction& charring_);

    //set material properties
    void setMaterials(Air& air_, solidMaterial& wetSolid_, solidMaterial& drySolid_, solidMaterial& Char_, solidMaterial& ash_);

    //computes mineral fraction from reaction parameters
    double getMineralFraction();

    //function initializing particle
    virtual void initialize(
                            double Temp_0,
                            double wetSolidVolFraction_0,
                            double drySolidVolFraction_0,
                            double charVolFraction_0,
                            double ashVolFraction_0,
                            double particleO2MassFraction_0,
                            double particlePressure_0
                            ) = 0;

    //step forward one (global) time step given local gas temperature, velocity, radiant flux and oxygen fraction
    virtual void stepForward(
                            const double globalTimeStepSize, 
                            const double externalTemperature, 
                            const double externalVelocity,
                            const double externalO2MassFrac,
                            const double externalPressure,
                            const double externalIrradiation
                            ) = 0;

    //clean up particle and destroy any dynamically allocated memory
    virtual void destroy() = 0;

    //write data at global time
    virtual void writeCoordLine(ofstream& outfile) = 0;     //write one line of coordinate output to outfile ()
    virtual void writeTempLine(ofstream& outfile) = 0;      //write one line of temperature output to outfile (temperature of each cell)
    virtual void writeWetSolidLine(ofstream& outfile) = 0;  //write one line of wet solid output to outfile (wet solid vol. frac. of each cell)
    virtual void writeDrySolidLine(ofstream& outfile) = 0;  //write one line of dry solid vol. fraction output to outfile (dry solid vol. frac. of each cell)
    virtual void writeCharLine(ofstream& outfile) = 0;      //write one line of char vol. fraction output to outfile (char vol. frac. of each cell)
    virtual void writeAshLine(ofstream& outfile) = 0;       //write one line of ash vol. fraction output to outfile (char vol. frac. of each cell)
    virtual void writeO2Line(ofstream& outfile) = 0;        //write one line of O2 mass fraction output to outfile (O2 mass frac. of each cell)
    virtual void writePressureLine(ofstream& outfile) = 0;  //write one line of pressure output to outfile (gauge pressure of each cell)


protected:

    static constexpr double sigma = 5.67e-8;     //Stefan-Boltzmann constant [W/m2/K4]
    static constexpr double R = 8.3144598;       //Gas constant [J/mol/K]

    //calculate minimum chemical timescale from used reaction mechanism
    double getChemicalTimescale(double Temp_); 

    //calculate minimum diffusion timescale
    double getDiffusionTimescale(double cellSize_, double Temp_); 

private:


};

#endif // PARTICLE_H
