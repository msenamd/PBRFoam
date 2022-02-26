#include "particle.h"

Particle::Particle()
{
    state = empty;
    burningRateModel = uniform;

    name = "";    
    shapeName = "";
    shape = NULL;

    pCurrentTime = NULL;

    wetSolid = NULL;
    drySolid = NULL;
    Char = NULL;
    ash = NULL;
    air = NULL;

    R1 = NULL;
    R2 = NULL;
    R3 = NULL;
    R4 = NULL;

    particleSize = 0.0;                
    particleVol = 0.0;                 
    particleSurfToVolRatio = 0.0;   
    projectedAreaRatio = 0.0;

    particleMass = 0.0;

    globalMassReleased = 0.0;          
    globalGasFuelReleased = 0.0; 
    globalMoistureReleased = 0.0;
    globalCO2Released = 0.0;
    globalHeatReleased = 0.0;

    globalMassLossRate = 0.0;   
    globalGasFuelReleaseRate = 0.0;
    globalMoistureReleaseRate = 0.0;
    globalCO2ReleaseRate = 0.0;
    globalHeatReleaseRate = 0.0;

    surfaceTemp = 0.0;    
    coreTemp = 0.0;
    surfaceO2MassFrac = 0.0;
    surfaceHeatFlux = 0.0;
    surfaceMassFlux = 0.0;
    h_conv = 0.0;
    h_mass = 0.0;
    dragCoeff = 0.0;
    surfaceEmissivity = 0.0;
    ignitionTime = 0.0;                
    outTime = 0.0;                     
    ignitionTemp = 0.0;                
    flamingFraction = 0.0;             
    residenceTime = 0.0;

    localTimeStepSize = 0.0;
}

Particle::~Particle()
{
    if(shape)
        delete shape;
    shape = NULL;
}

Particle::Particle(const Particle& rhs)
{
    //copy ctor

    state = rhs.state;
    burningRateModel = rhs.burningRateModel;

    name = rhs.name;
    shapeName = rhs.shapeName;
    if(shape)
        delete shape;
    shape = rhs.shape->clone();

    pCurrentTime = rhs.pCurrentTime;

    if (wetSolid)
        delete wetSolid;
    wetSolid = new solidMaterial(*rhs.wetSolid);

    if (drySolid)
        delete drySolid;
    drySolid = new solidMaterial(*rhs.drySolid);

    if (Char)
        delete Char;
    Char = new solidMaterial(*rhs.Char);

    if (ash)
        delete ash;
    ash = new solidMaterial(*rhs.ash);

    if (R1)
        delete R1;
    R1 = new solidReaction(*rhs.R1);

    if (R2)
        delete R2;
    R2 = new solidReaction(*rhs.R2);

    if (R3)
        delete R3;
    R3 = new solidReaction(*rhs.R3);

    if (R4)
        delete R4;
    R4 = new solidReaction(*rhs.R4);

    particleSize = rhs.particleSize;
    particleVol = rhs.particleVol;
    particleSurfToVolRatio = rhs.particleSurfToVolRatio;
    projectedAreaRatio = rhs.projectedAreaRatio;

    particleMass = rhs.particleMass;

    globalMassReleased = rhs.globalMassReleased;
    globalGasFuelReleased = rhs.globalGasFuelReleased;
    globalMoistureReleased = rhs.globalMoistureReleased;
    globalCO2Released = rhs.globalCO2Released;
    globalHeatReleased = rhs.globalHeatReleased;

    globalMassLossRate = rhs.globalMassLossRate;
    globalGasFuelReleaseRate = rhs.globalGasFuelReleaseRate;
    globalMoistureReleaseRate = rhs.globalMoistureReleaseRate;
    globalCO2ReleaseRate = rhs.globalCO2ReleaseRate;
    globalHeatReleaseRate = rhs.globalHeatReleaseRate;

    surfaceTemp = rhs.surfaceTemp;
    coreTemp = rhs.coreTemp;
    surfaceO2MassFrac = rhs.surfaceO2MassFrac;
    surfaceHeatFlux = rhs.surfaceHeatFlux;
    surfaceMassFlux = rhs.surfaceMassFlux;
    h_conv = rhs.h_conv;
    h_mass = rhs.h_mass;
    dragCoeff = rhs.dragCoeff;
    surfaceEmissivity = rhs.surfaceEmissivity;
    ignitionTime =rhs.ignitionTime;
    outTime = rhs.outTime;
    ignitionTemp = rhs.ignitionTemp;
    flamingFraction = rhs.flamingFraction;
    residenceTime = rhs.residenceTime;

    localTimeStepSize = rhs.localTimeStepSize;
}

Particle& Particle::operator=(const Particle& rhs)
{
    if (&rhs != this){

        state = rhs.state;
        burningRateModel = rhs.burningRateModel;

        name = rhs.name;
        shapeName = rhs.shapeName;
        if (shape)
            delete shape;
        shape = rhs.shape->clone();

        pCurrentTime = rhs.pCurrentTime;

        if (wetSolid)
            delete wetSolid;
        wetSolid = new solidMaterial(*rhs.wetSolid);

        if (drySolid)
            delete drySolid;
        drySolid = new solidMaterial(*rhs.drySolid);

        if (Char)
            delete Char;
        Char = new solidMaterial(*rhs.Char);

        if (ash)
            delete ash;
        ash = new solidMaterial(*rhs.ash);

        if (R1)
            delete R1;
        R1 = new solidReaction(*rhs.R1);

        if (R2)
            delete R2;
        R2 = new solidReaction(*rhs.R2);

        if (R3)
            delete R3;
        R3 = new solidReaction(*rhs.R3);

        if (R4)
            delete R4;
        R4 = new solidReaction(*rhs.R4);

        particleSize = rhs.particleSize;
        particleVol = rhs.particleVol;
        particleSurfToVolRatio = rhs.particleSurfToVolRatio;
        projectedAreaRatio = rhs.projectedAreaRatio;

        particleMass = rhs.particleMass;

        globalMassReleased = rhs.globalMassReleased;
        globalGasFuelReleased = rhs.globalGasFuelReleased;
        globalMoistureReleased = rhs.globalMoistureReleased;
        globalCO2Released = rhs.globalCO2Released;
        globalHeatReleased = rhs.globalHeatReleased;

        globalMassLossRate = rhs.globalMassLossRate;
        globalGasFuelReleaseRate = rhs.globalGasFuelReleaseRate;
        globalMoistureReleaseRate = rhs.globalMoistureReleaseRate;
        globalCO2ReleaseRate = rhs.globalCO2ReleaseRate;
        globalHeatReleaseRate = rhs.globalHeatReleaseRate;

        surfaceTemp = rhs.surfaceTemp;
        coreTemp = rhs.coreTemp;
        surfaceO2MassFrac = rhs.surfaceO2MassFrac;
        surfaceHeatFlux = rhs.surfaceHeatFlux;
        surfaceMassFlux = rhs.surfaceMassFlux;
        h_conv = rhs.h_conv;
        h_mass = rhs.h_mass;
        dragCoeff = rhs.dragCoeff;
        surfaceEmissivity = rhs.surfaceEmissivity;
        ignitionTime = rhs.ignitionTime;
        outTime = rhs.outTime;
        ignitionTemp = rhs.ignitionTemp;
        flamingFraction = rhs.flamingFraction;
        residenceTime = rhs.residenceTime;

        localTimeStepSize = rhs.localTimeStepSize;

    } ; // handle self assignment
    //assignment operator
    return *this;
}

/**
* Estimate time-step based on heterogeneous chemistry.
* @param Temp_ cell temperature [K].
* @return chemical time-scale.
*/
double Particle::getChemicalTimescale(double Temp_){

    double tau_R1_min, tau_R2_min, tau_R3_min, tau_R4_min;

    double dt = 0.0;

    tau_R1_min = 1.0 / max(R1->get_A() * exp(-R1->get_Ta() / Temp_), 1e-6); 
    tau_R2_min = 1.0 / max(R2->get_A() * exp(-R2->get_Ta() / Temp_), 1e-6);
    tau_R3_min = 1.0 / max(R3->get_A() * exp(-R3->get_Ta() / Temp_), 1e-6);
    tau_R4_min = 1.0 / max(R4->get_A() * exp(-R4->get_Ta() / Temp_), 1e-6);

    dt = 0.1*min({ tau_R1_min, tau_R2_min, tau_R3_min, tau_R4_min });

    return dt;
}

/**
* Estimate time-step based on heat diffusion inside solid phase.
* @param cellSize_ cell size [m].
* @param Temp_ cell temperature [K].
* @return diffusion time-scale.
*/
double Particle::getDiffusionTimescale(double cellSize_, double Temp_) {

    double dt = 0.0;

    // Enforce Fourier number FO = diffusivity*dt/dx^2 ~0.5
    double maxDiffusivity = max({
        wetSolid->get_conductivity(Temp_) / wetSolid->get_bulkDensity(Temp_) * (1- wetSolid->get_porosity(Temp_)) / wetSolid->get_specificHeat(Temp_),
        drySolid->get_conductivity(Temp_) / drySolid->get_bulkDensity(Temp_) * (1 - drySolid->get_porosity(Temp_)) / drySolid->get_specificHeat(Temp_),
        Char->get_conductivity(Temp_) / Char->get_bulkDensity(Temp_) * (1 - Char->get_porosity(Temp_)) / Char->get_specificHeat(Temp_),
        ash->get_conductivity(Temp_) / ash->get_bulkDensity(Temp_) * (1 - ash->get_porosity(Temp_)) / ash->get_specificHeat(Temp_)
        });

    dt = 0.5 * pow(cellSize_, 2.0) / maxDiffusivity;

    return dt;
}

/**
* Sets particle to be a specific geomerty.
* @param geometry_ geometric shape.
* @param radius_ Radius or half thickness (m).
* @param length_ Length of cylinder or slab (m).
* @param width_ Width of slab (m).
* @return
*/
void Particle::setGeometry(std::string shapeName_, double& radius_, double& length_, double& width_)
{
    if (shapeName_ == "slab")
    {
        shape = new slab(radius_, length_, width_);
    }
    else if (shapeName_ == "cylinder")
    {
        shape = new cylinder(radius_, length_);
    }
    else if (shapeName_ == "sphere")
    {
        shape = new sphere(radius_);
    }
    else
    {
        cout << "Error: unknown shape name. Available shapes are: slab-cylinder-sphere" << endl;
        return;
    }

    shapeName = shapeName_;
}