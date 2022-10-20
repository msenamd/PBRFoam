#include "particle_lumped.h"

particle_lumped::particle_lumped() : Particle()
{
	localTime = 0.0;
    initialSolidVolFraction    = 0.0;
    currentLocalTimeStep = 0;
    numLocalTimeSteps  = 0;
	Temp = 0.0;
    volume = 0.0;
    initialVolume = 0.0;
	h_conv = 0.0;
	rho_p = 0.0;
    charVolFraction = 0.0;
    moistureVolFraction = 0.0;
    solidVolFraction = 0.0;
	q_surf = 0.0;
	eps_surf = 0.0;
    massLossRate  = 0.0;
	MLRpuv  = 0.0;
    gasFuelReleasedForGlobalTimeStep =0.0;

    volume_old = 0.0;
    radius_old = 0.0;
    Temp_old = 0.0;
    rho_p_old = 0.0;
    moistureVolFraction_old = 0.0;
    solidVolFraction_old =0.0;
    charVolFraction_old = 0.0;
    q_surf_old = 0.0;
}

particle_lumped::~particle_lumped()
{
    destroy();
}

/**
* Makes a copy of a particle_lumped allocated with new and returns a pointer to it.
* @param rhs particle_lumped to copy.
* @return
*/
particle_lumped* particle_lumped::clone() const
{
   return new particle_lumped(*this);
}


particle_lumped::particle_lumped(const particle_lumped& rhs) : Particle(rhs)
{
    //copy ctor
	localTime = rhs.localTime;
    initialSolidVolFraction   = rhs.initialSolidVolFraction;
    currentLocalTimeStep = rhs.currentLocalTimeStep;
    numLocalTimeSteps = rhs.numLocalTimeSteps;
	Temp = rhs.Temp;
    volume = rhs.volume;
    initialVolume = rhs.initialVolume;
	h_conv = rhs.h_conv;
	rho_p = rhs.rho_p;
    charVolFraction = rhs.charVolFraction;
    moistureVolFraction = rhs.moistureVolFraction;
    solidVolFraction = rhs.solidVolFraction;
	q_surf = rhs.q_surf;
	eps_surf = rhs.eps_surf;
    massLossRate = rhs.massLossRate;
	MLRpuv = rhs.MLRpuv;
    gasFuelReleasedForGlobalTimeStep = rhs.gasFuelReleasedForGlobalTimeStep;

    volume_old = rhs.volume_old;
    radius_old = rhs.radius_old;
    Temp_old = rhs.Temp_old;
    rho_p_old = rhs.rho_p_old;
    moistureVolFraction_old = rhs.moistureVolFraction_old;
    solidVolFraction_old = rhs.solidVolFraction_old;
    charVolFraction_old = rhs.charVolFraction_old;
    q_surf_old = rhs.q_surf_old;
}

particle_lumped& particle_lumped::operator=(const particle_lumped& rhs)
{
    if (&rhs != this)
    {
        localTime = rhs.localTime;
        initialSolidVolFraction   = rhs.initialSolidVolFraction;
        currentLocalTimeStep = rhs.currentLocalTimeStep;
        numLocalTimeSteps = rhs.numLocalTimeSteps;
        Temp = rhs.Temp;
        volume = rhs.volume;
        initialVolume = rhs.initialVolume;
        h_conv = rhs.h_conv;
        rho_p = rhs.rho_p;
        charVolFraction = rhs.charVolFraction;
        moistureVolFraction = rhs.moistureVolFraction;
        solidVolFraction = rhs.solidVolFraction;
        q_surf = rhs.q_surf;
        eps_surf = rhs.eps_surf;
        massLossRate = rhs.massLossRate;
        MLRpuv = rhs.MLRpuv;
        gasFuelReleasedForGlobalTimeStep = rhs.gasFuelReleasedForGlobalTimeStep;

        volume_old = rhs.volume_old;
        radius_old = rhs.radius_old;
        Temp_old = rhs.Temp_old;
        rho_p_old = rhs.rho_p_old;
        moistureVolFraction_old = rhs.moistureVolFraction_old;
        solidVolFraction_old = rhs.solidVolFraction_old;
        charVolFraction_old = rhs.charVolFraction_old;
        q_surf_old = rhs.q_surf_old;
    };

    return *this;
}


/**
* Initialize variables.
* @param initialTemp_ Initial particle temperature (K).
* @param initialPercentMoistureContent_ Initial percent moisture content based on a dry mass fraction (%).
* @param initialPercentCharContent_ Initial percent char content based on a dry mass fraction (%).
* @return
*/
void particle_lumped::initialize(double initialTemp_, double initialPercentMoistureContent_, double initialPercentCharContent_)
{

	state = burning;
    shape->air = air;

    Temp  = initialTemp_;
    surfaceTemp = initialTemp_;
    coreTemp = initialTemp_;

    initialVolume= shape->get_volume();
    volume=initialVolume;

    double initialDryMass, initialVirginSolidMass, initialCharMass, initialMoistureMass, initialParticleMass,
            initialVirginSolidMassFraction, initialCharMassFraction, initialMoistureMassFraction, initialParticleMassDensity;

    initialPercentMoistureContent = initialPercentMoistureContent_;

    initialVirginSolidMass = (1.0-initialPercentCharContent_/100.0)*fuel->get_density(initialTemp_) * initialVolume;
    initialCharMass = initialPercentCharContent_/100.0*Char->get_density(initialTemp_)*initialVolume;
    initialDryMass = initialVirginSolidMass + initialCharMass;
    initialMoistureMass = initialPercentMoistureContent/100.0*initialDryMass;
    initialParticleMass = initialDryMass + initialMoistureMass;

    initialVirginSolidMassFraction = initialVirginSolidMass/initialParticleMass;
    initialVirginSolidMassFraction = max(0.0,min(1.0,initialVirginSolidMassFraction)); // Safety: enforce 0 <= x_c <= 1

    initialCharMassFraction = initialCharMass/initialParticleMass;
    initialCharMassFraction = max(0.0,min(1.0,initialCharMassFraction)); // Safety: enforce 0 <= x_c <= 1

    initialMoistureMassFraction = initialMoistureMass/initialParticleMass;
    initialMoistureMassFraction = max(0.0,min(1.0,initialMoistureMassFraction)); // Safety: enforce 0 <= x_c <= 1

    initialParticleMassDensity = initialParticleMass/initialVolume;

    initialSolidVolFraction = initialParticleMassDensity/fuel->get_density(initialTemp_) * initialVirginSolidMassFraction;
    initialSolidVolFraction   = max(0.0,min(1.0,initialSolidVolFraction)); // Safety: enforce 0 <= x_c <= 1
    solidVolFraction  = initialSolidVolFraction;

    moistureVolFraction   = initialParticleMassDensity/moisture->get_density(initialTemp_) * initialMoistureMassFraction;
    moistureVolFraction   = max(0.0,min(1.0,moistureVolFraction)); // Safety: enforce 0 <= x_c <= 1

    charVolFraction   = initialParticleMassDensity/Char->get_density(initialTemp_) * initialCharMassFraction;
    charVolFraction   = max(0.0,min(1.0,charVolFraction)); // Safety: enforce 0 <= x_c <= 1

    rho_p = initialParticleMassDensity;
    //rho_p = moisture->get_density(surfaceTemp)*moistureVolFraction + fuel->get_density(surfaceTemp)*solidVolFraction + Char->get_density(surfaceTemp)*charVolFraction;

    localTimeStepSize    = getChemicalTimescale() * 0.1; //Enforce (tauChemical/dt) >=10

    eps_surf = moisture->get_emissivity(Temp)*moistureVolFraction + fuel->get_emissivity(Temp)*solidVolFraction + Char->get_emissivity(Temp)*charVolFraction;

        //some check prints
    cout << "Particle fuel is: " << fuel->get_name() <<endl;
    cout << "Particle geometry is: " << shape << endl;
    cout << "Particle initial radius or half thickness (m): " << shape->currentSize << endl;
    cout << "Particle initial mass (kg): " << (rho_p*volume) << endl;
	cout << "----------------------------" << endl;
}

//calculates surface area to volume rate
//double particle_lumped::getSurfaceAreaToVolumeRatio()
//{
//    double areaToVolume=0.0;

//    if(shape == eGeometry::rectangle)
//    {
//        areaToVolume=1.0/(2.0*radius);
//    }
//    else if(shape == eGeometry::cylinder)
//    {
//        areaToVolume=2.0/radius;
//    }
//    else if(shape == eGeometry::sphere)
//    {
//        areaToVolume=3.0/radius;
//    }
//    return areaToVolume;
//}

//destroy any dynamically allocated memory
void particle_lumped::destroy()
{

}

//writes particle temperature distribution to a line of the file
void particle_lumped::writeTempLine(FILE* outfile)
{
    fprintf(outfile, "%.6f\t%.6f\n", *pCurrentTime, Temp);
}

//writes vol fraction of moisture to a line of the file
void particle_lumped::writeMoistureLine(FILE* outfile)
{
    fprintf(outfile, "%.6f\t%.6f\n", *pCurrentTime, moistureVolFraction);
}

//writes vol fraction of virgin solid to a line of the file
void particle_lumped::writeSolidLine(FILE* outfile)
{
    fprintf(outfile, "%.6f\t%.6f\n", *pCurrentTime, solidVolFraction);
}

//writes vol fraction of char to a line of the file
void particle_lumped::writeCharLine(FILE* outfile)
{
    fprintf(outfile, "%.6f\t%.6f\n", *pCurrentTime, charVolFraction);
}

void particle_lumped::writeCoordLine(FILE* outfile)
{
    fprintf(outfile, "%.6f\t%.6f\n", *pCurrentTime, shape->currentSize);
}

/**
* Steps forward one time step, computing particle heating/cooling.
* @param timeStepSize Global time step size (s).
* @param gasTemperature Gas temperature around particle (K).
* @param gasVelocity Gas velocity (m/s).
* @param irradiation Irradiation on particle (w/m^2).
* @param oxygenVolFraction Oxygen concentration of gas as volume fraction.
* @return
*/
void particle_lumped::stepForward(const double globalTimeStepSize,
                                  const double gasTemperature, const double gasVelocity,
                                  const double irradiation , const double oxygenVolFraction)
{
    currentLocalTimeStep      = 0;
    localTimeStepSize     = min(localTimeStepSize,globalTimeStepSize);  //in case of timeStepSize>dt
    numLocalTimeSteps    = round(globalTimeStepSize/localTimeStepSize);

    if((numLocalTimeSteps*localTimeStepSize) !=globalTimeStepSize){    //in case of dt*n_f not equal timeStepSize
        numLocalTimeSteps = numLocalTimeSteps+1;
        localTimeStepSize  = globalTimeStepSize/numLocalTimeSteps;
    }

    localTime   = 0.0;
    gasFuelReleasedForGlobalTimeStep =0.0;
    massLossRate = 0.0;

    // local time loop
    while((currentLocalTimeStep!=numLocalTimeSteps) && (state!=burned))
    {
        currentLocalTimeStep = currentLocalTimeStep+1;
        localTime = localTime + localTimeStepSize;

        volume_old= volume;
        radius_old= shape->currentSize;

        Temp_old=Temp;
        rho_p_old=rho_p;
        moistureVolFraction_old=moistureVolFraction;
        solidVolFraction_old=solidVolFraction;
        charVolFraction_old=charVolFraction;
        h_conv = shape->get_convectiveHeatTransferCoefficient(gasTemperature,Temp_old,gasVelocity,radius_old);
        //h_conv = get_h(gasTemperature,Temp_old,gasVelocity,radius_old);
        q_surf_old = eps_surf*irradiation - eps_surf*sigma*pow(Temp_old, 4.0) + h_conv*(gasTemperature - Temp_old);
        surfaceTemp = Temp;
        coreTemp = Temp;


        //update species volume fractions, particle volume, and particle currentSize
        mass_conservation_TTC(oxygenVolFraction);

        // update particle temperature
        energy_conservation_TTC(oxygenVolFraction);

        //update particle density
        rho_p = moisture->get_density(Temp)*moistureVolFraction + fuel->get_density(Temp)*solidVolFraction + Char->get_density(Temp)*charVolFraction;

        //update net surface heat flux
        eps_surf = moisture->get_emissivity(Temp)*moistureVolFraction + fuel->get_emissivity(Temp)*solidVolFraction + Char->get_emissivity(Temp)*charVolFraction;
        h_conv = shape->get_convectiveHeatTransferCoefficient(gasTemperature,Temp,gasVelocity,shape->currentSize);
        //h_conv = get_h(gasTemperature,Temp,gasVelocity,radius);
        q_surf = eps_surf*irradiation - eps_surf*sigma*pow(Temp, 4.0) + h_conv*(gasTemperature - Temp);

        //calculate volumetric mass loss rate at this local time-step
        MLRpuv = moisture->get_density(Temp) *moistureVolFraction *R1.get_A()*exp(-R1.get_Ta()/Temp)
               + fuel->get_density(Temp)*solidVolFraction * R2.get_A()*exp(-R2.get_Ta()/Temp) *(1-charYield)
               + Char->get_density(Temp) *charVolFraction *oxygenVolFraction *R3.get_A()*exp(-R3.get_Ta()/Temp);

        massLossRate = MLRpuv * volume;  //[kg/s]

        //calculate total gaseous fuel released in kg (integration over simulation time)
        gasFuelReleasedForGlobalTimeStep = gasFuelReleasedForGlobalTimeStep
                        + localTimeStepSize*volume* ( fuel->get_density(Temp)*solidVolFraction * R2.get_A()*exp(-R2.get_Ta()/Temp) *(1-charYield)
                                  + Char->get_density(Temp) *charVolFraction *oxygenVolFraction *R3.get_A()*exp(-R3.get_Ta()/Temp) );

        //update particle state
        if((charYield!=0.0) && (solidVolFraction <= 0.001))
        {
            state=burned;
			cout << "STATE (0-burning 1-burned): " << state << endl;
        }
        if((charYield==0.0) && ((volume/initialVolume) <= 0.002))  //volume shrinks to 0.2%
        {
            state=burned;
			cout << "STATE (0-burning 1-burned): " << state << endl;
        }
       

    } //end time loop

    //correction to account for full rectangular particle (ie. both top and bottom sides)
    gasFuelReleasedForGlobalTimeStep = shape->correctGasFuelReleasedForShape(gasFuelReleasedForGlobalTimeStep);

    //shape->currentSize = radius;
    outTime = localTime;

    cout << "particle size (m): " << shape->currentSize << endl;
    cout << "gaseous fuel released (kg): " << gasFuelReleasedForGlobalTimeStep << endl;
	cout << "---------------------------------" << endl;
}

void particle_lumped::mass_conservation_TTC(double x_O2_g)
{

    double xm_times_dV, xvs_times_dV, store , xc_times_dV;

    xm_times_dV  = (moistureVolFraction_old*volume_old)*exp( -localTimeStepSize*R1.get_A()*exp(-R1.get_Ta()/Temp_old) );
    xvs_times_dV = (solidVolFraction_old*volume_old)/(1+localTimeStepSize*R2.get_A()*exp(-R2.get_Ta()/Temp_old));

    if( charYield !=0.0 )
    {
        store = xvs_times_dV
              *localTimeStepSize*R2.get_A()*exp(-R2.get_Ta()/Temp_old)*(charYield*fuel->get_density(Temp_old)/Char->get_density(Temp_old));

        xc_times_dV = (charVolFraction_old*volume_old + store)
                     /(1+localTimeStepSize*x_O2_g*R3.get_A()*exp(-R3.get_Ta()/Temp_old));

    }
    else
    {
        xc_times_dV = 0.0;
    }


    volume = xm_times_dV + xvs_times_dV + xc_times_dV;
    moistureVolFraction  = xm_times_dV /volume;
    moistureVolFraction  = max(0.0,min(1.0,moistureVolFraction));    // Safety: enforce 0 <= x_m  <= 1
    solidVolFraction = xvs_times_dV/volume;
    solidVolFraction = max(0.0,min(1.0,solidVolFraction));   // Safety: enforce 0 <= x_vs <= 1
    charVolFraction  = 1-moistureVolFraction-solidVolFraction;
    charVolFraction  = max(0.0,min(1.0,charVolFraction));    //Safety: enforce 0 <= x_c  <= 1

    shape->currentSize = shape->computeSizeFromVolume(volume);
}


void particle_lumped::energy_conservation_TTC(double x_O2_g)
{
    double rho_times_cp, Qdotp;

    rho_times_cp = moisture->get_density(Temp_old)*moisture->get_specificHeat(Temp_old)*moistureVolFraction_old
                 + fuel->get_density(Temp_old)*fuel->get_specificHeat(Temp_old)*solidVolFraction_old
                 + Char->get_density(Temp_old)*Char->get_specificHeat(Temp_old)*charVolFraction_old;

    // Volumetric rate of heat production/consumption [W/m3]
    Qdotp = moisture->get_density(Temp_old) * moistureVolFraction_old*R1.get_A()*exp(-R1.get_Ta()/Temp_old)*R1.get_deltaH()
          + fuel->get_density(Temp_old)*solidVolFraction_old*R2.get_A()*exp(-R2.get_Ta()/Temp_old)*R2.get_deltaH()*(1-charYield)
          + Char->get_density(Temp_old) * charVolFraction_old*x_O2_g *R3.get_A()*exp(-R3.get_Ta()/Temp_old)*R3.get_deltaH();

    // Energy conservation statement
    Temp = Temp_old + (Qdotp * localTimeStepSize/rho_times_cp)
        + (q_surf_old * shape->get_surfaceArea(radius_old) * localTimeStepSize/rho_times_cp/volume_old);
}



