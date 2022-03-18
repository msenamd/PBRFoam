#include "particle_1D.h"

particle_1D::particle_1D() : Particle()
{
    numCells = 0;
    localTimeStepIndex = 0;
    finalTimeStepIndex = 0;
    initParticleSize = 0.0;
    surfaceConductivity = 0.0;
    surfacePorosity = 0.0;
    surfacePermeability = 0.0;
    surfaceDiffusivity = 0.0;
    surfaceGridSpacing = 0.0;
    localTime = 0.0;
    FO_L = 0.0;
    FO_R = 0.0;
    CFL_L = 0.0;
    CFL_R = 0.0;
    h_rad = 0.0;
    Bi = 0.0;

    iter = 0;
    iterMaxError = 0.0;
    maxIter = 1000;
    defaultMeshResolution = 100e-6;
    timeStepThreshold = 0.1;
    temperatureThreshold = 0.1;
    solidSpeciesThreshold = 0.01;
    O2Threshold = 0.01;
    pressureThreshold = 0.1;
    temperatureURF = 0.01;
    solidSpeciesURF = 0.0;
    O2URF = 0.01;
    pressureURF = 0.0;

    TempDifference_iter = 0.0;     
    wetSolidDifference_iter = 0.0;
    drySolidDifference_iter = 0.0;
    charDifference_iter = 0.0;
    ashDifference_iter = 0.0;
    O2Difference_iter = 0.0;
    pressureDifference_iter = 0.0;

    TempDifference_time = 0.0;
    wetSolidDifference_time = 0.0;
    drySolidDifference_time = 0.0;
    charDifference_time = 0.0;
    ashDifference_time = 0.0;
    O2Difference_time = 0.0;
    pressureDifference_time = 0.0;

    TempError = 0.0;
    wetSolidError = 0.0;
    drySolidError = 0.0;
    charError = 0.0;
    ashError = 0.0;
    O2Error = 0.0;
    pressureError = 0.0;

    numCells_old = 0;
    localTimeStepSize_old = 0.0;
    surfaceTemp_old = 0.0;

    surfaceTemp_oldIter = 0.0;

    surfaceTemp_newIter = 0.0;

    localMassLossRate = 0.0;
    localGasFuelReleaseRate = 0.0;
    localMoistureReleaseRate = 0.0;
    localCO2ReleaseRate = 0.0;
    localHeatReleaseRate = 0.0;

}

particle_1D::~particle_1D()
{
	destroy();
}

/**
* Makes a copy of a particle_1D allocated with new and returns a pointer to it.
* @param rhs particle_1D to copy.
* @return
*/
particle_1D* particle_1D::clone() const
{
    return new particle_1D(*this);
}

particle_1D::particle_1D(const particle_1D& rhs) : Particle(rhs)

{
	//copy ctor

    numCells = rhs.numCells;
    localTimeStepIndex = rhs.localTimeStepIndex;
    finalTimeStepIndex = rhs.finalTimeStepIndex;
    initParticleSize = rhs.initParticleSize;
    surfaceConductivity = rhs.surfaceConductivity;
    surfacePorosity = rhs.surfacePorosity;
    surfacePermeability = rhs.surfacePermeability;
    surfaceDiffusivity = rhs.surfaceDiffusivity;
    surfaceGridSpacing = rhs.surfaceGridSpacing;
    localTime = rhs.localTime;
    FO_L = rhs.FO_L;
    FO_R = rhs.FO_R;
    CFL_L = rhs.CFL_L;
    CFL_R = rhs.CFL_R;
    h_rad = rhs.h_rad;
    Bi = rhs.Bi;

    iter = rhs.iter;
    iterMaxError = rhs.iterMaxError;
    maxIter = rhs.maxIter;
    defaultMeshResolution = rhs.defaultMeshResolution;
    timeStepThreshold = rhs.timeStepThreshold;
    temperatureThreshold = rhs.temperatureThreshold;
    solidSpeciesThreshold = rhs.solidSpeciesThreshold;
    O2Threshold = rhs.O2Threshold;
    pressureThreshold = rhs.pressureThreshold;
    temperatureURF = rhs.temperatureURF;
    solidSpeciesURF = rhs.solidSpeciesURF;
    O2URF = rhs.O2URF;
    pressureURF = rhs.pressureURF;

    TempDifference_iter = rhs.TempDifference_iter;
    wetSolidDifference_iter = rhs.wetSolidDifference_iter;
    drySolidDifference_iter = rhs.drySolidDifference_iter;
    charDifference_iter = rhs.charDifference_iter;
    ashDifference_iter = rhs.ashDifference_iter;
    O2Difference_iter = rhs.O2Difference_iter;
    pressureDifference_iter = rhs.pressureDifference_iter;

    TempDifference_time = rhs.TempDifference_time;
    wetSolidDifference_time = rhs.wetSolidDifference_time;
    drySolidDifference_time = rhs.drySolidDifference_time;
    charDifference_time = rhs.charDifference_time;
    ashDifference_time = rhs.ashDifference_time;
    O2Difference_time = rhs.O2Difference_time;
    pressureDifference_time = rhs.pressureDifference_time;

    TempError = rhs.TempError;
    wetSolidError = rhs.wetSolidError;
    drySolidError = rhs.drySolidError;
    charError = rhs.charError;
    ashError = rhs.ashError;
    O2Error = rhs.O2Error;
    pressureError = rhs.pressureError;

    numCells_old = rhs.numCells_old;
    localTimeStepSize_old = rhs.localTimeStepSize_old;
    surfaceTemp_old = rhs.surfaceTemp_old;

    surfaceTemp_oldIter = rhs.surfaceTemp_oldIter;

    surfaceTemp_newIter = rhs.surfaceTemp_newIter;

    localMassLossRate = rhs.localMassLossRate;
    localGasFuelReleaseRate = rhs.localGasFuelReleaseRate;
    localMoistureReleaseRate = rhs.localMoistureReleaseRate;
    localCO2ReleaseRate = rhs.localCO2ReleaseRate;
    localHeatReleaseRate = rhs.localHeatReleaseRate;

    initWetSolidMass = rhs.initWetSolidMass;
    initDrySolidMass = rhs.initDrySolidMass;
    initCharMass = rhs.initCharMass;
}

particle_1D& particle_1D::operator=(const particle_1D& rhs)
{
    if (&rhs != this)
    {
        Particle::operator=(rhs);

        numCells = rhs.numCells;
        localTimeStepIndex = rhs.localTimeStepIndex;
        finalTimeStepIndex = rhs.finalTimeStepIndex;
        initParticleSize = rhs.initParticleSize;
        surfaceConductivity = rhs.surfaceConductivity;
        surfacePorosity = rhs.surfacePorosity;
        surfacePermeability = rhs.surfacePermeability;
        surfaceDiffusivity = rhs.surfaceDiffusivity;
        surfaceGridSpacing = rhs.surfaceGridSpacing;
        localTime = rhs.localTime;
        FO_L = rhs.FO_L;
        FO_R = rhs.FO_R;
        CFL_L = rhs.CFL_L;
        CFL_R = rhs.CFL_R;
        h_rad = rhs.h_rad;
        Bi = rhs.Bi;

        iter = rhs.iter;
        iterMaxError = rhs.iterMaxError;
        maxIter = rhs.maxIter;
        defaultMeshResolution = rhs.defaultMeshResolution;
        timeStepThreshold = rhs.timeStepThreshold;
        temperatureThreshold = rhs.temperatureThreshold;
        solidSpeciesThreshold = rhs.solidSpeciesThreshold;
        O2Threshold = rhs.O2Threshold;
        pressureThreshold = rhs.pressureThreshold;
        temperatureURF = rhs.temperatureURF;
        solidSpeciesURF = rhs.solidSpeciesURF;
        O2URF = rhs.O2URF;
        pressureURF = rhs.pressureURF;

        TempDifference_iter = rhs.TempDifference_iter;
        wetSolidDifference_iter = rhs.wetSolidDifference_iter;
        drySolidDifference_iter = rhs.drySolidDifference_iter;
        charDifference_iter = rhs.charDifference_iter;
        ashDifference_iter = rhs.ashDifference_iter;
        O2Difference_iter = rhs.O2Difference_iter;
        pressureDifference_iter = rhs.pressureDifference_iter;

        TempDifference_time = rhs.TempDifference_time;
        wetSolidDifference_time = rhs.wetSolidDifference_time;
        drySolidDifference_time = rhs.drySolidDifference_time;
        charDifference_time = rhs.charDifference_time;
        ashDifference_time = rhs.ashDifference_time;
        O2Difference_time = rhs.O2Difference_time;
        pressureDifference_time = rhs.pressureDifference_time;

        TempError = rhs.TempError;
        wetSolidError = rhs.wetSolidError;
        drySolidError = rhs.drySolidError;
        charError = rhs.charError;
        ashError = rhs.ashError;
        O2Error = rhs.O2Error;
        pressureError = rhs.pressureError;

        numCells_old = rhs.numCells_old;
        localTimeStepSize_old = rhs.localTimeStepSize_old;
        surfaceTemp_old = rhs.surfaceTemp_old;

        surfaceTemp_oldIter = rhs.surfaceTemp_oldIter;

        surfaceTemp_newIter = rhs.surfaceTemp_newIter;

        localMassLossRate = rhs.localMassLossRate;
        localGasFuelReleaseRate = rhs.localGasFuelReleaseRate;
        localMoistureReleaseRate = rhs.localMoistureReleaseRate;
        localCO2ReleaseRate = rhs.localCO2ReleaseRate;
        localHeatReleaseRate = rhs.localHeatReleaseRate;

        initWetSolidMass = rhs.initWetSolidMass;
        initDrySolidMass = rhs.initDrySolidMass;
        initCharMass = rhs.initCharMass;
    }
    return *this;
}

/**
* Set the parameters that control solution accuracy and stability.
* @param maxIter_: Maximum number of iterations allowed [-].
* @param timeStepThreshold_: Max. allowed time-step size for the particle [s].
* @param temperatureThreshold_: Max. value of temperature variation during dt [K].
* @param solidSpeciesThreshold_: Max. value of vol. variation fractions during dt [-].
* @param O2Threshold_: Max. value of vol. O2 mass fraction variation during dt [-].
* @param pressureThreshold_: Max. value of pressure variation during dt [Pa].
* @param temperatureURF_: Under-relaxation parameter for temperature [-].
* @param solidSpeciesURF_: Under-relaxation parameter for vol. fractions [-].
* @param O2URF_: Under-relaxation parameter for oxygen mass fraction [-].
* @param pressureURF_: Under-relaxation parameter for pressure [-].
* @return
*/
void particle_1D::setSolutionControl(
    const int maxIter_,
    const double defaultMeshResolution_,
    const double timeStepThreshold_,
    const double temperatureThreshold_,
    const double solidSpeciesThreshold_,
    const double O2Threshold_,
    const double pressureThreshold_,
    const double temperatureURF_,
    const double solidSpeciesURF_,
    const double O2URF_,
    const double pressureURF_
)
{
    maxIter = maxIter_;
    defaultMeshResolution = defaultMeshResolution_;
    timeStepThreshold = timeStepThreshold_;
    temperatureThreshold = temperatureThreshold_;
    solidSpeciesThreshold = solidSpeciesThreshold_;
    O2Threshold = O2Threshold_;
    pressureThreshold = pressureThreshold_;
    temperatureURF = temperatureURF_;
    solidSpeciesURF = solidSpeciesURF_;
    O2URF = O2URF_;
    pressureURF = pressureURF_;
}

/**
* Set computational grid (works for rectlinear, cylindrical or spherical grids)
   sets: downwind face coord - cell center coord - cell volume, 
   overwrites these to minimize memory usage and save time
* @param sizeCell: half thickness or radius (m)
* @return
*/
void particle_1D::setMesh(
    int& numCells, std::vector<double>& sizeCell)
{
    xFacePositive.assign(numCells, 0.0);    //downwind face location
    xCellCenter.assign(numCells, 0.0);    //cell center location
    cellVolume.assign(numCells, 0.0);   //cell volume
	double xL1 = 0.0;            //first upwind face location (always at 0.0)

    xFacePositive[0] = sizeCell[0];
	for (int i = 1; i < numCells; i++)
	{
        xFacePositive[i] = sizeCell[i] + xFacePositive[i-1];
	}

    xCellCenter[0] = 0.5 * (xL1 + xFacePositive[0]);
	for (int i = 1; i < numCells; i++)
	{
        xCellCenter[i] = 0.5 * (xFacePositive[i-1] + xFacePositive[i]);
	}

    shape->set_cellVolumes(numCells, xFacePositive, cellVolume);
}

/**
* Set the particle temperature, composition, etc.
* @param Temp_ particle: temperature (K).
* @param wetSolidVolFraction_: volume fraction of wet solid (-).
* @param drySolidVolFraction_: volume fraction of dry solid (-).
* @param charVolFraction_: volume fraction of char (-).
* @param ashVolFraction_: volume fraction of ash (-).
* @param particleO2MassFraction_: mass fraction of O2 inside the pores (-).
* @param particlePressure_: gauge pressure inside the pores (pa).
* @return
*/
void particle_1D::preStepForward(
                                double localTimeStepSize_,
                                double particleSize_,
                                std::vector<double> Temp_,
                                std::vector<double> wetSolidVolFraction_,
                                std::vector<double> drySolidVolFraction_,
                                std::vector<double> charVolFraction_,
                                std::vector<double> ashVolFraction_,
                                std::vector<double> particleO2MassFraction_,
                                std::vector<double> particlePressure_,
                                std::vector<double> integralWetSolidMass_,
                                std::vector<double> integralDrySolidMass_,
                                std::vector<double> integralCharMass_
)
{
    localTimeStepSize = localTimeStepSize_;

    numCells = Temp_.size();
    Temp = Temp_;
    wetSolidVolFraction = wetSolidVolFraction_;
    drySolidVolFraction = drySolidVolFraction_;
    charVolFraction = charVolFraction_;
    ashVolFraction = ashVolFraction_;
    particleO2MassFraction = particleO2MassFraction_;
    particlePressure = particlePressure_;
    integralWetSolidMass =  integralWetSolidMass_;
    integralDrySolidMass = integralDrySolidMass_;
    integralCharMass = integralCharMass_;

    surfaceTemp = Temp.back();
    coreTemp = Temp.front();

    shape->currentSize = particleSize_;
    cellSize.assign(numCells, shape->currentSize /numCells);
    setMesh(numCells, cellSize);
}

/**
* Set the particle temperature, composition, etc (at time 0)
* @return
*/
void particle_1D::initialize(
                            double Temp_0,
                            double wetSolidVolFraction_0,
                            double drySolidVolFraction_0,
                            double charVolFraction_0,
                            double ashVolFraction_0,
                            double particleO2MassFraction_0,
                            double particlePressure_0
)
{

    state = burning;

    // Set particle size
    initParticleSize = shape->currentSize;
    shape->air = air;

    //set spatial resolution ~ 100-micron or at least 5 cells
    numCells = round(initParticleSize / defaultMeshResolution);
    numCells = max(5, numCells);
    Temp.assign(numCells, Temp_0);
    wetSolidVolFraction.assign(numCells, wetSolidVolFraction_0);
    drySolidVolFraction.assign(numCells, drySolidVolFraction_0);
    charVolFraction.assign(numCells, charVolFraction_0);
    ashVolFraction.assign(numCells, ashVolFraction_0);
    particleO2MassFraction.assign(numCells, particleO2MassFraction_0);
    particlePressure.assign(numCells, particlePressure_0);

    surfaceTemp = Temp.back();
    coreTemp = Temp.front();

    // Set the mesh
    cellSize.assign(numCells, initParticleSize /numCells);
    setMesh(numCells, cellSize);

    // Estimate the local time-step size based on the characteristic timescales
    localTimeStepSize = min(getChemicalTimescale(Temp.back()), getDiffusionTimescale(cellSize.back(), Temp.back()));

    // Calculate initial values of solid masses
    initWetSolidMass.assign(numCells, 0.0);
    initDrySolidMass.assign(numCells, 0.0);
    initCharMass.assign(numCells, 0.0);

    for (int j = 0; j < numCells; j++)
    {
        initWetSolidMass[j] = max(wetSolid->get_bulkDensity(Temp[j]) * wetSolidVolFraction[j] * cellVolume[j], 1e-14);
        initDrySolidMass[j] = max(drySolid->get_bulkDensity(Temp[j]) * drySolidVolFraction[j] * cellVolume[j], 1e-14);
        initCharMass[j] = max(Char->get_bulkDensity(Temp[j]) * charVolFraction[j] * cellVolume[j], 1e-14);
    }

    // Initialize arrays for integration of solid masses
    integralWetSolidMass.assign(numCells, 0.0);
    integralDrySolidMass.assign(numCells, 0.0);
    integralCharMass.assign(numCells, 0.0);
}

/**
* De-allocate memory
* @param
* @return
*/
void particle_1D::destroy()
{

}

/**
* Writes coordinates of cell centroid to a line of the file
* @param outfile: output file
* @return
*/
void particle_1D::writeCoordLine(FILE* outfile)
{
    fprintf(outfile, "%.6f", *pCurrentTime);
    for (int i = 0; i < numCells; i++)
    {
        fprintf(outfile, "\t%.6f", xCellCenter[i]);
    }
    fprintf(outfile, "\n");
}

/**
* Writes particle temperature distribution to a line of the file
* @param outfile: output file
* @return
*/
void particle_1D::writeTempLine(FILE* outfile)
{
	fprintf(outfile, "%.6f", *pCurrentTime);
    for (int i = 0; i < numCells; i++)
    {
        fprintf(outfile, "\t%.6f", Temp[i]);
    }
	fprintf(outfile, "\n");
}

/**
* Writes vol fraction of wet solid to a line of the file
* @param outfile: output file
* @return
*/
void particle_1D::writeWetSolidLine(FILE* outfile)
{
	fprintf(outfile, "%.6f", *pCurrentTime);
    for (int i = 0; i < numCells; i++)
    {
        fprintf(outfile, "\t%.6f", wetSolidVolFraction[i]);
    }
	fprintf(outfile, "\n");
}

/**
* Writes vol fraction of dry solid to a line of the file
* @param outfile: output file
* @return
*/
void particle_1D::writeDrySolidLine(FILE* outfile)
{
	fprintf(outfile, "%.6f", *pCurrentTime);
    for (int i = 0; i < numCells; i++)
    {
        fprintf(outfile, "\t%.6f", drySolidVolFraction[i]);
    }
	fprintf(outfile, "\n");
}

/**
* Writes vol fraction of char to a line of the file
* @param outfile: output file
* @return
*/
void particle_1D::writeCharLine(FILE* outfile)
{
	fprintf(outfile, "%.6f", *pCurrentTime);
    for (int i = 0; i < numCells; i++)
    {
        fprintf(outfile, "\t%.6f", charVolFraction[i]);
    }
	fprintf(outfile, "\n");
}

/**
* Writes vol fraction of ash to a line of the file
* @param outfile: output file
* @return
*/
void particle_1D::writeAshLine(FILE* outfile)
{
    fprintf(outfile, "%.6f", *pCurrentTime);
    for (int i = 0; i < numCells; i++)
    {
        fprintf(outfile, "\t%.6f", ashVolFraction[i]);
    }
    fprintf(outfile, "\n");
}

/**
* Writes mass fraction of oxygen to a line of the file
* @param outfile: output file
* @return
*/
void particle_1D::writeO2Line(FILE* outfile)
{
    fprintf(outfile, "%.6f", *pCurrentTime);
    for (int i = 0; i < numCells; i++)
    {
        fprintf(outfile, "\t%.6f", particleO2MassFraction[i]);
    }
    fprintf(outfile, "\n");
}

/**
* Writes gauge pressure to a line of the file
* @param outfile: output file
* @return
*/
void particle_1D::writePressureLine(FILE* outfile)
{
    fprintf(outfile, "%.6f", *pCurrentTime);
    for (int i = 0; i < numCells; i++)
    {
        fprintf(outfile, "\t%.6f", particlePressure[i]);
    }
    fprintf(outfile, "\n");
}

/**
* Steps forward one time step, computing particle heating/cooling.
* @param timeStepSize: Global time step size (s).
* @param externalTemperature: Gas temperature around particle (K).
* @param externalVelocity: Gas velocity around the particle (m/s).
* @param externalO2MassFrac: Oxygen mass fraction in ambient gas around the particle [-].
* @param externalPressure: Gauge pressure of the gas outside the particle.
* @param externalIrradiation: Irradiation on particle surface (w/m^2).
* @return
*/
void particle_1D::stepForward(
    const double globalTimeStepSize,
    const double externalTemperature,
    const double externalVelocity,
    const double externalO2MassFrac,
    const double externalPressure,
    const double externalIrradiation
)
{
    // Set time
    localTime = 0.0;
    localTimeStepIndex = 0;
    getNumLocalTimeSteps(globalTimeStepSize);

    // Reset the time accumulated output quantities
    globalMassReleased = 0.0;
    globalGasFuelReleased = 0.0;
    globalMoistureReleased = 0.0;
    globalCO2Released = 0.0;
    globalHeatReleased = 0.0;

    // Local time loop
    while ((localTimeStepIndex < finalTimeStepIndex) && (state != burned))
    {
        localTimeStepIndex++;
        localTime = localTime + localTimeStepSize;

        // Save solution at previous local time-step
        numCells_old = numCells;
        cellVolume_old = cellVolume;
        cellSize_old = cellSize;
        Temp_old = Temp;
        wetSolidVolFraction_old = wetSolidVolFraction;
        drySolidVolFraction_old = drySolidVolFraction;
        charVolFraction_old = charVolFraction;
        ashVolFraction_old = ashVolFraction;
        particleO2MassFraction_old = particleO2MassFraction;
        particlePressure_old = particlePressure;
        surfaceTemp_old = surfaceTemp;

        // Calculate the reaction & thermophyscial properties at old time-step
        calcReactionThermo();

        // Accumulate time integration of solid species masses (used by the reaction rates)
        accumulateSolidMass();

        // Calculate the surface areas of cell boundary faces at old time-step
        calcCellFaceArea();

        // Calculate the heat transfer coefficient based on external flow state
        h_conv = shape->get_convectiveHeatTransferCoefficient(externalTemperature, Temp.back(), externalVelocity, shape->currentSize);
        correctForBlowing();

        // Initial guess for the iterative loop
        Temp_newIter = Temp_old;
        wetSolidVolFraction_newIter = wetSolidVolFraction_old;
        drySolidVolFraction_newIter = drySolidVolFraction_old;
        charVolFraction_newIter = charVolFraction_old;
        ashVolFraction_newIter = ashVolFraction_old;
        particleO2MassFraction_newIter = particleO2MassFraction_old;
        particlePressure_newIter = particlePressure_old;
        surfaceTemp_newIter = surfaceTemp_old;

        Temp_oldIter = Temp_old;
        wetSolidVolFraction_oldIter = wetSolidVolFraction_old;
        drySolidVolFraction_oldIter = drySolidVolFraction_old;
        charVolFraction_oldIter = charVolFraction_old;
        ashVolFraction_oldIter = ashVolFraction_old;
        particleO2MassFraction_oldIter = particleO2MassFraction_old;
        particlePressure_oldIter = particlePressure_old;
        surfaceTemp_oldIter = surfaceTemp_old;

        // Begin the iterative loop
        iter = 0;
        iterMaxError = 1.0;
        while ((iter < maxIter) && (iterMaxError > 0.001))
        {
            iter++;

            // Save solution at previous iteration
            for (int i = 0; i < numCells_old; i++)
            { 
                //Safety: limit temperature and O2 mass fraction
                particleO2MassFraction_newIter[i] = max(0.0, min(particleO2MassFraction_newIter[i], 1.0));
                Temp_newIter[i] = max(273.0, min(Temp_newIter[i], 3000.0));

                //Use the average of old and new iterations
                Temp_oldIter[i] = 0.5*Temp_newIter[i] + 0.5*Temp_oldIter[i];
                wetSolidVolFraction_oldIter[i] = 0.5*wetSolidVolFraction_newIter[i] + 0.5*wetSolidVolFraction_oldIter[i];
                drySolidVolFraction_oldIter[i] = 0.5*drySolidVolFraction_newIter[i] + 0.5*drySolidVolFraction_oldIter[i];
                charVolFraction_oldIter[i] = 0.5*charVolFraction_newIter[i] + 0.5*charVolFraction_oldIter[i];
                ashVolFraction_oldIter[i] = 0.5*ashVolFraction_newIter[i] + 0.5*ashVolFraction_oldIter[i];
                particleO2MassFraction_oldIter[i] = 0.5*particleO2MassFraction_newIter[i] + 0.5*particleO2MassFraction_oldIter[i];
                particlePressure_oldIter[i] = 0.5*particlePressure_newIter[i] + 0.5*particlePressure_oldIter[i];
                surfaceTemp_oldIter = 0.5*surfaceTemp_newIter + 0.5*surfaceTemp_oldIter;
            }

            // Calculate reaction rates from information at old iteration
            calcReaction_oldIter();

            // Calculate temperature at new time step, at new iteration
            vectA.assign(numCells_old, 0.0);
            vectB.assign(numCells_old, 0.0);
            vectC.assign(numCells_old, 0.0);
            vectD.assign(numCells_old, 0.0);

            energy_conservation(externalTemperature, externalIrradiation);
            Temp_newIter = solveTriDiag(vectA, vectB, vectC, vectD);

            vectA.clear();
            vectB.clear();
            vectC.clear();
            vectD.clear();

            // Calculate solid species volume fractions and cell volumes at new iteration
            mass_conservation();

            // Calculate O2 mass fraction at new iteration
            vectA.assign(numCells_old, 0.0);
            vectB.assign(numCells_old, 0.0);
            vectC.assign(numCells_old, 0.0);
            vectD.assign(numCells_old, 0.0);

            O2_mass_conservation(externalO2MassFrac);
            particleO2MassFraction_newIter = solveTriDiag(vectA, vectB, vectC, vectD);

            vectA.clear();
            vectB.clear();
            vectC.clear();
            vectD.clear();

            // Calculate pressure at new iteration
            vectA.assign(numCells_old, 0.0);
            vectB.assign(numCells_old, 0.0);
            vectC.assign(numCells_old, 0.0);
            vectD.assign(numCells_old, 0.0);

            pressure_equation(externalPressure);
            particlePressure_newIter = solveTriDiag(vectA, vectB, vectC, vectD);

            vectA.clear();
            vectB.clear();
            vectC.clear();
            vectD.clear();

            // Calculate the numerical error of the iterative loop
            calcIterError();

            iterMaxError = max({
                                    TempError,
                                    wetSolidError, drySolidError, charError, ashError,
                                    O2Error,
                                    pressureError
                                });

        } //end iteration loop

        if (iter >= maxIter)
        {
            cout << "**** Error. Iterative loop did not converge after " << iter << " iterations" << endl;
            cout << "**** Iterative loop maximum error = " << iterMaxError << endl;
            cout.precision(12);
            cout << "--local time (s) = " << localTime << endl;
            cout << "   --local time step size (s) = " << localTimeStepSize << endl;
            cout << "   --local time-step index= " << localTimeStepIndex << endl;
            cout << "   --iterative loop #iterations = " << iter << endl;
            cout << "   --iterative loop maximum error = " << iterMaxError << endl;
            cout << "   --temperature error = " << TempError << endl;
            cout << "   --wetSolid error = " << wetSolidError << endl;
            cout << "   --drySolid error = " << drySolidError << endl;
            cout << "   --char error = " << charError << endl;
            cout << "   --ash error = " << ashError << endl;
            cout << "   --O2 error = " << O2Error << endl;
            cout << "   --pressure error = " << pressureError << endl;
            cout << "   --max temperature = " << *std::max_element(Temp.begin(), Temp.end()) << endl;
            cout << "   --max pressure = " << *std::max_element(particlePressure.begin(), particlePressure.end()) << endl;
            throw exception();
        }

        // Update solution at new time
        Temp = Temp_newIter;
        wetSolidVolFraction = wetSolidVolFraction_newIter;
        drySolidVolFraction = drySolidVolFraction_newIter;
        charVolFraction = charVolFraction_newIter;
        ashVolFraction = ashVolFraction_newIter;
        particleO2MassFraction = particleO2MassFraction_newIter;
        particlePressure = particlePressure_newIter;

        // Update computational grid (in case of volume change)
        moveMesh(xFacePositive, xCellCenter, cellSize, numCells_old, cellVolume);

        // Update the conditions at the exposed surface
        updateExposedSurface(externalTemperature, externalVelocity,
                             externalO2MassFrac, externalIrradiation);

        // Update particle state and check burnout
        state = checkState();

        // Update model outputs
        updateOutputs();

        // Re-mesh the particle due to size change
        remeshing();

        // check prints
        #if defined _DEBUG
        cout.precision(12);
        cout << "--local time (s) = " << localTime << endl;
        cout << "   --local time step size (s) = " << localTimeStepSize << endl;
        cout << "   --local time-step index= " << localTimeStepIndex << endl;
        cout << "   --iterative loop #iterations = " << iter << endl;
        cout << "   --iterative loop maximum error = " << iterMaxError << endl;
        cout << "   --temperature error = " << TempError << endl;
        cout << "   --wetSolid error = " << wetSolidError << endl;
        cout << "   --drySolid error = " << drySolidError << endl;
        cout << "   --char error = " << charError << endl;
        cout << "   --ash error = " << ashError << endl;
        cout << "   --O2 error = " << O2Error << endl;
        cout << "   --pressure error = " << pressureError << endl;
        #endif      
    } //end time loop

    // Time-step restriction
    adjustTimeStep(globalTimeStepSize);
}

/**
* Calculation of reaction rates and thermophysical properties of the porous medium
* @param
* @return
*/
void particle_1D::calcReactionThermo()
{
    // All variables must be resized due to mesh changes
    R1reactionRate.assign(numCells_old, 0.0);
    R2reactionRate.assign(numCells_old, 0.0);
    R3reactionRate.assign(numCells_old, 0.0);
    R4reactionRate.assign(numCells_old, 0.0);

    porosity.assign(numCells_old, 0.0);
    effectiveConductivity.assign(numCells_old, 0.0);
    effectiveVolHeatCapacity.assign(numCells_old, 0.0);
    diffusivity.assign(numCells_old, 0.0);
    permeability.assign(numCells_old, 0.0);

    for (int i = 0; i < numCells_old; i++)
    {
        //calculate reaction rates

        R1reactionRate[i] = pow(wetSolid->get_bulkDensity(Temp_old[i]) * wetSolidVolFraction_old[i] * cellVolume_old[i], R1->get_n())
                            * pow(integralWetSolidMass[i] + initWetSolidMass[i] , 1.0 - R1->get_n())
                            * R1->get_A()
                            * exp(-R1->get_Ta() / Temp_old[i]);

        R2reactionRate[i] = pow(drySolid->get_bulkDensity(Temp_old[i]) * drySolidVolFraction_old[i] * cellVolume_old[i], R2->get_n())
                            * pow(integralDrySolidMass[i] + initDrySolidMass[i] , 1.0 - R2->get_n())
                            * R2->get_A()
                            * exp(-R2->get_Ta() / Temp_old[i]);

        R3reactionRate[i] = pow(drySolid->get_bulkDensity(Temp_old[i]) * drySolidVolFraction_old[i] * cellVolume_old[i], R3->get_n())
                            * pow(integralDrySolidMass[i] + initDrySolidMass[i] , 1.0 - R3->get_n())
                            * pow(particleO2MassFraction_old[i] / 0.226 , R3->get_nO2())
                            * R3->get_A()
                            * exp(-R3->get_Ta() / Temp_old[i]);

        R4reactionRate[i] = pow(Char->get_bulkDensity(Temp_old[i]) * charVolFraction_old[i] * cellVolume_old[i], R4->get_n())
                            * pow(integralCharMass[i] + initCharMass[i] , 1.0 - R4->get_n())
                            * pow(particleO2MassFraction_old[i] / 0.226 , R4->get_nO2())
                            * R4->get_A()
                            * exp(-R4->get_Ta() / Temp_old[i]);

        //porosity
        porosity[i] = wetSolid->get_porosity(Temp_old[i]) * wetSolidVolFraction_old[i]
                    + drySolid->get_porosity(Temp_old[i]) * drySolidVolFraction_old[i]
                    + Char->get_porosity(Temp_old[i]) * charVolFraction_old[i]
                    + ash->get_porosity(Temp_old[i]) * ashVolFraction_old[i];

        //effective thermal conductivity, porous medium treatment:
        //conductivity in solid (including radiation) + porosity x conductivity in gas
        effectiveConductivity[i] = wetSolid->get_conductivity(Temp_old[i]) * wetSolidVolFraction_old[i]
                                 + drySolid->get_conductivity(Temp_old[i]) * drySolidVolFraction_old[i]
                                 + Char->get_conductivity(Temp_old[i]) * charVolFraction_old[i]
                                 + ash->get_conductivity(Temp_old[i]) * ashVolFraction_old[i]
                                 + sigma * pow(Temp_old[i], 3.0) *
                                     ( wetSolid->get_radConductivity(Temp_old[i]) * wetSolidVolFraction_old[i]
                                     + drySolid->get_radConductivity(Temp_old[i]) * drySolidVolFraction_old[i]
                                     + Char->get_radConductivity(Temp_old[i]) * charVolFraction_old[i]
                                     + ash->get_radConductivity(Temp_old[i]) * ashVolFraction_old[i] 
                                     )
                                 +  porosity[i] * air->get_k(Temp_old[i]);

        //effective volumetric heat capacity, porous medium treatment:
        //density x specific heat of solid + porosity x density x specific heat of gas
        effectiveVolHeatCapacity[i] = wetSolid->get_bulkDensity(Temp_old[i]) * wetSolid->get_specificHeat(Temp_old[i]) * wetSolidVolFraction_old[i]
                                 + drySolid->get_bulkDensity(Temp_old[i]) * drySolid->get_specificHeat(Temp_old[i]) * drySolidVolFraction_old[i]
                                 + Char->get_bulkDensity(Temp_old[i]) * Char->get_specificHeat(Temp_old[i]) * charVolFraction_old[i]
                                 + ash->get_bulkDensity(Temp_old[i]) * ash->get_specificHeat(Temp_old[i]) * ashVolFraction_old[i]
                                 + porosity[i] * air->get_rho(Temp_old[i]) * air->get_cSubP(Temp_old[i]);

        //permeability
        permeability[i] = wetSolid->get_permeability(Temp_old[i]) * wetSolidVolFraction_old[i]
                        + drySolid->get_permeability(Temp_old[i]) * drySolidVolFraction_old[i]
                        + Char->get_permeability(Temp_old[i]) * charVolFraction_old[i]
                        + ash->get_permeability(Temp_old[i]) * ashVolFraction_old[i];

        //diffusivity of the gas inside the pores (Schmidt number = 1)
        diffusivity[i] = air->get_v(Temp_old[i]);
    }

    //effective thermal conductivity at the particle surface
    surfaceConductivity = effectiveConductivity.back();

    //emissivity at the particle surface
    surfaceEmissivity = wetSolid->get_emissivity(Temp_old.back()) * wetSolidVolFraction_old.back()
        + drySolid->get_emissivity(Temp_old.back()) * drySolidVolFraction_old.back()
        + Char->get_emissivity(Temp_old.back()) * charVolFraction_old.back()
        + ash->get_emissivity(Temp_old.back()) * ashVolFraction_old.back();

    //grid spaceing at the particle surface
    surfaceGridSpacing = xFacePositive.back() - xCellCenter.back();

}

/**
* Accumulation of solid species masses (used by the reaction rate equation)
* @param
* @return
*/
void particle_1D::accumulateSolidMass()
{
    for (int i = 0; i < numCells_old; i++)
    {
        //accumulate time integral of solid species mass
        integralWetSolidMass[i] = integralWetSolidMass[i];
        integralDrySolidMass[i] = integralDrySolidMass[i] + localTimeStepSize * R1->get_productYield() * R1reactionRate[i];
        integralCharMass[i] = integralCharMass[i] + localTimeStepSize * (R2->get_productYield() * R2reactionRate[i]
                                                                        + R3->get_productYield() * R3reactionRate[i]);
    }
}


/**
* Calculation of reaction rates from informations at the old iteration
* @param
* @return
*/
void particle_1D::calcReaction_oldIter()
{
    for (int i = 0; i < numCells; i++)
    {

        R1reactionRate[i] = pow(wetSolid->get_bulkDensity(Temp_oldIter[i]) * wetSolidVolFraction_oldIter[i]
                            * cellVolume_old[i], R1->get_n())
                            * pow(integralWetSolidMass[i] + initWetSolidMass[i] , 1.0 - R1->get_n())
                            * R1->get_A()
                            * exp(-R1->get_Ta() / Temp_oldIter[i]);

        R2reactionRate[i] = pow(drySolid->get_bulkDensity(Temp_oldIter[i]) * drySolidVolFraction_oldIter[i]
                            * cellVolume_old[i], R2->get_n())
                            * pow(integralDrySolidMass[i] + initDrySolidMass[i] , 1.0 - R2->get_n())
                            * R2->get_A()
                            * exp(-R2->get_Ta() / Temp_oldIter[i]);

        R3reactionRate[i] = pow(drySolid->get_bulkDensity(Temp_oldIter[i]) * drySolidVolFraction_oldIter[i]
                            * cellVolume_old[i], R3->get_n())
                            * pow(integralDrySolidMass[i] + initDrySolidMass[i] , 1.0 - R3->get_n())
                            * pow(particleO2MassFraction_oldIter[i] / 0.226 , R3->get_nO2())
                            * R3->get_A()
                            * exp(-R3->get_Ta() / Temp_oldIter[i]);

        R4reactionRate[i] = pow(Char->get_bulkDensity(Temp_oldIter[i]) * charVolFraction_oldIter[i]
                            * cellVolume_old[i], R4->get_n())
                            * pow(integralCharMass[i] + initCharMass[i] , 1.0 - R4->get_n())
                            * pow(particleO2MassFraction_oldIter[i] / 0.226, R4->get_nO2())
                            * R4->get_A()
                            * exp(-R4->get_Ta() / Temp_old[i]);
    }

}

/**
* Solves mass conservation equation
* @param
* @return
*/
void particle_1D::mass_conservation()
{
    double wetSolidVolume_old = 0.0;
    double wetSolidVolume_oldIter = 0.0;
    double wetSolidVolume_newIter = 0.0;

    double drySolidVolume_old = 0.0;
    double drySolidVolume_oldIter = 0.0;
    double drySolidVolume_newIter = 0.0;

    double charVolume_old = 0.0;
    double charVolume_oldIter = 0.0;
    double charVolume_newIter = 0.0;

    double ashVolume_old = 0.0;
    double ashVolume_oldIter = 0.0;
    double ashVolume_newIter = 0.0;

    double RHSp = 0.0;
    double RHSm = 0.0;

    for (int i = 0; i < numCells_old; i++)
    {
        // Update wet solid volume in the cell
        wetSolidVolume_old = wetSolidVolFraction_old[i] * cellVolume_old[i];
        wetSolidVolume_oldIter = wetSolidVolFraction_oldIter[i] * cellVolume_old[i];

        RHSp = 0.0;
        RHSm = -R1reactionRate[i] / wetSolid->get_bulkDensity(Temp_oldIter[i]);

        wetSolidVolume_newIter = 
                (wetSolidVolume_old + localTimeStepSize * RHSp + solidSpeciesURF * wetSolidVolume_oldIter)
                / (1 - localTimeStepSize * RHSm / max(wetSolidVolume_oldIter, 1e-14) + solidSpeciesURF);


        // Update dry solid volume in the cell
        drySolidVolume_old = drySolidVolFraction_old[i] * cellVolume_old[i];
        drySolidVolume_oldIter = drySolidVolFraction_oldIter[i] * cellVolume_old[i];

        RHSp = R1reactionRate[i] * R1->get_productYield() / drySolid->get_bulkDensity(Temp_oldIter[i]);
        RHSm = -(R2reactionRate[i] + R3reactionRate[i]) / drySolid->get_bulkDensity(Temp_oldIter[i]);

        drySolidVolume_newIter =
            (drySolidVolume_old + localTimeStepSize * RHSp + solidSpeciesURF * drySolidVolume_oldIter)
            / (1 - localTimeStepSize * RHSm / max(drySolidVolume_oldIter, 1e-14) + solidSpeciesURF);


        // Update char volume in the cell
        charVolume_old = charVolFraction_old[i] * cellVolume_old[i];
        charVolume_oldIter = charVolFraction_oldIter[i] * cellVolume_old[i];

        RHSp = R2reactionRate[i] * R2->get_productYield() / Char->get_bulkDensity(Temp_oldIter[i])
             + R3reactionRate[i] * R3->get_productYield() / Char->get_bulkDensity(Temp_oldIter[i]);
        RHSm = -R4reactionRate[i] / Char->get_bulkDensity(Temp_oldIter[i]);

        charVolume_newIter =
            (charVolume_old + localTimeStepSize * RHSp + solidSpeciesURF * charVolume_oldIter)
            / (1 - localTimeStepSize * RHSm / max(charVolume_oldIter, 1e-14) + solidSpeciesURF);


        // Update ash volume in the cell
        ashVolume_old = ashVolFraction_old[i] * cellVolume_old[i];
        ashVolume_oldIter = ashVolFraction_oldIter[i] * cellVolume_old[i];

        RHSp =  R4reactionRate[i] * R4->get_productYield() / ash->get_bulkDensity(Temp_oldIter[i]);
        RHSm = 0.0;

        ashVolume_newIter =
            (ashVolume_old + localTimeStepSize * RHSp + solidSpeciesURF * ashVolume_oldIter)
            / (1 - localTimeStepSize * RHSm / max(ashVolume_oldIter, 1e-14) + solidSpeciesURF);


        // Update the cell volume
        cellVolume[i] = wetSolidVolume_newIter + drySolidVolume_newIter
                      + charVolume_newIter + ashVolume_newIter;


        // Update the volume fractions, Note: Sum(VolFrac) = 1 by construction
        wetSolidVolFraction_newIter[i] = wetSolidVolume_newIter / cellVolume[i];
        drySolidVolFraction_newIter[i] = drySolidVolume_newIter / cellVolume[i];
        charVolFraction_newIter[i] = charVolume_newIter / cellVolume[i];
        ashVolFraction_newIter[i] = ashVolume_newIter / cellVolume[i];

        // Safety: enforce 0 <= VolFraction <= 1
        wetSolidVolFraction_newIter[i] = max(0.0, min(1.0, wetSolidVolFraction_newIter[i]));
        drySolidVolFraction_newIter[i] = max(0.0, min(1.0, drySolidVolFraction_newIter[i]));
        charVolFraction_newIter[i] = max(0.0, min(1.0, charVolFraction_newIter[i]));
        ashVolFraction_newIter[i] = max(0.0, min(1.0, ashVolFraction_newIter[i]));
    }
}

/**
* Sets matrix coefficient of the energy equation
* @param
* @return
*/
void particle_1D::energy_conservation(const double externalTemperature, const double externalIrradiation)
{
    std::vector<double> dRRTemp, bRRTemp;
    dRRTemp.assign(numCells_old, 0.0); //contribute to RHS
    bRRTemp.assign(numCells_old, 0.0); //contribute to diagonal

    double Qdotp = 0.0;   //positive heat of reactions
    double Qdotm = 0.0;   //negative heat of reactions

    for (int i = 0; i < numCells_old; i++)
    {
        Qdotp = (1.0 - R1->get_productYield()) * R1reactionRate[i] * max(0.0, R1->get_deltaH())
              + (1.0 - R2->get_productYield()) * R2reactionRate[i] * max(0.0, R2->get_deltaH())
              + (1.0 - R3->get_productYield()) * R3reactionRate[i] * max(0.0, R3->get_deltaH())
              + (1.0 - R4->get_productYield()) * R4reactionRate[i] * max(0.0, R4->get_deltaH());

        Qdotm = (1.0 - R1->get_productYield()) * R1reactionRate[i] * min(0.0, R1->get_deltaH())
              + (1.0 - R2->get_productYield()) * R2reactionRate[i] * min(0.0, R2->get_deltaH())
              + (1.0 - R3->get_productYield()) * R3reactionRate[i] * min(0.0, R3->get_deltaH())
              + (1.0 - R4->get_productYield()) * R4reactionRate[i] * min(0.0, R4->get_deltaH());

        dRRTemp[i] = Qdotp * localTimeStepSize / effectiveVolHeatCapacity[i] / cellVolume_old[i];
        bRRTemp[i] = -Qdotm * localTimeStepSize / effectiveVolHeatCapacity[i] / Temp_oldIter[i] / cellVolume_old[i];
    }

	// Constructing Temperature matrix (filling vectors A,B,C,D)
	int i;

    // ---- Innermost surface boundary cell
    i = 0;
    FO_L = 0;
    set_heatFO_R(i, localTimeStepSize);
    CFL_L = 0;
    set_heatCFL_R(i, localTimeStepSize);

    vectA[i] = 0;
    vectB[i] = 1 + 0.5 * FO_R - 0.5 * CFL_R + bRRTemp[i];
    vectC[i] = -0.5 * FO_R + 0.5 * CFL_R;
    vectD[i] = (1 - 0.5 * FO_R + 0.5 * CFL_R) * Temp_old[i]
                + (0.5 * FO_R - 0.5 * CFL_R) * Temp_old[i+1] + dRRTemp[i];

    vectB[i] += temperatureURF; //apply under-relaxation
    vectD[i] += temperatureURF*Temp_oldIter[i];

    // ---- Exposed surface boundary cell
    i = numCells_old - 1;
    FO_R = 0;
    set_heatFO_L(i, localTimeStepSize);
    CFL_R = 0;
    set_heatCFL_L(i, localTimeStepSize);

    vectA[i] = -0.5 * FO_L - 0.5 * CFL_L;
    vectB[i] = 1 + 0.5 * FO_L + 0.5 * CFL_L + bRRTemp[i];
    vectC[i] = 0;
    vectD[i] = (0.5 * FO_L + 0.5 * CFL_L) * Temp_old[i-1]
                + (1 - 0.5 * FO_L - 0.5 * CFL_L) * Temp_old[i] + dRRTemp[i];

    h_rad = surfaceEmissivity * sigma * pow(surfaceTemp_old, 3.0);

	Bi = (h_conv + h_rad) * surfaceGridSpacing / surfaceConductivity;

    surfaceTemp_newIter = (Temp_oldIter.back() + (h_conv * externalTemperature + surfaceEmissivity * externalIrradiation) 
                          * surfaceGridSpacing / surfaceConductivity) 
                          / (1 + Bi);

	vectB[i] = vectB[i] + ((h_conv + h_rad) / (1 + Bi))
             * localTimeStepSize / effectiveVolHeatCapacity[i] * areaFacePositive[i]/cellVolume[i];
    
    vectD[i] = vectD[i] + ((h_conv * externalTemperature + surfaceEmissivity * externalIrradiation) / (1 + Bi))
             * localTimeStepSize / effectiveVolHeatCapacity[i] * areaFacePositive[i]/cellVolume[i];

    vectB[i] += temperatureURF; //apply under-relaxation
    vectD[i] += temperatureURF * Temp_oldIter[i];

    // ---- Interior cells
    for (i = 1; i < (numCells_old - 1); i++)
    {
        set_heatFO_L(i, localTimeStepSize);
        set_heatFO_R(i, localTimeStepSize);
        set_heatCFL_L(i, localTimeStepSize);
        set_heatCFL_R(i, localTimeStepSize);

        vectA[i] = -0.5 * FO_L - 0.5 * CFL_L;
        vectB[i] = 1 + 0.5 * FO_L + 0.5 * FO_R + 0.5 * CFL_L - 0.5 * CFL_R + bRRTemp[i];
        vectC[i] = -0.5 * FO_R + 0.5 * CFL_R;
        vectD[i] = (0.5 * FO_L + 0.5 * CFL_L) * Temp_old[i-1]
                    + (1 - 0.5 * FO_L - 0.5 * FO_R - 0.5 * CFL_L + 0.5 * CFL_R) * Temp_old[i]
                    + (0.5 * FO_R - 0.5 * CFL_R) * Temp_old[i+1] + dRRTemp[i];

        vectB[i] += temperatureURF; //apply under-relaxation
        vectD[i] += temperatureURF * Temp_oldIter[i];
    }
}

/**
* Sets matrix coefficient of the oxygen mass conservation equation
* @param
* @return
*/
void particle_1D::O2_mass_conservation(const double externalO2MassFrac)
{
    std::vector<double> bRRO2;
    bRRO2.assign(numCells_old, 0.0); //contribute to diagonal

    for (int i = 0; i < numCells_old; i++)
    {
        bRRO2[i] = (
                      (1.0 - R1->get_productYield()) * R1reactionRate[i] * particleO2MassFraction_oldIter[i]
                    + (1.0 - R2->get_productYield()) * R2reactionRate[i] * particleO2MassFraction_oldIter[i]
                    + (1.0 - R3->get_productYield()) * R3reactionRate[i] * particleO2MassFraction_oldIter[i]
                    + (1.0 - R4->get_productYield()) * R4reactionRate[i] * particleO2MassFraction_oldIter[i]
                    + R3->get_O2Yield() * R3reactionRate[i]
                    + R4->get_O2Yield() * R4reactionRate[i]
                    ) / cellVolume_old[i];

        bRRO2[i] *= localTimeStepSize / air->get_rho(Temp_oldIter[i]) / porosity[i] 
                    / max(particleO2MassFraction_oldIter[i], 1e-14);
    }

    // Constructing Temperature matrix (filling vectors A,B,C,D)
    int i;

    // ---- Innermost surface boundary cell
    i = 0;
    FO_L = 0;
    set_massFO_R(i, localTimeStepSize);
    CFL_L = 0;
    set_massCFL_R(i, localTimeStepSize);

    vectA[i] = 0;
    vectB[i] = 1 + 0.5 * FO_R - 0.5 * CFL_R + bRRO2[i];
    vectC[i] = -0.5 * FO_R + 0.5 * CFL_R;
    vectD[i] = (1 - 0.5 * FO_R + 0.5 * CFL_R) * particleO2MassFraction_old[i]
             + (0.5 * FO_R - 0.5 * CFL_R) * particleO2MassFraction_old[i+1];

    vectB[i] += O2URF; //apply under-relaxation
    vectD[i] += O2URF * particleO2MassFraction_oldIter[i];

    // ---- Exposed surface boundary cell
    i = numCells_old - 1;
    FO_R = 0;
    set_massFO_L(i, localTimeStepSize);
    CFL_R = 0;
    set_massCFL_L(i, localTimeStepSize);

    vectA[i] = -0.5 * FO_L - 0.5 * CFL_L;
    vectB[i] = 1 + 0.5 * FO_L + 0.5 * CFL_L + bRRO2[i];
    vectC[i] = 0;
    vectD[i] = (0.5 * FO_L + 0.5 * CFL_L) * particleO2MassFraction_old[i-1]
             + (1 - 0.5 * FO_L - 0.5 * CFL_L) * particleO2MassFraction_old[i];

    h_mass = h_conv / air->get_cSubP(Temp_oldIter[i]);

    Bi = h_mass * surfaceGridSpacing 
        / (porosity[i] * air->get_rho(Temp_oldIter[i]) * diffusivity[i]);

    vectB[i] = vectB[i] + (h_mass / (1 + Bi))
             * localTimeStepSize / air->get_rho(Temp_oldIter[i]) /porosity[i] * areaFacePositive[i] / cellVolume[i];

    vectD[i] = vectD[i] + ((h_mass * externalO2MassFrac) / (1 + Bi))
             * localTimeStepSize / air->get_rho(Temp_oldIter[i]) / porosity[i] * areaFacePositive[i] / cellVolume[i];

    vectB[i] += O2URF; //apply under-relaxation
    vectD[i] += O2URF * particleO2MassFraction_oldIter[i];

    // ---- Interior cells
    for (i = 1; i < (numCells_old - 1); i++)
    {
        set_massFO_L(i, localTimeStepSize);
        set_massFO_R(i, localTimeStepSize);
        set_massCFL_L(i, localTimeStepSize);
        set_massCFL_R(i, localTimeStepSize);

        vectA[i] = -0.5 * FO_L - 0.5 * CFL_L;
        vectB[i] = 1 + 0.5 * FO_L + 0.5 * FO_R + 0.5 * CFL_L - 0.5 * CFL_R + bRRO2[i];
        vectC[i] = -0.5 * FO_R + 0.5 * CFL_R;
        vectD[i] = (0.5 * FO_L + 0.5 * CFL_L) * particleO2MassFraction_old[i-1]
                    + (1 - 0.5 * FO_L - 0.5 * FO_R - 0.5 * CFL_L + 0.5 * CFL_R) * particleO2MassFraction_old[i]
                    + (0.5 * FO_R - 0.5 * CFL_R) * particleO2MassFraction_old[i+1];

        vectB[i] += O2URF; //apply under-relaxation
        vectD[i] += O2URF * particleO2MassFraction_oldIter[i];
    }
}

/**
* Sets matrix coefficient of the pressure equation
* @param
* @return
*/
void particle_1D::pressure_equation(const double externalPressure)
{
    std::vector<double> dRRpressure;
    dRRpressure.assign(numCells_old, 0.0); //contribute to RHS

    for (int i = 0; i < numCells_old; i++)
    {
        dRRpressure[i] = ((1.0 - R1->get_productYield()) * R1reactionRate[i]
                        + (1.0 - R2->get_productYield()) * R2reactionRate[i]
                        + (1.0 - R3->get_productYield()) * R3reactionRate[i]
                        + (1.0 - R4->get_productYield()) * R4reactionRate[i]
                        );
    }

    // Constructing Temperature matrix (filling vectors A,B,C,D)
    int i;

    // ---- Innermost surface boundary cell
    i = 0;

    vectA[i] = 0;
    vectC[i] = -0.5 * (permeability[i+1] / diffusivity[i+1] + permeability[i] / diffusivity[i])
                * areaFacePositive[i] / (xCellCenter[i+1] - xCellCenter[i]);
    vectB[i] = -vectC[i];
    vectD[i] = dRRpressure[i];

    vectB[i] += pressureURF; //apply under-relaxation
    vectD[i] += pressureURF * particlePressure_oldIter[i];

    // ---- Exposed surface boundary cell
    i = numCells_old - 1;

    vectA[i] = -0.5 * (permeability[i] / diffusivity[i] + permeability[i-1] / diffusivity[i-1])
                * areaFaceNegative[i] / (xCellCenter[i] - xCellCenter[i-1]);
    vectB[i] = -vectA[i] + permeability[i] / diffusivity[i] * areaFacePositive[i] / surfaceGridSpacing;
    vectC[i] = 0.0;
    vectD[i] = permeability[i] / diffusivity[i] * areaFacePositive[i] / surfaceGridSpacing * externalPressure
                + dRRpressure[i];

    vectB[i] += pressureURF; //apply under-relaxation
    vectD[i] += pressureURF * particlePressure_oldIter[i];

    // ---- Interior cells
    for (i = 1; i < (numCells_old - 1); i++)
    {
        vectA[i] = -0.5 * (permeability[i] / diffusivity[i] + permeability[i-1] / diffusivity[i-1])
                    * areaFaceNegative[i] / (xCellCenter[i] - xCellCenter[i-1]);        
        vectC[i] = -0.5 * (permeability[i+1] / diffusivity[i+1] + permeability[i] / diffusivity[i])
                    * areaFacePositive[i] / (xCellCenter[i+1] - xCellCenter[i]);  

        vectB[i] = -vectA[i] - vectC[i];
        vectD[i] = dRRpressure[i];

        vectB[i] += pressureURF; //apply under-relaxation
        vectD[i] += pressureURF * particlePressure_oldIter[i];
    }
}

/**
* Calculate heat transfer's Fourier number (FO), for left face
* @param i: cell index [-]
* @param dt_: time-step size [s]
* @return
*/
void particle_1D::set_heatFO_L(const int &i , const double &dt_)
{
    FO_L = 0.5 * (effectiveConductivity[i] + effectiveConductivity[i-1]) / effectiveVolHeatCapacity[i] 
         * areaFaceNegative[i] * dt_ / (cellVolume[i] * (xCellCenter[i] - xCellCenter[i-1]));
}
/**
* Calculate mass transfer's Fourier number (FO), for left face
* @param i: cell index [-]
* @param dt_: time-step size [s]
* @return
*/
void particle_1D::set_massFO_L(const int& i, const double& dt_)
{
    FO_L = 0.5 * (porosity[i] * diffusivity[i] * air->get_rho(Temp_oldIter[i]) + porosity[i-1] * diffusivity[i-1] * air->get_rho(Temp_oldIter[i-1])) 
         / (porosity[i] * air->get_rho(Temp_oldIter[i])) 
         * areaFaceNegative[i] * dt_ / (cellVolume[i] * (xCellCenter[i] - xCellCenter[i-1]));
}

/**
* Calculate heat transfer's Fourier number (FO), for right face
* @param i: cell index [-]
* @param dt_: time-step size [s]
* @return
*/void particle_1D::set_heatFO_R(const int &i , const double &dt_)
{
    FO_R = 0.5 * (effectiveConductivity[i] + effectiveConductivity[i+1]) / effectiveVolHeatCapacity[i] 
         * areaFacePositive[i] * dt_ / (cellVolume[i] * (xCellCenter[i+1] - xCellCenter[i]));
}
/**
* Calculate mass transfer's Fourier number (FO), for right face
* @param i: cell index [-]
* @param dt_: time-step size [s]
* @return
*/
void particle_1D::set_massFO_R(const int& i, const double& dt_)
{
    FO_R = 0.5 * (porosity[i] * diffusivity[i] * air->get_rho(Temp_oldIter[i]) + porosity[i+1] * diffusivity[i+1] * air->get_rho(Temp_oldIter[i+1]))
         / (porosity[i] * air->get_rho(Temp_oldIter[i])) 
         * areaFacePositive[i] * dt_ / (cellVolume[i] * (xCellCenter[i+1] - xCellCenter[i]));
}

/**
* Calculate heat transfer's Courant–Friedrichs–Lewy (CFL) number, for left face
* @param i: cell index [-]
* @param dt_: time-step size [s]
* @return
*/
void particle_1D::set_heatCFL_L(const int& i, const double& dt_)
{
    // Case for which mass flux > 0, CFL_L > 0
    if ( 
        ((i == numCells_old - 1) && (particlePressure_old[i-1] > particlePressure_old[i])) 
        ||
        ((i != numCells_old - 1) && (particlePressure_old[i-1] > particlePressure_old[i])
          && (particlePressure_old[i] > particlePressure_old[i+1])) 
        )
    {
        CFL_L = dt_ * permeability[i] / diffusivity[i]
                * air->get_cSubP(Temp_old[i]) / effectiveVolHeatCapacity[i]
                * (particlePressure_old[i-1] - particlePressure_old[i])
                / pow(xCellCenter[i] - xCellCenter[i-1], 2.0);
    }
    // Case for which dp/dx ~ 0
    else
    {
        CFL_L = 0.0;
    }
}
/**
* Calculate mass transfer's Courant–Friedrichs–Lewy (CFL) number, for left face
* @param i: cell index [-]
* @param dt_: time-step size [s]
* @return
*/
void particle_1D::set_massCFL_L(const int& i, const double& dt_)
{
    // Case for which mass flux > 0, CFL_L > 0
    if (
        ((i == numCells_old - 1) && (particlePressure_old[i-1] > particlePressure_old[i]))
        ||
        ((i != numCells_old - 1) && (particlePressure_old[i-1] > particlePressure_old[i])
            && (particlePressure_old[i] > particlePressure_old[i+1]))
        )
    {
        CFL_L = dt_ * permeability[i] / diffusivity[i]
                * 1.0 / (air->get_rho(Temp_old[i]) * porosity[i])
                * (particlePressure_old[i-1] - particlePressure_old[i])
                / pow(xCellCenter[i] - xCellCenter[i-1], 2.0);
    }
    // Case for which dp/dx ~ 0
    else
    {
        CFL_L = 0.0;
    }
}

/**
* Calculate heat transfer's Courant–Friedrichs–Lewy (CFL) number, for right face
* @param i: cell index [-]
* @param dt_: time-step size [s]
* @return
*/void particle_1D::set_heatCFL_R(const int& i, const double& dt_)
{
    // Case for which mass flux < 0, CFL_R < 0
    if (
        ((i == 0) && (particlePressure_old[i] < particlePressure_old[i+1]))
        ||
        ((i != 0) && (particlePressure_old[i-1] < particlePressure_old[i])
            && (particlePressure_old[i] < particlePressure_old[i+1]))
        )
    {
        CFL_R = dt_ * permeability[i] / diffusivity[i]
                * air->get_cSubP(Temp_old[i]) / effectiveVolHeatCapacity[i]
                * (particlePressure_old[i] - particlePressure_old[i+1])
                / pow(xCellCenter[i+1] - xCellCenter[i], 2.0);
    }
    // Case for which dp/dx ~ 0
    else
    {
        CFL_R = 0.0;
    }
}
/**
* Calculate mass transfer's Courant–Friedrichs–Lewy (CFL) number, for right face
* @param i: cell index [-]
* @param dt_: time-step size [s]
* @return
*/
void particle_1D::set_massCFL_R(const int& i, const double& dt_)
{
    // Case for which mass flux < 0, CFL_R < 0
    if (
        ((i == 0) && (particlePressure_old[i] < particlePressure_old[i+1]))
        ||
        ((i != 0) && (particlePressure_old[i-1] < particlePressure_old[i])
            && (particlePressure_old[i] < particlePressure_old[i+1]))
        )
    {
        CFL_R = dt_ * permeability[i] / diffusivity[i]
                * 1.0 / (air->get_rho(Temp_old[i]) * porosity[i])
                * (particlePressure_old[i] - particlePressure_old[i+1])
                / pow(xCellCenter[i+1] - xCellCenter[i], 2.0);
    }
    // Case for which dp/dx ~ 0
    else
    {
        CFL_R = 0.0;
    }
}

/**
* Reconstructs the grid(in case of volume change)
* @param xR: array of coordinates of right face [m]
* @param xC: array of coordinates of cell center [m]
* @param deltaX: array of cells sizes [m]
* @param numCells: number of mesh cells [m]
* @param vol: array of cells volumes [m]
* @return
*/
void particle_1D::moveMesh(std::vector<double>& xR, std::vector<double>& xC, std::vector<double>& deltaX,
    const int& numCells, std::vector<double>& vol)
{
    double accumulatedVolume = 0.0;

    accumulatedVolume += vol[0] * shape->correctForShape();
    deltaX[0] = shape->computeSizeFromVolume(accumulatedVolume);
    xR[0] = deltaX[0];
    xC[0] = 0.5 * xR[0];
    for (int i = 1; i < numCells; i++)
    {
        accumulatedVolume += vol[i] * shape->correctForShape();
        xR[i] = shape->computeSizeFromVolume(accumulatedVolume);
        deltaX[i] = xR[i] - xR[i-1];
        xC[i] = 0.5 * (xR[i-1] + xR[i]);
    }

    shape->currentSize = xR.back();
}
 
/**
* Calculation of surface areas of cell boundary faces
* @param
* @return
*/
void particle_1D::calcCellFaceArea()
{
    areaFacePositive.assign(numCells_old, 0.0);
    areaFaceNegative.assign(numCells_old, 0.0);

    areaFacePositive[0] = shape->get_surfaceArea(xFacePositive[0]);
    areaFaceNegative[0] = shape->get_surfaceArea(0.0);
    for (int i = 1; i < numCells_old; i++)
    {
        areaFacePositive[i] = shape->get_surfaceArea(xFacePositive[i]);
        areaFaceNegative[i] = shape->get_surfaceArea(xFacePositive[i-1]);
    }
}

/**
* Calculation of the maximum absolute difference between two vectors
* @param vectX: elements of the first vector
* @param vectY: elements of the second vector
* @return value of the maximum absolute element-wise difference
*/
double particle_1D::maxAbsDifference(std::vector<double> vectX, std::vector<double> vectY)
{
    std::vector<double> difference;
    difference.resize(vectX.size(), 0.0);

    for (unsigned int i = 0; i < vectX.size(); i++)
    {
        difference[i] = abs(vectX[i] - vectY[i]);
    }

    return *std::max_element(difference.begin(), difference.end());
}

/**
* Calculation of the numerical error of the iterative loop of a variable
* @param
* @return
*/
void particle_1D::calcIterError()
{
    TempDifference_iter = maxAbsDifference(Temp_oldIter, Temp_newIter);
    wetSolidDifference_iter = maxAbsDifference(wetSolidVolFraction_oldIter, wetSolidVolFraction_newIter);
    drySolidDifference_iter = maxAbsDifference(drySolidVolFraction_oldIter, drySolidVolFraction_newIter);
    charDifference_iter = maxAbsDifference(charVolFraction_oldIter, charVolFraction_newIter);
    ashDifference_iter = maxAbsDifference(ashVolFraction_oldIter, ashVolFraction_newIter);
    O2Difference_iter = maxAbsDifference(particleO2MassFraction_oldIter, particleO2MassFraction_newIter);
    pressureDifference_iter = maxAbsDifference(particlePressure_oldIter, particlePressure_newIter);

    TempDifference_time = maxAbsDifference(Temp_old, Temp_newIter);
    wetSolidDifference_time = maxAbsDifference(wetSolidVolFraction_old, wetSolidVolFraction_newIter);
    drySolidDifference_time = maxAbsDifference(drySolidVolFraction_old, drySolidVolFraction_newIter);
    charDifference_time = maxAbsDifference(charVolFraction_old, charVolFraction_newIter);
    ashDifference_time = maxAbsDifference(ashVolFraction_old, ashVolFraction_newIter);
    O2Difference_time = maxAbsDifference(particleO2MassFraction_old, particleO2MassFraction_newIter);
    pressureDifference_time = maxAbsDifference(particlePressure_old, particlePressure_newIter);

    TempError = temperatureURF * TempDifference_iter / max(TempDifference_time, 1e-14);
    O2Error = O2URF * O2Difference_iter / max(O2Difference_time, 1e-14);
    
    wetSolidError = wetSolidDifference_iter / solidSpeciesThreshold;
    drySolidError = drySolidDifference_iter / solidSpeciesThreshold;
    charError = charDifference_iter / solidSpeciesThreshold;
    ashError = ashDifference_iter / solidSpeciesThreshold;

    pressureError = pressureDifference_iter / pressureThreshold;
}

/**
* Check the number of local time-steps is equal to the simulation time
* @param globalTimeStepSize: the global time-step size
* @return
*/
void particle_1D::getNumLocalTimeSteps(const double remianingTime)
{
    localTimeStepSize = min(localTimeStepSize, remianingTime);
    finalTimeStepIndex = round(remianingTime / localTimeStepSize);

    if ((finalTimeStepIndex * localTimeStepSize) != remianingTime)
    {
        finalTimeStepIndex = finalTimeStepIndex + 1;
        localTimeStepSize = remianingTime / finalTimeStepIndex;
    }
}

/**
* Update dt such that abs(Q(n + 1) - Q(n)) <= Threshold
* @param globalTimeStepSize: the global time-step size
* @return
*/
void particle_1D::adjustTimeStep(const double remainingTime)
{
    double dt_Temp, dt_wetSolid, dt_drySolid, dt_char, dt_ash, dt_O2, dt_pressure;

    localTimeStepSize_old = localTimeStepSize;

    // Limit the time step dt so that variations in Q are less than a user defined Threshold
    dt_Temp = localTimeStepSize_old * temperatureThreshold / max(TempDifference_time, 1e-14);
    dt_wetSolid = localTimeStepSize_old * solidSpeciesThreshold / max(wetSolidDifference_time, 1e-14);
    dt_drySolid = localTimeStepSize_old * solidSpeciesThreshold / max(drySolidDifference_time, 1e-14);
    dt_char = localTimeStepSize_old * solidSpeciesThreshold / max(charDifference_time, 1e-14);
    dt_ash = localTimeStepSize_old * solidSpeciesThreshold / max(ashDifference_time, 1e-14);
    dt_O2 = localTimeStepSize_old * O2Threshold / max(O2Difference_time, 1e-14);
    dt_pressure = localTimeStepSize_old * pressureThreshold / max(pressureDifference_time, 1e-14);

    // Limit the time step to a maximum of 10 % variations or a user defined threshold
    dt_Temp = max(0.9 * localTimeStepSize_old, min(1.1 * localTimeStepSize_old, dt_Temp));
    dt_wetSolid = max(0.9 * localTimeStepSize_old, min(1.1 * localTimeStepSize_old, dt_wetSolid));
    dt_drySolid = max(0.9 * localTimeStepSize_old, min(1.1 * localTimeStepSize_old, dt_drySolid));
    dt_char = max(0.9 * localTimeStepSize_old, min(1.1 * localTimeStepSize_old, dt_char));
    dt_ash = max(0.9 * localTimeStepSize_old, min(1.1 * localTimeStepSize_old, dt_ash));
    dt_O2 = max(0.9 * localTimeStepSize_old, min(1.1 * localTimeStepSize_old, dt_O2));
    dt_pressure = max(0.9 * localTimeStepSize_old, min(1.1 * localTimeStepSize_old, dt_pressure));

    localTimeStepSize = min({ timeStepThreshold,
                              dt_Temp,
                              dt_wetSolid, dt_drySolid, dt_char, dt_ash,
                              dt_O2,
                              dt_pressure
                            });

}

/**
* Update the conditions at the exposed surface
* @param externalTemperature: Gas temperature around particle (K).
* @param externalVelocity: Gas velocity around the particle (m/s).
* @param externalO2MassFrac: Oxygen mass fraction in ambient gas around the particle [-].
* @param externalIrradiation: Irradiation on particle surface (w/m^2).
* @return
*/
void particle_1D::updateExposedSurface(const double externalTemperature, const double externalVelocity,
                                       const double externalO2MassFrac, const double externalIrradiation)

{
    // Update thermal properties at the exposed surface

    surfacePorosity = wetSolid->get_porosity(Temp.back()) * wetSolidVolFraction.back()
                    + drySolid->get_porosity(Temp.back()) * drySolidVolFraction.back()
                    + Char->get_porosity(Temp.back()) * charVolFraction.back()
                    + ash->get_porosity(Temp.back()) * ashVolFraction.back();

    surfaceConductivity = wetSolid->get_conductivity(Temp.back()) * wetSolidVolFraction.back()
                        + drySolid->get_conductivity(Temp.back()) * drySolidVolFraction.back()
                        + Char->get_conductivity(Temp.back()) * charVolFraction.back()
                        + ash->get_conductivity(Temp.back()) * ashVolFraction.back()
                        + sigma * pow(Temp.back(), 3.0) *
                             (wetSolid->get_radConductivity(Temp.back()) * wetSolidVolFraction.back()
                            + drySolid->get_radConductivity(Temp.back()) * drySolidVolFraction.back()
                            + Char->get_radConductivity(Temp.back()) * charVolFraction.back()
                            + ash->get_radConductivity(Temp.back()) * ashVolFraction.back()
                            )
                        + surfacePorosity * air->get_k(Temp.back());

    surfaceEmissivity = wetSolid->get_emissivity(Temp.back()) * wetSolidVolFraction.back()
                      + drySolid->get_emissivity(Temp.back()) * drySolidVolFraction.back()
                      + Char->get_emissivity(Temp.back()) * charVolFraction.back()
                      + ash->get_emissivity(Temp.back()) * ashVolFraction.back();

    surfacePermeability = wetSolid->get_permeability(Temp.back()) * wetSolidVolFraction.back()
                        + drySolid->get_permeability(Temp.back()) * drySolidVolFraction.back()
                        + Char->get_permeability(Temp.back()) * charVolFraction.back()
                        + ash->get_permeability(Temp.back()) * ashVolFraction.back();

    surfaceDiffusivity = air->get_v(Temp.back());

    surfaceGridSpacing = xFacePositive.back() - xCellCenter.back();


    // Calculate the heat flux and the temperature at the exposed surface

    h_conv = shape->get_convectiveHeatTransferCoefficient(externalTemperature, surfaceTemp_old, externalVelocity, shape->currentSize);
    correctForBlowing();

    h_rad = surfaceEmissivity * sigma * pow(surfaceTemp_newIter, 3.0);
    Bi = (h_conv + h_rad) * surfaceGridSpacing / surfaceConductivity;

    surfaceTemp = (Temp.back() +
                    (h_conv * externalTemperature + surfaceEmissivity * externalIrradiation)
                    * surfaceGridSpacing / surfaceConductivity) / (1.0 + Bi);

    surfaceTemp = max(273.0, min(surfaceTemp, 3000.0));

    surfaceHeatFluxConv = h_conv * (externalTemperature - surfaceTemp);

    surfaceHeatFluxRad = surfaceEmissivity * externalIrradiation
                        - surfaceEmissivity * sigma * pow(surfaceTemp, 4.0);

    surfaceHeatFlux = surfaceHeatFluxConv + surfaceHeatFluxRad;


    // Calculate the oxygen mass flux at the exposed surface

    h_mass = h_conv / air->get_cSubP(Temp.back());
    Bi = h_mass * surfaceGridSpacing
        / (surfacePorosity * air->get_rho(Temp.back()) * surfaceDiffusivity);

    surfaceO2MassFrac = (Bi * externalO2MassFrac + particleO2MassFraction.back())
                        / (1 + Bi);

    surfaceO2MassFrac = max(0.0, min(surfaceO2MassFrac, 1.0));

    surfaceMassFlux = h_mass * (externalO2MassFrac - surfaceO2MassFrac);

    // Drag coefficient 
    dragCoeff = shape->get_dragCoefficient(externalTemperature, externalVelocity, xFacePositive.back());

}

/**
* Correction of heat transfer coeff. due to blowing 
* using a Couette approximation
* @param
* @return
*/
void particle_1D::correctForBlowing()
{

    double outletMassFlux = 0.0;

    if (particlePressure[numCells-2] > particlePressure[numCells-1])
    {
        outletMassFlux = permeability.back() / diffusivity.back()
                          * (particlePressure[numCells-2]  - particlePressure[numCells-1]) 
                          / (xCellCenter[numCells-1] - xCellCenter[numCells-2]);

        h_conv = outletMassFlux * air->get_cSubP(surfaceTemp) 
                / (exp(outletMassFlux * air->get_cSubP(surfaceTemp) / h_conv) - 1.0 );
    }
}

/**
* Update particle state and check burnout
* @param
* @return the enum of the particle state
*/
Particle::eState particle_1D::checkState()
{
    if ((R3->get_A() == 0) && (R4->get_A() == 0)) //no oxidation
    {
        if ( (R2->get_productYield() == 0) && (shape->currentSize / initParticleSize < 0.01) )
        {
            #if defined _DEBUG
            cout << "particle is completely consumed to zero size" << endl;
            #endif

            return burned;
        }
        else if ((R2->get_productYield() != 0) && (charVolFraction.front() > 0.999))
        {
            #if defined _DEBUG
            cout << "particle is completely transformed into char" << endl;
            #endif

            return burned;
        }
        else
        {
            return burning;
        }
    }
    else
    {
        if ( (R4->get_productYield() == 0) && (shape->currentSize / initParticleSize < 0.01) )
        {
            #if defined _DEBUG
            cout << "particle is completely consumed to zero size" << endl;
            #endif

            return burned;
        }
        else if ((R4->get_productYield() != 0) && (ashVolFraction.front() > 0.999))
        {
            #if defined _DEBUG
            cout << "particle is completely transformed into ash" << endl;
            #endif

            return burned;
        }
        else
        {
            return burning;
        }
    }

    return burned;
}

/**
* Update model output quantities for the global time-step
* @param
* @return
*/
void particle_1D::updateOutputs()
{
    //update particle size
    particleSize = shape->currentSize;
    particleVol = shape->get_volume();
    particleSurfToVolRatio = shape->get_surfaceAreaToVolumeRatio();
    projectedAreaRatio = shape->get_projectedAreaRatio();

    //update temperature of the particle core
    coreTemp = Temp.front();

    //calculate cell based quantities
    localMass_cell.assign(numCells, 0.0);
    localMassLossRate_cell.assign(numCells, 0.0);
    localGasFuelReleaseRate_cell.assign(numCells, 0.0);
    localMoistureReleaseRate_cell.assign(numCells, 0.0);
    localCO2ReleaseRate_cell.assign(numCells, 0.0);
    localHeatReleaseRate_cell.assign(numCells, 0.0);

    for (int j = 0; j < numCells; j++)
    {
        localMass_cell[j] = (
                            wetSolid->get_bulkDensity(Temp[j]) * wetSolidVolFraction[j]
                            + drySolid->get_bulkDensity(Temp[j]) * drySolidVolFraction[j]
                            + Char->get_bulkDensity(Temp[j]) * charVolFraction[j]
                            + ash->get_bulkDensity(Temp[j]) * ashVolFraction[j]
                            )
                            * cellVolume[j];

        localMassLossRate_cell[j] = (1.0 - R1->get_productYield()) * R1reactionRate[j]
                                  + (1.0 - R2->get_productYield()) * R2reactionRate[j]
                                  + (1.0 - R3->get_productYield()) * R3reactionRate[j]
                                  + (1.0 - R4->get_productYield()) * R4reactionRate[j];

        localGasFuelReleaseRate_cell[j] = (1.0 - R2->get_productYield()) * R2reactionRate[j]
                                        + (1.0 - R3->get_productYield()) * R3reactionRate[j];

        localMoistureReleaseRate_cell[j] = (1.0 - R1->get_productYield()) * R1reactionRate[j];

        localCO2ReleaseRate_cell[j] = (1.0 + R4->get_O2Yield() - R4->get_productYield()) * R4reactionRate[j];

        localHeatReleaseRate_cell[j] = (1.0 - R1->get_productYield()) * R1reactionRate[j] * R1->get_deltaH()
                                     + (1.0 - R2->get_productYield()) * R2reactionRate[j] * R2->get_deltaH()
                                     + (1.0 - R3->get_productYield()) * R3reactionRate[j] * R3->get_deltaH()
                                     + (1.0 - R4->get_productYield()) * R4reactionRate[j] * R4->get_deltaH();
    }

    //accumulate cell based quantities
    particleMass = accumulate(localMass_cell.begin(), localMass_cell.end(), 0.0f)
                   * shape->correctForShape();

    localMassLossRate = accumulate(localMassLossRate_cell.begin(), localMassLossRate_cell.end(), 0.0f)
                        * shape->correctForShape();

    localGasFuelReleaseRate = accumulate(localGasFuelReleaseRate_cell.begin(), localGasFuelReleaseRate_cell.end(), 0.0f)
                                * shape->correctForShape();

    localMoistureReleaseRate = accumulate(localMoistureReleaseRate_cell.begin(), localMoistureReleaseRate_cell.end(), 0.0f)
                                * shape->correctForShape();

    localCO2ReleaseRate = accumulate(localCO2ReleaseRate_cell.begin(), localCO2ReleaseRate_cell.end(), 0.0f)
                            * shape->correctForShape();

    localHeatReleaseRate = accumulate(localHeatReleaseRate_cell.begin(), localHeatReleaseRate_cell.end(), 0.0f)
                             * shape->correctForShape();

    //accumulate global mass (or heat) produced/consumed in local time [kg or J]
    globalMassReleased     += localMassLossRate * localTimeStepSize;
    globalGasFuelReleased  += localGasFuelReleaseRate * localTimeStepSize;
    globalMoistureReleased += localMoistureReleaseRate * localTimeStepSize;
    globalCO2Released      += localCO2ReleaseRate * localTimeStepSize;
    globalHeatReleased     += localHeatReleaseRate * localTimeStepSize;

    //calculate global mass production/consumption rates in [kg/s]
    globalMassLossRate        = globalMassReleased / localTime;
    globalGasFuelReleaseRate  = globalGasFuelReleased / localTime;
    globalMoistureReleaseRate = globalMoistureReleased / localTime;
    globalCO2ReleaseRate      = globalCO2Released / localTime;

    //calculate global heat release rate in [J/s]
    globalHeatReleaseRate = globalHeatReleased / localTime;

    //calculate global reaction rates in [kg/s] (at the last local time for diagnostic objectives)
    globalR1reactionRate = accumulate(R1reactionRate.begin(), R1reactionRate.end(), 0.0f);
    globalR2reactionRate = accumulate(R2reactionRate.begin(), R2reactionRate.end(), 0.0f);
    globalR3reactionRate = accumulate(R3reactionRate.begin(), R3reactionRate.end(), 0.0f);
    globalR4reactionRate = accumulate(R4reactionRate.begin(), R4reactionRate.end(), 0.0f);

    //TODO: update firelab outputs
    outTime = localTime;
}

/**
* Re-generate the mesh due to particle shrinking or swelling
* @param
* @return
*/
void particle_1D::remeshing()
{
    // Remeshing step 1: update numCells, cellSize

    cellSize_old = cellSize;
    cellVolume_old = cellVolume;

    numCells = round(shape->currentSize / defaultMeshResolution); 	// Update the number of cells
    numCells = min(numCells, numCells_old);   				        // numCells <= numCells_old
    numCells = max(numCells, 5);						            // Keep minimum of 5 grid cells

    cellSize.assign(numCells, shape->currentSize / numCells);


    // Remeshing, step 2: update computational grid

    xCellCenter_old = xCellCenter;

    setMesh(numCells, cellSize);

    particleVol = shape->get_volume();


    // Remeshing, step 3: interpolate solution on new computational grid
    interpolateOnNewMesh();
}

/**
* Interpolate solution on new mesh based on linear interpolation
* @param
* @return
*/
void particle_1D::interpolateOnNewMesh()
{
    // save the results of the old mesh
    Temp_old = Temp;
    wetSolidVolFraction_old = wetSolidVolFraction;
    drySolidVolFraction_old = drySolidVolFraction;
    charVolFraction_old = charVolFraction;
    ashVolFraction_old = ashVolFraction;
    particleO2MassFraction_old = particleO2MassFraction;
    particlePressure_old = particlePressure;

    // save the time-accumulated quantities
    std::vector<double> integralWetSolidMass_old = integralWetSolidMass; 
    std::vector<double> integralDrySolidMass_old = integralDrySolidMass; 
    std::vector<double> integralCharMass_old = integralCharMass;

    std::vector<double> initWetSolidMass_old = initWetSolidMass;
    std::vector<double> initDrySolidMass_old = initDrySolidMass;
    std::vector<double> initCharMass_old = initCharMass;

    // allocate arrays for saving results on new mesh
    Temp.assign(numCells, 0.0);
    wetSolidVolFraction.assign(numCells, 0.0);
    drySolidVolFraction.assign(numCells, 0.0);
    charVolFraction.assign(numCells, 0.0);
    ashVolFraction.assign(numCells, 0.0);
    particleO2MassFraction.assign(numCells, 0.0);
    particlePressure.assign(numCells, 0.0);

    integralWetSolidMass.assign(numCells, 0.0);
    integralDrySolidMass.assign(numCells, 0.0);
    integralCharMass.assign(numCells, 0.0);

    initWetSolidMass.assign(numCells, 0.0);
    initDrySolidMass.assign(numCells, 0.0);
    initCharMass.assign(numCells, 0.0);

    // use linear interpolation to estimate the results on the new mesh
    int ii = -1;
    for (int i = 0; i < numCells; i++)
    {
        ii = -1;
        for (int j = 0; j < (numCells_old - 1); j++)
        {
            if (((xCellCenter_old[j] - xCellCenter[i]) <= 0.0) && ((xCellCenter[i] - xCellCenter_old[j + 1]) < 0.0))
            {
                ii = j;
            }

        }
        if ((ii == -1) && (xCellCenter[i] < xCellCenter_old[0]))
        {
            ii = 0;
        }
        if ((ii == -1) && (xCellCenter_old[numCells_old - 1] <= xCellCenter[i]))
        {
            ii = numCells_old - 2;
        }
        if (ii > -1)
        {
            Temp[i] = Temp_old[ii]
                + (Temp_old[ii + 1] - Temp_old[ii]) * (xCellCenter[i] - xCellCenter_old[ii])
                / (xCellCenter_old[ii + 1] - xCellCenter_old[ii]);

            wetSolidVolFraction[i] = wetSolidVolFraction_old[ii]
                + (wetSolidVolFraction_old[ii + 1] - wetSolidVolFraction_old[ii]) * (xCellCenter[i] - xCellCenter_old[ii])
                / (xCellCenter_old[ii + 1] - xCellCenter_old[ii]);

            drySolidVolFraction[i] = drySolidVolFraction_old[ii]
                + (drySolidVolFraction_old[ii + 1] - drySolidVolFraction_old[ii]) * (xCellCenter[i] - xCellCenter_old[ii])
                / (xCellCenter_old[ii + 1] - xCellCenter_old[ii]);

            charVolFraction[i] = charVolFraction_old[ii]
                + (charVolFraction_old[ii + 1] - charVolFraction_old[ii]) * (xCellCenter[i] - xCellCenter_old[ii])
                / (xCellCenter_old[ii + 1] - xCellCenter_old[ii]);

            ashVolFraction[i] = ashVolFraction_old[ii]
                + (ashVolFraction_old[ii + 1] - ashVolFraction_old[ii]) * (xCellCenter[i] - xCellCenter_old[ii])
                / (xCellCenter_old[ii + 1] - xCellCenter_old[ii]);

            particleO2MassFraction[i] = particleO2MassFraction_old[ii]
                + (particleO2MassFraction_old[ii + 1] - particleO2MassFraction_old[ii]) * (xCellCenter[i] - xCellCenter_old[ii])
                / (xCellCenter_old[ii + 1] - xCellCenter_old[ii]);

            particlePressure[i] = particlePressure_old[ii]
                + (particlePressure_old[ii + 1] - particlePressure_old[ii]) * (xCellCenter[i] - xCellCenter_old[ii])
                / (xCellCenter_old[ii + 1] - xCellCenter_old[ii]);

            integralWetSolidMass[i] = integralWetSolidMass_old[ii]
                + (integralWetSolidMass_old[ii + 1] - integralWetSolidMass_old[ii]) * (xCellCenter[i] - xCellCenter_old[ii])
                / (xCellCenter_old[ii + 1] - xCellCenter_old[ii]);

            integralDrySolidMass[i] = integralDrySolidMass_old[ii]
                + (integralDrySolidMass_old[ii + 1] - integralDrySolidMass_old[ii]) * (xCellCenter[i] - xCellCenter_old[ii])
                / (xCellCenter_old[ii + 1] - xCellCenter_old[ii]);

            integralCharMass[i] = integralCharMass_old[ii]
                + (integralCharMass_old[ii + 1] - integralCharMass_old[ii]) * (xCellCenter[i] - xCellCenter_old[ii])
                / (xCellCenter_old[ii + 1] - xCellCenter_old[ii]);

            initWetSolidMass[i] = initWetSolidMass_old[ii]
                + (initWetSolidMass_old[ii + 1] - initWetSolidMass_old[ii]) * (xCellCenter[i] - xCellCenter_old[ii])
                / (xCellCenter_old[ii + 1] - xCellCenter_old[ii]);

            initDrySolidMass[i] = initDrySolidMass_old[ii]
                + (initDrySolidMass_old[ii + 1] - initDrySolidMass_old[ii]) * (xCellCenter[i] - xCellCenter_old[ii])
                / (xCellCenter_old[ii + 1] - xCellCenter_old[ii]);

            initCharMass[i] = initCharMass_old[ii]
                + (initCharMass_old[ii + 1] - initCharMass_old[ii]) * (xCellCenter[i] - xCellCenter_old[ii])
                / (xCellCenter_old[ii + 1] - xCellCenter_old[ii]);


            // Safety: enforce 0 <= VolFraction (massFration) <= 1
            wetSolidVolFraction[i] = max(0.0, min(1.0, wetSolidVolFraction[i]));
            drySolidVolFraction[i] = max(0.0, min(1.0, drySolidVolFraction[i]));
            charVolFraction[i] = max(0.0, min(1.0, charVolFraction[i]));
            ashVolFraction[i] = max(0.0, min(1.0, ashVolFraction[i]));  
            particleO2MassFraction[i] = max(0.0, min(1.0, particleO2MassFraction[i])); 

            // Safety: disallow negative mass
            integralWetSolidMass[i] = max(1e-14, integralWetSolidMass[i]);
            integralDrySolidMass[i] = max(1e-14, integralDrySolidMass[i]);
            integralCharMass[i] = max(1e-14, integralCharMass[i]);

        }
        if (ii == -1)
        {
            cout << "*** Error in linear interpolation routine ***" << endl;
            break;
            throw exception();
        }
    }

}