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
#include <tuple>
#include <omp.h>

#include "solidMaterial.h"
#include "solidReaction.h"
#include "air.h"
#include "TDMA.h"
#include "particle.h"
#include "particle_lumped.h"
#include "particle_1D.h"

using namespace std;

/* NOTE:
This main program serves as an example of the interface between the particle model
and larger frameworks. The parameters coded here could be obtained from external
user input files or calculated and saved within OpenFOAM or Firelab codes
*/

int main()
{
    double startTimer;
    double endTimer;
    startTimer = omp_get_wtime();

    // Conditions of the ambient external gas (external input) 
    double simulationTime = 600;                // Total global simulation time [s]
    double timeStepSize = 1.0;                  // To mimic spread model step size [s]
    double externalGasTemp = 300.0;             // External gas temperature [K]
    double externalGasVel = 1.0;                // External gas velocity [m/s]
    double externalO2MassFrac = 0.233;          // External oxygen mass fraction [-]
    double externalGasPressure = 0.0;           // Gauge prssure of the external gas [pa]
    double externalIrradation = 40.0e3;         // Incident radiation flux [W/m2]

    // Particle solution control (external input)
    int maxIter = 1000;                         // Max. number of iterations [-]
    double meshResolution = 100e-6;             // Spatial resolution of the 1-D partice [m]
    double timeStepThreshold = 0.1;             // Max. allowed time-step size [s]
    double temperatureThreshold = 0.1;          // Max. value of temperature variation during dt [K]
    double solidSpeciesThreshold = 0.01;        // Max.value of solid species vol. fraction variation during dt [-]
    double O2Threshold = 0.01;                  // Max.value of gaseous variation during dt[-]
    double pressureThreshold = 0.1;             // Max.value of pres variation during dt[Pa]
    double temperatureURF = 0.1;                // Under relaxation factor for temperature [-]
    double solidSpeciesURF = 0.0;               // Under relaxation factor for solid species [-]
    double O2URF = 0.1;                         // Under relaxation factor for O2 [-]
    double pressureURF = 0.0;                   // Under relaxation factor for pressure [-]

    // Particle geometric properties (external input)
    std::string shapeName = "slab";          // Particle shape: slab/cylinder/sphere
    double radius = 3.8e-2;             // Half thickness or radius of the particle [m]
    double length = 1.0e-0;             // Length of the particle [m]
    double width = 1.0e-0;              // Width of the particle [m]

    // Particle physical properties (external input)
    //double moistureContent = 0.05;       // Moisture content [-]
    double wetSolidDensity = 380;        // Measured density of the wet solid [kg/m3]
    double drySolidDensity = 361;// wetSolidDensity / (1.0 + moistureContent); // Density of the dry solid [kg/m3]
    double charDensity = 73;             // Mass density of char [kg/m3]
    double ashDensity = 5.7;             // Mass density of ash [kg/m3]

    // User-defined conductivity, heat capacity, emissivity, porosity & permeabiltiy  
    // (external input)
    double k0_ws = 0.186;  // Conductivity of wet solid(at 300 K) [W/m/K]
    double k0_ds = 0.176;  // Conductivity of dry solid(at 300 K) [W/m/K]
    double k0_c = 0.065;   // Conductivity of char(at 300 K) [W/m/K]
    double k0_a = 0.058;   // Conductivity of ash(at 300 K) [W/m/K]
    double nk_ws = 0.185;  // Temperature exponent of conductivity of wet solid
    double nk_ds = 0.594;  // Temperature exponent of conductivity of dry solid
    double nk_c = 0.435;   // Temperature exponent of conductivity of char
    double nk_a = 0.353;   // Temperature exponent of conductivity of ash

    double gamma_ws = 0;    // Effective radiation conductivity of wet solid [m]
    double gamma_ds = 0;    // Effective radiation conductivity of dry solid [m]
    double gamma_c = 3.3e-3;    // Effective radiation conductivity of char [m]
    double gamma_a = 6.4e-3;    // Effective radiation conductivity of ash [m]

    double c0_ws = 1764;    // Heat capacity of wet solid [J/kg/K]
    double c0_ds = 1664;    // Heat capacity of dry solid [J/kg/K]
    double c0_c = 1219;     // Heat capacity of char [J/kg/K]
    double c0_a = 1244;     // Heat capacity of ash [J/kg/K]
    double nc_ws = 0.406;   // Temperature exponent of heat capacity of wet solid
    double nc_ds = 0.660;   // Temperature exponent of heat capacity of dry solid
    double nc_c = 0.283;    // Temperature exponent of heat capacity of char
    double nc_a = 0.315;    // Temperature exponent of heat capacity of ash

    double eps_ws = 0.757;  // Surface emissivity of wet solid [-]
    double eps_ds = 0.759;  // Surface emissivity of dry solid [-]
    double eps_c = 0.957;   // Surface emissivity of char [-]
    double eps_a = 0.955;   // Surface emissivity of ash [-]

    double psi_ws = 0.0256; // Porosity of the particle when pure wet solid [-]
    double psi_ds = 0.0744; // Porosity of the particle when pure solid [-]
    double psi_c = 0.8128;  // Porosity of the particle when pure char [-]
    double psi_a = 0.9854;  // Porosity of the particle when pure ash [-]

    double Kperm_ws = 1e-10;    // Permeability of the particle when pure wet solid [m2]
    double Kperm_ds = 1e-10;    // Permeability of the particle when pure solid [m2]
    double Kperm_c = 1e-10;     // Permeability of the particle when pure char [m2]
    double Kperm_a = 1e-10;     // Permeability of the particle when pure ash [m2]


    // Pass the user input material properties to a material class
    solidMaterial wetSolid, drySolid, Char, ash;
    Air air;

    wetSolid.set_bulkDensity(wetSolidDensity);
    drySolid.set_bulkDensity(drySolidDensity);
    Char.set_bulkDensity(charDensity);
    ash.set_bulkDensity(ashDensity);

    wetSolid.set_conductivity(k0_ws, nk_ws);
    drySolid.set_conductivity(k0_ds, nk_ds);
    Char.set_conductivity(k0_c, nk_c);
    ash.set_conductivity(k0_a, nk_a);

    wetSolid.set_radConductivity(gamma_ws);
    drySolid.set_radConductivity(gamma_ds);
    Char.set_radConductivity(gamma_c);
    ash.set_radConductivity(gamma_a);

    wetSolid.set_specificHeat(c0_ws, nc_ws);
    drySolid.set_specificHeat(c0_ds, nc_ds);
    Char.set_specificHeat(c0_c, nc_c);
    ash.set_specificHeat(c0_a, nc_a);

    wetSolid.set_emissivity(eps_ws);
    drySolid.set_emissivity(eps_ds);
    Char.set_emissivity(eps_c);
    ash.set_emissivity(eps_a);

    wetSolid.set_porosity(psi_ws);
    drySolid.set_porosity(psi_ds);
    Char.set_porosity(psi_c);
    ash.set_porosity(psi_a);

    wetSolid.set_permeability(Kperm_ws);
    drySolid.set_permeability(Kperm_ds);
    Char.set_permeability(Kperm_c);
    ash.set_permeability(Kperm_a);

    // Reaction model
    // The reaction yields are external inputs
    // The values selected here result in no shrinking or swelling
    // use R3 (or R4).set_reaction("passive reaction, 0.0); if the reaction is not used
    double drySolidYield = drySolidDensity / wetSolidDensity;
    double thermalCharYield = charDensity / drySolidDensity;
    double oxidativeCharYield = charDensity / drySolidDensity;
    double ashYield = ashDensity / charDensity;

    solidReaction R1, R2, R3, R4;
    R1.set_reaction("drying", drySolidYield);
    R2.set_reaction("thermal pyrolysis", thermalCharYield);
    R3.set_reaction("oxidative pyrolysis", oxidativeCharYield);
    R4.set_reaction("char oxidation", ashYield);

    // Current state of the particle or initial condition at t = 0
    //(Could be saved in the interface between the model and OpenFOAM/Firelab such that
    // these arrays are not uniform for a resumed simulation)
    std::vector<double> temperature(1, externalGasTemp);         // Initial temperature [K]
    std::vector<double> wetSolidVolFraction(1, 1.0);             // Initial volume fraction of wet solid [-]
    std::vector<double> drySolidVolFraction(1, 0.0);             // Initial volume fraction of dry solid [-]
    std::vector<double> charVolFraction(1, 0.0);                 // Initial volume fraction of char [-]
    std::vector<double> ashVolFraction(1, 0.0);                  // Initial volume fraction of ash [-]
    std::vector<double> O2MassFraction(1, externalO2MassFrac);   // Initial mass fraction of oxygen [-]
    std::vector<double> pressure(1, externalGasPressure);        // Initial gauge pressure [Pa]

    // Write output file headers
    remove("particleTemporalOut.csv");
    std::string header1 = "time (s),delta (m),Tsurf (k),Tcore (k),MLR(kg/s),qsurf(W/m2),msurf(kg/s/m2),h_conv(W/m2/k),YO2surf(-),HRR(W)";
    ofstream writeTemporal("particleTemporalOut.csv", ios::app);
    writeTemporal << header1 << endl;

    remove("temperature.dat");
    FILE* TempFile;
    TempFile = fopen("temperature.dat", "w+");

    remove("wetSolid.dat");
    FILE* wetSolidFile;
    wetSolidFile = fopen("wetSolid.dat", "w+");

    remove("drySolid.dat");
    FILE* drySolidFile;
    drySolidFile = fopen("drySolid.dat", "w+");

    remove("char.dat");
    FILE* charFile;
    charFile = fopen("char.dat", "w+");

    remove("ash.dat");
    FILE* ashFile;
    ashFile = fopen("ash.dat", "w+");

    remove("O2.dat");
    FILE* O2File;
    O2File = fopen("O2.dat", "w+");

    remove("pressure.dat");
    FILE* pressureFile;
    pressureFile = fopen("pressure.dat", "w+");

    remove("cellCenter.dat");
    FILE* coordFile;
    coordFile = fopen("cellCenter.dat", "w+");


    // Test run for a certain global time 
    //-----------------------------------------------------------

    particle_1D testParticle;

    testParticle.air = &air;
    testParticle.wetSolid = &wetSolid;
    testParticle.drySolid = &drySolid;
    testParticle.Char = &Char;
    testParticle.ash = &ash;

    testParticle.R1 = &R1;
    testParticle.R2 = &R2;
    testParticle.R3 = &R3;
    testParticle.R4 = &R4;

    testParticle.setGeometry(shapeName, radius, length, width);
  
    testParticle.initialize(
                        temperature, 
                        wetSolidVolFraction, 
                        drySolidVolFraction, 
                        charVolFraction, 
                        ashVolFraction,
                        O2MassFraction,
                        pressure
                    );

    testParticle.setSolutionControl(
                        maxIter,
                        meshResolution,
                        timeStepThreshold,
                        temperatureThreshold,
                        solidSpeciesThreshold,
                        O2Threshold,
                        pressureThreshold,
                        temperatureURF,
                        solidSpeciesURF,
                        O2URF,
                        pressureURF
                    );

    double globalTime = 0.0;

    while ((globalTime < simulationTime) && (testParticle.state==Particle::burning))
    {
        globalTime += timeStepSize;
        testParticle.pCurrentTime = &globalTime;

        cout << "----------------------------------------------------" << endl;
        cout << "----------------------------------------------------" << endl;
        cout << "Global time (s) = " << globalTime << endl;

        testParticle.stepForward(
                                    timeStepSize,
                                    externalGasTemp,
                                    externalGasVel,
                                    externalO2MassFrac,
                                    externalGasPressure,
                                    externalIrradation
                                );

        writeTemporal << globalTime << "," << testParticle.particleSize << "," << 
                testParticle.surfaceTemp << "," << testParticle.coreTemp << "," << 
                testParticle.globalMassLossRate << "," << testParticle.surfaceHeatFlux << "," << 
                testParticle.surfaceMassFlux << "," << testParticle.h_conv << "," << 
                testParticle.surfaceO2MassFrac << "," << testParticle.globalHeatReleaseRate << endl;

        testParticle.writeTempLine(TempFile);
        testParticle.writeWetSolidLine(wetSolidFile);
        testParticle.writeDrySolidLine(drySolidFile);
        testParticle.writeCharLine(charFile);
        testParticle.writeAshLine(ashFile);
        testParticle.writeO2Line(O2File);
        testParticle.writePressureLine(pressureFile);
        testParticle.writeCoordLine(coordFile);
    }

    writeTemporal.close();
    fclose(TempFile);
    fclose(wetSolidFile);
    fclose(drySolidFile);
    fclose(charFile);
    fclose(ashFile);
    fclose(pressureFile);
    fclose(coordFile);

    endTimer = omp_get_wtime();
    printf("CPU time was %f seconds\n", endTimer - startTimer);

}


