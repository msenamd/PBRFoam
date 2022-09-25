#ifndef PARTICLE_1D_H
#define PARTICLE_1D_H

#include "particle.h"
#include "TDMA.h"
#include <numeric>      // std::accumulate

using namespace std;

class particle_1D: public virtual Particle
{
public:

 // Constructors
        particle_1D();
        virtual particle_1D* clone() const;
        particle_1D(const particle_1D& rhs);
        particle_1D& operator=(const particle_1D& rhs);

// Default destructor
        virtual ~particle_1D();

// Public member functions
        void initialize(
                            double Temp_0,
                            double wetSolidVolFraction_0,
                            double drySolidVolFraction_0,
                            double charVolFraction_0,
                            double ashVolFraction_0,
                            double particleO2MassFraction_0,
                            double particlePressure_0
                        );

        void preStepForward(
                            eState state_,
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
                        );

        void stepForward(
                            const double globalTimeStepSize,
                            const double externalTemperature,
                            const double externalVelocity,
                            const double externalO2MassFrac,
                            const double externalPressure,
                            const double externalIrradiation
                        );

        void destroy();

        //set convergence critera and relaxation parameters
        void setSolutionControl(
            const int maxIter_,
            const double defaultMeshresolution_,
            const double timeStepThreshold_,
            const double temperatureThreshold_,
            const double solidSpeciesThreshold_,
            const double O2Threshold_,
            const double pressureThreshold_,
            const double temperatureURF_,
            const double solidSpeciesURF_,
            const double O2URF_,
            const double pressureURF_,
            const bool flagRemeshing_
        );

        //write data at global time
        void writeCoordLine(ofstream& outfile);
        void writeTempLine(ofstream& outfile);
        void writeWetSolidLine(ofstream& outfile);
        void writeDrySolidLine(ofstream& outfile);
	    void writeCharLine(ofstream& outfile);
        void writeAshLine(ofstream& outfile);
        void writeO2Line(ofstream& outfile);
        void writePressureLine(ofstream& outfile);

// Public data
        int numCells;                                // number of cells [-]
        std::vector<double> Temp;                    // particle temperature [K]
        std::vector<double> wetSolidVolFraction;     // volume fractions of wet solid [-]
        std::vector<double> drySolidVolFraction;     // volume fractions of dry solid [-]
        std::vector<double> charVolFraction;         // volume fractions of char [-]
        std::vector<double> ashVolFraction;          // volume fractions of ash [-]
        std::vector<double> particleO2MassFraction;  // O2 mass fraction inside the particle pores [-]
        std::vector<double> particlePressure;        // Gauge pressure inside the particle pores [pa]
        std::vector<double> integralWetSolidMass;    // Time integration of wet solid mass production [kg]
        std::vector<double> integralDrySolidMass;    // Time integration of dry solid mass production [kg]
        std::vector<double> integralCharMass;        // Time integration of char mass production [kg]

protected:

private:

        int localTimeStepIndex;         // local time loop index [-]
        int finalTimeStepIndex;         // total number of particle local time-steps [-]
        double initParticleSize;        // initial particle size before time loop [m]
		double surfaceConductivity;     // surface thermal conductivity [W/m/K]
        double surfacePorosity;         // surface porosity [-]
        double surfacePermeability;     // surface permeability [m2]
        double surfaceDiffusivity;      // surface mass diffusivity [m2/s2]
        double surfaceGridSpacing;      // grid spacing at the surface [m]
		double localTime;               // for local time loop of the particle (always starts from 0) [s]
        double FO_L;                    // Fourier number (FO), for left face [-]
        double FO_R;                    // Fourier number (FO), for right face [-]
        double CFL_L;                   // CourantM-^VFriedrichsM-^VLewy number (CFL), for left face [-]
        double CFL_R;                   // CourantM-^VFriedrichsM-^VLewy number (CFL), for right face [-]
        double h_rad;                   // radiative heat transfer coefficient [W/m2/K]
        double Bi;                      // Biot number [-]

        int iter;                       // iterative loop index [-]
        double iterMaxError;            // numerical error from all variables [-]
        int maxIter;                    // Max. number of iterations [-]
        double defaultMeshResolution;   // Default cell size [m]
        double timeStepThreshold;       // Max. allowed time-step size [s]
        double temperatureThreshold;    // Max. value of temperature variation during dt [K]
        double solidSpeciesThreshold;   // Max.value of solid species vol. fraction variation during dt [-]
        double O2Threshold;             // Max.value of gaseous variation during dt[-]
        double pressureThreshold;       // Max.value of pres variation during dt[Pa]
        double temperatureURF;          // Under relaxation factor for temperature [-]
        double solidSpeciesURF;         // Under relaxation factor for solid species [-]
        double O2URF;                   // Under relaxation factor for O2 [-]
        double pressureURF;             // Under relaxation factor for pressure [-]
        bool flagRemeshing;             // Flag for optional remeshing

        double TempDifference_iter;     // temperature difference of old and new iterations [K]
        double wetSolidDifference_iter; // wetSolid vol. frac. difference of old and new iterations [-]
        double drySolidDifference_iter; // drySolid vol. frac. difference of old and new iterations [-]
        double charDifference_iter;     // char vol. frac. difference of old and new iterations [-]
        double ashDifference_iter;      // ash vol. frac. difference of old and new iterations [-]
        double O2Difference_iter;       // oxygen mass frac. difference of old and new iterations [-]
        double pressureDifference_iter; // gauge pressure difference of old and new iterations [Pa]

        double TempDifference_time;     // temperature difference of last iteration and prevous time [K]
        double wetSolidDifference_time; // wetSolid vol. frac. difference of last iteration and prevous time [-]
        double drySolidDifference_time; // drySolid vol. frac. difference of last iteration and prevous times [-]
        double charDifference_time;     // char vol. frac. difference of last iteration and prevous time [-]
        double ashDifference_time;      // ash vol. frac. difference of last iteration and prevous time [-]
        double O2Difference_time;       // oxygen mass frac. difference last iteration and prevous time [-]
        double pressureDifference_time; // gauge pressure difference of last iteration and prevous time [Pa]

        double TempError;               // numerical error of temperature equation [-]
        double wetSolidError;           // numerical error of wetSolid mass conservation equation [-]
        double drySolidError;           // numerical error of drySolid mass conservation equation [-]
        double charError;               // numerical error of char mass conservation equation [-]
        double ashError;                // numerical error of ash mass conservation equation [-]
        double O2Error;                 // numerical error of oxygen mass conservation equation [-]
        double pressureError;           // numerical error of pressure equation [-]

        std::vector<double> initWetSolidMass;        // Initial wet solid mass [kg]
        std::vector<double> initDrySolidMass;        // Initial dry solid mass [kg]
        std::vector<double> initCharMass;            // Initial char mass [kg]

        std::vector<double> cellSize;                // cell size [m] (length of rectlinear grid or radius of shperical and cylindrical shells)
        std::vector<double> cellVolume;              // cell volume [m3]
        std::vector<double> xFacePositive;           // location of downwind face (considered positive face) [m]
        std::vector<double> xCellCenter;             // location of cell center [m]
        std::vector<double> areaFacePositive;        // downwind face area [m2]
        std::vector<double> areaFaceNegative;        // upwind face area [m2]

        std::vector<double> R1reactionRate;          // drying raction rate [kg/s]
        std::vector<double> R2reactionRate;          // thermal pyrolysis raction rate [kg/s]
        std::vector<double> R3reactionRate;          // oxidative pyrolysis raction rate [kg/s]
        std::vector<double> R4reactionRate;          // char oxidation raction rate [kg/s]
        
        std::vector<double> porosity;                // Porosity of the particle [-]
        std::vector<double> effectiveConductivity;   // Effective thermal conductivity of porous particle [W/m/K]
        std::vector<double> effectiveVolHeatCapacity;// Effective volumetric heat capacity porous particle [kg/m3 * J/kg/K]
        std::vector<double> diffusivity;             // Diffusivity of the porous particle [m2/s]
        std::vector<double> permeability;            // Permeability of the porous particle [m2]

		//variables at the old local timestep
        int numCells_old;
        double surfaceTemp_old;
        std::vector<double> xCellCenter_old;
        std::vector<double> cellVolume_old;
        std::vector<double> cellSize_old;
		std::vector<double> Temp_old;
        std::vector<double> wetSolidVolFraction_old;
        std::vector<double> drySolidVolFraction_old;
        std::vector<double> charVolFraction_old;
        std::vector<double> ashVolFraction_old;
        std::vector<double> particleO2MassFraction_old;
        std::vector<double> particlePressure_old;

        //variables at the old iteration
        double surfaceTemp_oldIter;
        std::vector<double> Temp_oldIter;
        std::vector<double> wetSolidVolFraction_oldIter;
        std::vector<double> drySolidVolFraction_oldIter;
        std::vector<double> charVolFraction_oldIter;
        std::vector<double> ashVolFraction_oldIter;
        std::vector<double> particleO2MassFraction_oldIter;
        std::vector<double> particlePressure_oldIter;

        //variables at the new iteration
        double surfaceTemp_newIter;
        std::vector<double> Temp_newIter;
        std::vector<double> wetSolidVolFraction_newIter;
        std::vector<double> drySolidVolFraction_newIter;
        std::vector<double> charVolFraction_newIter;
        std::vector<double> ashVolFraction_newIter;
        std::vector<double> particleO2MassFraction_newIter;
        std::vector<double> particlePressure_newIter;

		//dummy vectors for TDMA conversion
		std::vector<double> vectA;
		std::vector<double> vectB;
		std::vector<double> vectC;
		std::vector<double> vectD;

        //Quantities for calculating model output
        std::vector<double> localMass_cell;                 // Total mass of each cell during local time-step [kg]
        std::vector<double> localWetSolidMass_cell;         // Mass of wet solid in each cell during local time-step [kg]
        std::vector<double> localDrySolidMass_cell;         // Mass of dry solid in each cell during local time-step [kg]
        std::vector<double> localCharMass_cell;             // Mass of char in each cell during local time-step [kg]

        std::vector<double> localMassLossRate_cell;          // Mass Loss Rate of each cell during local time-step [kg/s]
        std::vector<double> localGasFuelReleaseRate_cell;    // Gaseous fuel release rate from each cell during local time-step [kg/s]     
        std::vector<double> localMoistureReleaseRate_cell;   // Moisturel release rate from each cell during local time-step [kg/s]      
        std::vector<double> localCO2ReleaseRate_cell;        // CO2 release rate from each cell during local time-step [kg/s]        
        std::vector<double> localHeatReleaseRate_cell;       // Heat release rate from each cell during local time-step [J/s] 
        double localMassLossRate;                       // Mass Loss Rate during local time-step [kg/s]
        double localGasFuelReleaseRate;                 // Gaseous fuel release rate during local time-step [kg/s]   
        double localMoistureReleaseRate;                // Moisture release rate during local time-step [kg/s]        
        double localCO2ReleaseRate;                     // CO2 release rate during local time-step [kg/s] 
        double localHeatReleaseRate;                    // Heat release rate during local time-step [J/s] 

		// Set computational grid (works for rectlinear, cylindrical or spherical grids)
        void setMesh(int& numCells, std::vector<double>& sizeCell);

		// Reconstructs the grid (in case of volume change)
        void moveMesh(std::vector<double>& xR, std::vector<double>& xC, std::vector<double>& deltaX,
                      const int& numCells, std::vector<double>& vol);

        // Calculation of thermophysical properties of the porous medium
        void calcThermo();

        // Accumulation of solid species masses (used by the reaction rate equation)
        void accumulateSolidMass();

        // Calculation of reaction rates from given temperature and composition ant certain time (or iter)
        void calcReaction(const std::vector<double>& Temp_ ,
                          const std::vector<double>& wetSolidVolFraction_,
                          const std::vector<double>& drySolidVolFraction_,
                          const std::vector<double>& charVolFraction_,
                          const std::vector<double>& ashVolFraction_,
                          const std::vector<double>& particleO2MassFraction_);

		// Solid Species mass conservation
		void mass_conservation();

		// Energy conservation
		void energy_conservation(const double gasTemperature, const double irradiation);

        // Oxygen mass conservation inside the particle
        void O2_mass_conservation(const double externalO2MassFrac);

        // Pressure equation
        void pressure_equation(const double externalPressure);
        
        // Calculate Fourier number (FO), for left face
        void set_heatFO_L(const int &i , const double &dt_);
        void set_massFO_L(const int& i, const double& dt_);

        // Calculate Fourier number (FO), for right face
        void set_heatFO_R(const int &i , const double &dt_);
        void set_massFO_R(const int& i, const double& dt_);

        // Calculate CourantM-^VFriedrichsM-^VLewy number (CFL), for left face
        void set_heatCFL_L(const int& i, const double& dt_);
        void set_massCFL_L(const int& i, const double& dt_);

        // Calculate CourantM-^VFriedrichsM-^VLewy number (CFL), for right face
        void set_heatCFL_R(const int& i, const double& dt_);
        void set_massCFL_R(const int& i, const double& dt_);

        // Calculation of surface areas of cell boundary faces
        void calcCellFaceArea();
  
        // Calculation of the numerical error of the iterative loop
        void calcIterError();

        // Auxilary function to return maximum absolute difference of two vectors
        double maxAbsDifference(std::vector<double> vectX, std::vector<double> vectY);

        // Check the summation of local time-steps is equal to the global time step size
        void getNumLocalTimeSteps(const double remainingTime);

        // Adjust the local time-step of particle integration
		void adjustTimeStep(const double remianingTime);

        // Update the conditions at the exposed surface
        void updateExposedSurface(const double externalTemperature, const double externalVelocity,
                                  const double externalO2MassFrac, const double externalIrradiation);

        // Blowing effect at the exposed surface
        void correctForBlowing();

        // Update particle state and check burnout
        Particle::eState checkState();

        // Update model output quantities for the global time-step
        void updateOutputs();

		// Remeshing: Reduce number of cells while particle is burning
	    void remeshing();
        void interpolateOnNewMesh();
};

#endif // PARTICLE_1D_H
