#ifndef PARTICLE_LUMPED_H
#define PARTICLE_LUMPED_H

#include "particle.h"
using namespace std;

class particle_lumped: public virtual Particle
{
public:
        particle_lumped();
        virtual ~particle_lumped();
        virtual particle_lumped* clone() const;
        particle_lumped(const particle_lumped& rhs);
        particle_lumped& operator=(const particle_lumped& rhs);

        void initialize(double initialTemp_, double initialPercentMoistureContent_, double initialPercentCharContent_);
        void destroy();

        //write data at global time
        void writeTempLine(FILE* outfile);
        void writeMoistureLine(FILE* outfile);
        void writeSolidLine(FILE* outfile);
	    void writeCharLine(FILE* outfile);
        void writeCoordLine(FILE* outfile);

        //step forward in global time (solve particle for local number of time steps)
        void stepForward(const double timeStepSize, const double gasTemperature, const double gasVelocity,const double irradiation, const double oxygenVolFraction);


protected:

private:

        double localTime;               // for local time loop of the particle (always starts from 0)
        int currentLocalTimeStep;              // local time iteration counter
        int numLocalTimeSteps;                        // total number of particle local time-steps
        double initialSolidVolFraction;                  // initial volume fraction of solid
        double volume;                      // particle volume
        double initialVolume;                    // particle initial volume
        double Temp;                    // particle temperature
        double moistureVolFraction, solidVolFraction, charVolFraction;          // volume fractions of moisture, virgin solid and char
        double rho_p;                   // particle mass density
        double h_conv;                  // convective heat transfer coefficient
        double q_surf;                  // Net surface heat flux[W / m2]
        double eps_surf;                // Particle surface emissivity
        double MLRpuv;                  // Mass Loss Rate per unit volume [kg/s/m3]

        double volume_old;
        double radius_old;
        double Temp_old;
        double rho_p_old;
        double moistureVolFraction_old;
        double solidVolFraction_old;
        double charVolFraction_old;
        double q_surf_old;

        //solve mass conservation
		void mass_conservation_TTC(double x_O2_g);

		//solve energy conservation
		void energy_conservation_TTC(double x_O2_g);
};

#endif // PARTICLE_LUMPED_H
