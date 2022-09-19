/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::biomassParticle

Description
    Solid biomass particle class with two-way coupling with the
    continuous phase.

SourceFiles
    biomassParticleI.H
    biomassParticle.C
    biomassParticleIO.C

\*---------------------------------------------------------------------------*/

#ifndef biomassParticle_H
#define biomassParticle_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "interpolationCellPoint.H"
#include "interpolationCell.H"
#include "contiguous.H"
#include "SLGThermo.H"

#include <vector>
#include "particle_1D.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace std;

namespace Foam
{

class biomassCloud;


// Forward declaration of friend functions and operators

class biomassParticle;

Ostream& operator<<(Ostream&, const biomassParticle&);


/*---------------------------------------------------------------------------*\
                           Class biomassParticle Declaration
\*---------------------------------------------------------------------------*/

class biomassParticle
:
    public particle
{

private:

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;

        // Stored particle data

            // Particle ID for diagnostics
            label particleID_;

            // Particle state 
            label particleState_;

            // Number of particles per super particle (i.e., multiples of the calculated one) 
            label nParticlesPerSuperParticle_;

            // Particle time step size
            scalar particledt_;

            // Radius or half thickness [m]
            scalar particleSize_;

            // Velocity of particle [m/s]
            vector particleVelo_;

public:
            //- UMD particle model declaration
            particle_1D p1D;

            // Initial conditions at time Zero
            scalar particleSize0_;
            scalar particleTemp0_;                        
            scalar particlePressure0_;
            scalar particleO2MassFraction0_;
            scalar wetSolidVolFraction0_;
            scalar drySolidVolFraction0_;
            scalar charVolFraction0_;
            scalar ashVolFraction0_;

private:
            // Temperautre distribution [K]
            scalarField particleTemp_;

            // Pressure distribution [pa]
            scalarField particlePressure_;

            // Mass fraction of O2 distribution [-]
            scalarField particleO2MassFraction_;

            // Chemical composition [-]
            scalarField wetSolidVolFraction_;
            scalarField drySolidVolFraction_;
            scalarField charVolFraction_;
            scalarField ashVolFraction_;

            // Integral mass [kg]
            scalarField integralWetSolidMass_;
            scalarField integralDrySolidMass_;
            scalarField integralCharMass_;

            // Particle mass [kg]
            scalar particleMass_;

            // Particle volume [m3]
            scalar particleVol_;

            // Surface temperature [K]
            scalar surfaceTemp_;

            // Core temperature [K]
            scalar coreTemp_;

            // Convective heat flux [W/m2]
            scalar convFlux_;

            // Radiative heat flux [W/m2]
            scalar radFlux_;

            // Convective heat flux [kg/m2]
            scalar massFlux_;

            // Surface O2 mass fraction -
            scalar surfaceO2MassFrac_;

            // Particle convection heat trans. coeff [W/m2/K]
            scalar hConv_;

            // Particle drag coeff [-]
            scalar CD_;

            // drying reaction rate [kg/s]
            scalar dryingRate_;        

            // pyrolysis reaction rate [kg/s]
            scalar pyrolysisRate_;  

            // oxidative pyrolysis reaction rate [kg/s]
            scalar oxidPyrolysisRate_;  

            // char oxidation reaction rate [kg/s]
            scalar charOxidRate_;  

            // total mass loss rate [kg/s]
            scalar massLossRate_;  

            // total mass loss rate [W]
            scalar heatReleaseRate_; 

            // Auxilary vectors for data transfer
            std::vector<double> particleTemp_std; 
            std::vector<double> wetSolidVolFraction_std;
            std::vector<double> drySolidVolFraction_std;
            std::vector<double> charVolFraction_std;
            std::vector<double> ashVolFraction_std;
            std::vector<double> particleO2MassFraction_std;
            std::vector<double> particlePressure_std;
            std::vector<double> integralWetSolidMass_std;
            std::vector<double> integralDrySolidMass_std;
            std::vector<double> integralCharMass_std;

            // Output path
            word outputPath;
public:

    friend class Cloud<biomassParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::TrackingData<biomassCloud>
    {
        // Interpolators for continuous phase fields

            const interpolationCell<scalar>& rhoInterp_;
            const interpolationCell<scalar>& YO2Interp_;
            const interpolationCell<scalar>& GInterp_;
            const interpolationCellPoint<vector>& UInterp_;
            const interpolationCell<scalar>& TInterp_;
            const interpolationCell<scalar>& pInterp_;

        //- Local gravitational or other body-force acceleration
        const vector& g_;

    public:

        // Constructors

            inline trackingData
            (
                biomassCloud& bC,
                const interpolationCell<scalar>& rhoInterp,
                const interpolationCell<scalar>& YO2Interp,
                const interpolationCell<scalar>& GInterp,
                const interpolationCellPoint<vector>& UInterp,
                const interpolationCell<scalar>& TInterp,
                const interpolationCell<scalar>& pInterp,
                const vector& g
            );


        // Member functions

            inline const interpolationCell<scalar>& rhoInterp() const;
            inline const interpolationCell<scalar>& YO2Interp() const;
            inline const interpolationCell<scalar>& GInterp() const;
            inline const interpolationCellPoint<vector>& UInterp() const;
            inline const interpolationCell<scalar>& TInterp() const;
            inline const interpolationCell<scalar>& pInterp() const;

            inline const vector& g() const;
    };



    // Constructors

        //- Construct from components
        inline biomassParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            const label tetFacei,
            const label tetPtI,
            const label particleID,
            const label particleState,
            const label nParticlesPerSuperParticle,            
            const scalar particledt,
            const scalar particleSize,
            const vector& particleVelo,
            const scalarField particleTemp,
            const scalarField particlePressure,
            const scalarField particleO2MassFraction,
            const scalarField wetSolidVolFraction,
            const scalarField drySolidVolFraction,
            const scalarField charVolFraction,
            const scalarField ashVolFraction,
            const scalarField integralWetSolidMass,
            const scalarField integralDrySolidMass,
            const scalarField integralCharMass
        );

        //- Construct from Istream
        biomassParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new biomassParticle(*this));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<biomassParticle> operator()(Istream& is) const
            {
                return autoPtr<biomassParticle>
                (
                    new biomassParticle(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Tracking

            //- Move contains all tracking functions
            bool move(trackingData&, const scalar);


        // Patch interactions

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            bool hitPatch
            (
                const polyPatch&,
                trackingData& td,
                const label patchi,
                const scalar trackFraction,
                const tetIndices& tetIs
            );

            //- Overridable function to handle the particle hitting a
            //  processorPatch
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                trackingData& td
            );

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch
            (
                const wallPolyPatch&,
                trackingData& td,
                const tetIndices&
            );

            //- Overridable function to handle the particle hitting a polyPatch
            void hitPatch
            (
                const polyPatch&,
                trackingData& td
            );

            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
            virtual void transformProperties(const tensor& T);

            //- Transform the physical properties of the particle
            //  according to the given separation vector
            virtual void transformProperties(const vector& separation);

            //- The nearest distance to a wall that
            //  the particle can be in the n direction
            virtual scalar wallImpactDistance(const vector& n) const;

    // I-O

        static void readFields(Cloud<biomassParticle>& c);

        static void writeFields(const Cloud<biomassParticle>& c);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const biomassParticle&);


private:

    // Private Member Functions

        //- Update particle variables
        void updateParticle
        (
            trackingData& td,
            const scalar dt_,
            const label celli
        );
};


template<>
inline bool contiguous<biomassParticle>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "biomassParticleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
