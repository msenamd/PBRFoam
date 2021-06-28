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
    Simple solid spherical particle class with one-way coupling with the
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

    // Private data

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;

        //- Radius or half thickness
        scalar delta_;

        //- Velocity of parcel
        vector U_;

        // particle state
        bool status_;

        //- Temperautre distribution
        scalarField T_;

        //- Chemical Composition
        scalarField x_m_;
        scalarField x_vs_;

        //output properties:
        //- Surface Temperature
        scalar Tsurf_;

        //- Total volumetric production rate (gaseous fuel+moisture)
        scalar volProdRate_tot_;

        //- Volumetric production rate of gaseous fuel
        scalar volProdRate_fuel_;

        //- Volumetric production rate of H2O
        scalar volProdRate_H2O_;

        //- Voumetric heat release rate (from char oxidation)
        scalar volHRR_;

        // - Surface to Volume Ratio
        scalar surfToVol_;

        // - Particle mass
        scalar mass_;

        // - Particle Volume
        scalar vol_;

        // - Particle convection heat trans. coeff
        scalar hConv_;

        // - Particles drag coeff
        scalar CD_;        
public:

    friend class Cloud<biomassParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::TrackingData<biomassCloud>
    {
        // Interpolators for continuous phase fields

            const interpolationCellPoint<scalar>& rhoInterp_;
            const interpolationCellPoint<vector>& UInterp_;
            const interpolationCellPoint<scalar>& muInterp_;
            const interpolationCellPoint<scalar>& TInterp_;
            const interpolationCellPoint<scalar>& GInterp_;

        //- Local gravitational or other body-force acceleration
        const vector& g_;


    public:

        // Constructors

            inline trackingData
            (
                biomassCloud& spc,
                const interpolationCellPoint<scalar>& rhoInterp,
                const interpolationCellPoint<vector>& UInterp,
                const interpolationCellPoint<scalar>& muInterp,
                const interpolationCellPoint<scalar>& TInterp,
                const interpolationCellPoint<scalar>& GInterp,
                const vector& g
            );


        // Member functions

            inline const interpolationCellPoint<scalar>& rhoInterp() const;
            inline const interpolationCellPoint<vector>& UInterp() const;
            inline const interpolationCellPoint<scalar>& muInterp() const;
            inline const interpolationCellPoint<scalar>& TInterp() const;
            inline const interpolationCellPoint<scalar>& GInterp() const;

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
            const label tetPti,
            const scalar delta,
            const vector& U,
            const bool status,
            scalarField T,
            scalarField x_m,
            scalarField x_vs
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

        // Access

            //- Return radius/half width
            inline scalar delta() const;

            //- Return velocity
            inline const vector& U() const;

            //- Return particle state ("true" or "false")
            inline bool status() const;

            //- Return Temperature
            inline scalarField T() const;

            //- Return composition
            inline scalarField x_m() const;
            inline scalarField x_vs() const;

            //- Return outputs
        	inline scalar Tsurf() const;
        	inline scalar volProdRate_tot() const;
			inline scalar volProdRate_fuel() const;
			inline scalar volProdRate_H2O() const;
			inline scalar volHRR() const;                       


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

        //- Update parcel properties over the time interval
        void calcAll
        (
            trackingData& td,
            const scalar dt_,
            const label celli
        );

        //- Calculate new particle velocity and forces
        vector calcVelocity
        (
            trackingData& td,
            const word  movePartice,
            const scalar dt_,          // timestep
            const scalar mass_p,       // particle mass
            const scalar delta_p,      // particle radius or half thickness
            const scalar rhoc,         // gas density
            const scalar muc,          // gas viscosity
            const vector Uc,           // gas velocity
            vector& dUTrans            // momentum transfer to gas phase
        );


        //- Calculate new particle temperature,size and compisition
        void calcDegradation
        (
            trackingData& td, 
            const scalar dt_,
            const scalar Tc,            //local gas temperature
            const vector Uc,            //local gas velocity
            const scalar Gc,            //local irradation
            const scalar YO2c           //local oxygen mass fraction
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
