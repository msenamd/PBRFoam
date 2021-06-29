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
    Foam::biomassCloud

Description
    A Cloud of solid particles

SourceFiles
    biomassCloudI.H
    biomassCloud.C
    biomassCloudIO.C

\*---------------------------------------------------------------------------*/

#ifndef biomassCloud_H
#define biomassCloud_H


#include "Cloud.H"
#include "biomassParticle.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "fvMesh.H"
#include "volFields.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class biomassCloud Declaration
\*---------------------------------------------------------------------------*/

class biomassCloud
:
    public Cloud<biomassParticle>
{
    // Private data

        const fvMesh& mesh_;

        IOdictionary particleProperties_;

        // particle properties
        scalar nParticles_;           //number of particles in a super particle
        
        word fuel_;                   //composition of released fuel
        word firebrands_;             //motion flag
        scalar e_;                    //particle elasticity
        scalar mu_;                   //particle viscosity

        word dragModel_;              //drag model
        scalar dragCoeff_;            //particle drag coefficient

        word geometry_;               //particle geometry
        scalar resolution_;            //particle grid resolution
		scalar cylinderLength_;       //length of cylindeical particles
		scalar rectangleLength_;      //length of rectangular particles
        scalar rectangleWidth_;       //width of rectangular particles

		scalar eta_c_;                //particle char yield

		scalar rho_m_;                //moisture density
		scalar rho_vs_;               //solid density
		scalar rho_c_;                //char density  

		scalar c_m_;                  //moisture heat capacity
		scalar c_vs_;              	  //solid heat capacity
		scalar c_c_;                  //char heat capacity 

		scalar k_m_;                  //moisture heat conductivity
		scalar k_vs_;              	  //solid heat conductivity
		scalar k_c_;                  //char heat conductivity 

		scalar e_m_;                  //moisture emssivity
		scalar e_vs_;              	  //solid heat emssivity
		scalar e_c_;                  //char heat emssivity 


        // References to the carrier gas fields
            //- Density [kg/m^3]
            const volScalarField& rho_g_;

            //- Velocity [m/s]
            const volVectorField& U_g_;

            //- Dynamic viscosity [Pa.s]
            const volScalarField& mu_g_;

            //- Temperature[K]
            const volScalarField& T_g_;

            //- Radiation Flux [W/m2]
            volScalarField G_g_;

            //- Gravity
            const dimensionedVector& g_;

            //- Thermo
            const SLGThermo& thermo_;

        // Packing Ratio
            autoPtr<volScalarField> packingRatio_;

        // Surface temperature (K)
            autoPtr<volScalarField> surfaceTemp_;            

        // Surface to Volume Ratio (1/m)
            autoPtr<volScalarField> surfaceToVolumeRatio_;

        // Momentum (kg m/s)
            autoPtr<volVectorField> momentum_;

        // Sources

            //- Mass Transfer - Total (kg/m3)
            autoPtr<volScalarField::Internal> rhoTrans_;

            //- Mass Transfer - species
            PtrList<volScalarField::Internal> rhoYTrans_;

            //- Momentum Transfer (force/m3)
            autoPtr<volVectorField::Internal> UTrans_;

            //- Energy Transfer (J/m3)
            autoPtr<volScalarField::Internal> Qr3Trans_;
            autoPtr<volScalarField::Internal> QconvTrans_;

            //- Radiation temperature^4
            autoPtr<volScalarField::Internal> absorption_;
            autoPtr<volScalarField::Internal> emissionTrans_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        biomassCloud(const biomassCloud&);

        //- Disallow default bitwise assignment
        void operator=(const biomassCloud&);


public:

    // Constructors

        //- Construct given mesh
        biomassCloud
        (
            const fvMesh&,
            const word& cloudName,
            const volScalarField& rho_g,
            const volVectorField& U_g,
            const dimensionedVector& g,
            const SLGThermo& thermo,
            bool readFields = true
        );


    // Member Functions

        // Access

            bool hasWallImpactDistance() const;

            inline const fvMesh& mesh() const;

            inline scalar nParticles() const;

            inline word fuel() const;

            inline word firebrands() const;

            inline scalar e() const;
            inline scalar mu() const;

            inline word dragModel() const;
            inline scalar dragCoeff() const;

            inline word geometry() const;
            inline scalar resolution() const;
            inline scalar cylinderLength() const;
            inline scalar rectangleLength() const;
            inline scalar rectangleWidth() const;

            inline scalar eta_c() const;
            inline scalar rho_m() const;
            inline scalar rho_vs() const;
            inline scalar rho_c() const;
            inline scalar c_m() const;
            inline scalar c_vs() const;
            inline scalar c_c() const;
            inline scalar k_m() const;
            inline scalar k_vs() const;
            inline scalar k_c() const;
            inline scalar e_m() const;
            inline scalar e_vs() const;
            inline scalar e_c() const;
            

            // Thermo class
            inline const SLGThermo& thermo() const;

            // Packing Ratio
            inline volScalarField::Internal& packingRatio();
            inline const volScalarField::Internal& packingRatio() const;

            // Surface Temperatre
            inline volScalarField::Internal& surfaceTemp();
            inline const volScalarField::Internal& surfaceTemp() const;

            // Surface to Volume Ratio
            inline volScalarField::Internal& surfaceToVolumeRatio();
            inline const volScalarField::Internal& surfaceToVolumeRatio() const;

            // momentum
            inline volVectorField::Internal& momentum();
            inline const volVectorField::Internal& momentum() const;


            // Sources

            	// Mass

                    //- Return const access to mass source fields
                    inline const volScalarField::Internal&
                        rhoTrans() const;

                    //- Return reference to mass source fields
                    inline volScalarField::Internal&
                        rhoTrans();

                    //- Return total mass source term [kg/s/m3]
                    inline tmp<volScalarField::Internal> Srho() const;

                // Species Mass

                    //- Return reference to mass source for field i
                    inline volScalarField::Internal&
                        rhoYTrans(const label i);

                    //- Return const access to mass source fields
                    inline const PtrList<volScalarField::Internal>&
                        rhoYTrans() const;

                    //- Return reference to mass source fields
                    inline PtrList<volScalarField::Internal>&
                        rhoYTrans();

                    //- Return mass source term for specie i - specie eqn                        
                    inline tmp<volScalarField::Internal> SYi(const label i) const;                 

                // Momentum

                    //- Return reference to momentum source
                    inline volVectorField::Internal& 
                        UTrans();

                    //- Return const reference to momentum source
                    inline const volVectorField::Internal&
                        UTrans() const;

                    //- Return momentum source term [force/m3]
                    inline tmp<volVectorField::Internal> SU() const;


                // Energy

                    //- Char oxidation (R3) heat transfer [W/m3.s]
                    inline const volScalarField::Internal& Qr3Trans() const;
                    inline volScalarField::Internal& Qr3Trans();

                    //- Return Char oxidation (R3) heat transfer [W/m3.s]
                    inline tmp<volScalarField::Internal> Qr3() const; 

                    //- Particle-Gas convective heat transfer [J/m3]
                    inline const volScalarField::Internal& QconvTrans() const;
                    inline volScalarField::Internal& QconvTrans();

                    //- Return Char oxidation (R3) heat transfer [W/m3]
                    inline tmp<volScalarField::Internal> Qconv() const; 


                // Radiation

                    inline const volScalarField::Internal& absorption() const;
                    inline volScalarField::Internal& absorption();
                    inline const volScalarField::Internal& emissionTrans() const;
                    inline volScalarField::Internal& emissionTrans();

                    //- Return emission term
                    inline tmp<volScalarField> Ep() const;

                    //- Return absorption coefficient
                    inline tmp<volScalarField> ap() const;


        // Cloud Evolution Functions

            //- Move the particles, solve thermo-chemical decomposition
            void evolve();

            //- Reset the cloud source terms
            void resetSourceTerms();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "biomassCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
