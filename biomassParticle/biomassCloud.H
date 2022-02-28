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
    A Cloud of solid biomass particles

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

        //- Reference to external gas density [kg/m^3]
        const volScalarField& rho_;

        //- Reference to external oxygen mass fraction [-]
        const volScalarField& YO2_;

        //- Reference to external radiation Flux [W/m2]
        const volScalarField& G_;

        //- Reference to external gas velocity [m/s]
        const volVectorField& U_;

        //- Reference to Thermo properties
        const SLGThermo& thermo_;

        //- Reference to external gas temperature [K]
        const volScalarField& T_;

        //- Reference to external gas pressure [Pa]
        const volScalarField& p_;

        //- Reference to gravity [m/s2]
        const dimensionedVector& g_;

    // Sources

        //- Mass Transfer - Total [kg/m3]
        autoPtr<DimensionedField<scalar, volMesh>> massTrans_;

        //- Mass Transfer - one per carrier phase specie [kg/m3]
        PtrList<DimensionedField<scalar, volMesh>> speciesMassTrans_;

        //- Momentum Transfer [kg m/s /m3 = N.s/m3]
        autoPtr<DimensionedField<vector, volMesh>> momentumTrans_;

        //- Convective Energy Transfer [J/m3]
        autoPtr<DimensionedField<scalar, volMesh>> QconvTrans_;

        //- RTE Absorption Coefficient [1/m]
        autoPtr<DimensionedField<scalar, volMesh>> absorptionCoeff_;

        //- RTE Emission Transfer [J/m3]
        autoPtr<DimensionedField<scalar, volMesh>> emissionTrans_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        biomassCloud(const biomassCloud&);

        //- Disallow default bitwise assignment
        void operator=(const biomassCloud&);

public:

    // Public data

        // Reactions and Material classes
            particle_1D testParticle;
            solidMaterial wetSolid, drySolid, Char, ash;
            Air air;
            solidReaction R1, R2, R3, R4;

        // Vegetation bed Properties
            scalar nParticles;           //number of particles in a super particle      
            word fuel;                   //composition of released fuel

        // Motion properties
            bool firebrands;             //motion flag
            scalar particleElasticity;   //particle elasticity [Pa]
            scalar particleViscosity;    //particle viscosity [kg/m.s]
            word dragModel;              //drag model
            scalar dragCoeff;            //particle drag coefficient [-]

        // Particle solution control
            label maxIter;                       // Max. number of iterations [-]
            scalar meshResolution;               // Spatial resolution of the 1-D partice [m]
            scalar timeStepThreshold;            // Max. allowed time-step size [s]
            scalar temperatureThreshold;         // Max. value of temperature variation during dt [K]
            scalar solidSpeciesThreshold;        // Max.value of solid species vol. fraction variation during dt [-]
            scalar O2Threshold;                  // Max.value of gaseous variation during dt[-]
            scalar pressureThreshold;            // Max.value of pres variation during dt[Pa]
            scalar temperatureURF;               // Under relaxation factor for temperature [-]
            scalar solidSpeciesURF;              // Under relaxation factor for solid species [-]
            scalar O2URF;                        // Under relaxation factor for O2 [-]
            scalar pressureURF;                  // Under relaxation factor for pressure [-]

        // Particle geometric properties
            word shapeName;                     // Particle shape: slab/cylinder/sphere
            scalar length;                      // Length of the particle [m]
            scalar width;                       // Width of the particle [m]

        // Particle components density and yield

            scalar moistureContent;             // Moisture content [-]
            scalar wetSolidDensity;             // Measured density of the wet solid [kg/m3]
            scalar drySolidDensity;             // Density of the dry solid [kg/m3]
            scalar charDensity;                 // Mass density of char [kg/m3]
            scalar ashDensity;                  // Mass density of ash [kg/m3]   

            scalar drySolidYield;               // Mass yield of dry solid in drying reaction [-]
            scalar thermalCharYield;            // Mass yield of char in thermal pyolysis [-]
            scalar oxidativeCharYield;          // Mass yield of char in oxidative pyolysis [-]
            scalar ashYield;                    // Mass yield of ash in char oxidation [-]

        // User-defined conductivity, heat capacity, emissivity, porosity & permeabiltiy  

            scalar k0_ws;       // Conductivity of wet solid (at 300 K) [W/m/K]
            scalar k0_ds;       // Conductivity of dry solid (at 300 K) [W/m/K]
            scalar k0_c;        // Conductivity of char (at 300 K) [W/m/K]
            scalar k0_a;        // Conductivity of ash (at 300 K) [W/m/K]
            scalar nk_ws;       // Temperature exponent of conductivity of wet solid
            scalar nk_ds;       // Temperature exponent of conductivity of dry solid
            scalar nk_c;        // Temperature exponent of conductivity of char
            scalar nk_a;        // Temperature exponent of conductivity of ash

            scalar gamma_ws;    // Effective radiation conductivity of wet solid [m]
            scalar gamma_ds;    // Effective radiation conductivity of dry solid [m]
            scalar gamma_c;     // Effective radiation conductivity of char [m]
            scalar gamma_a;     // Effective radiation conductivity of ash [m]

            scalar c0_ws;       // Heat capacity of wet solid [J/kg/K]
            scalar c0_ds;       // Heat capacity of dry solid [J/kg/K]
            scalar c0_c;        // Heat capacity of char [J/kg/K]
            scalar c0_a;        // Heat capacity of ash [J/kg/K]
            scalar nc_ws;       // Temperature exponent of heat capacity of wet solid
            scalar nc_ds;       // Temperature exponent of heat capacity of dry solid
            scalar nc_c;        // Temperature exponent of heat capacity of char
            scalar nc_a;        // Temperature exponent of heat capacity of ash

            scalar eps_ws;      // Surface emissivity of wet solid [-]
            scalar eps_ds;      // Surface emissivity of dry solid [-]
            scalar eps_c;       // Surface emissivity of char [-]
            scalar eps_a;       // Surface emissivity of ash [-]

            scalar psi_ws;      // Porosity of the particle when pure wet solid [-]
            scalar psi_ds;      // Porosity of the particle when pure solid [-]
            scalar psi_c;       // Porosity of the particle when pure char [-]
            scalar psi_a;       // Porosity of the particle when pure ash [-]

            scalar Kperm_ws;    // Permeability of the particle when pure wet solid [m2]
            scalar Kperm_ds;    // Permeability of the particle when pure solid [m2]
            scalar Kperm_c;     // Permeability of the particle when pure char [m2]
            scalar Kperm_a;     // Permeability of the particle when pure ash [m2]


           // Diagnostic Fields  

            volScalarField packingRatio;        // Packing Ratio [-]
            volScalarField gasVolRatio;         // 1-Packing Ratio [-]
            volScalarField surfaceTemp;         // Surface temperature [K]         
            volScalarField surfaceO2MassFrac;   // Surface O2 mass fraction [-]
            volVectorField momentum;            // Momentum [kg m/s]          


    // Constructors

        //- Construct given mesh
        biomassCloud
        (
            const fvMesh&,
            const word& cloudName,
            const volScalarField& rho,
            const volScalarField& YO2,
            const volScalarField& G,
            const volVectorField& U,
            const dimensionedVector& g,
            const SLGThermo& thermo,
            bool readFields = true
        );

    // Member Functions

        // Access

            bool hasWallImpactDistance() const;

            inline const fvMesh& mesh() const;

            inline const SLGThermo& thermo() const;
            

            // Sources

            	// Mass

                    //- Return reference to mass transfer fields [kg/m3]
                    inline DimensionedField<scalar, volMesh>& massTrans();

                    //- Return const access to mass transfer fields [kg/m3]
                    inline const DimensionedField<scalar, volMesh>& massTrans() const;

                    //- Return total mass source term [kg/s/m3]
                    // - fully explicit    
                    inline tmp<DimensionedField<scalar, volMesh>> Srho() const;

                // Species Mass

                    //- Return reference to mass transfer for field i [kg/m3]
                    inline DimensionedField<scalar, volMesh>& speciesMassTrans(const label i);

                    //- Return reference to mass transfer fields [kg/m3]
                    inline PtrList<DimensionedField<scalar, volMesh>>& speciesMassTrans();

                    //- Return const access to mass transfer fields [kg/m3]
                    inline const PtrList<DimensionedField<scalar, volMesh>>&  speciesMassTrans() const;

                    //- Return mass source term for specie i - specie eqn [kg/s/m3]   
                    // - fully explicit                        
                    inline tmp<DimensionedField<scalar, volMesh>> SYi(const label i) const;                 

                // Momentum

                    //- Return reference to momentum source [kg m/s /m3]
                    inline DimensionedField<vector, volMesh>& momentumTrans();

                    //- Return const reference to momentum source [kg m/s /m3]
                    inline const DimensionedField<vector, volMesh>& momentumTrans() const;

                    //- Return momentum source term [kg m/s^2 /m3 = N/m3]
                    inline tmp<DimensionedField<vector, volMesh>> SU() const;

                // Energy

                    //- Return reference to convective heat transfer [J/m3]
                    inline DimensionedField<scalar, volMesh>& QconvTrans();

                    //- Return const reference to convective heat transfer [J/m3]
                    inline const DimensionedField<scalar, volMesh>&QconvTrans() const;

                    //- Return conevtive heat transfer source term [W/m3]
                    inline tmp<DimensionedField<scalar, volMesh>> Qconv() const; 

                // Radiation

                    //- Return reference to RTE absorption coefficient [1/m]
                    inline DimensionedField<scalar, volMesh>& absorptionCoeff();

                    //- Return const reference to RTE absorption coefficient [1/m]
                    inline const DimensionedField<scalar, volMesh>& absorptionCoeff() const;

                    //- Return reference to RTE emission [J/m3]
                    inline DimensionedField<scalar, volMesh>& emissionTrans();

                    //- Return const reference to RTE emission [J/m3]                    
                    inline const DimensionedField<scalar, volMesh>& emissionTrans() const;

                    //- Return RTE absorption coefficient [1/m]
                    inline tmp<volScalarField> ap() const;

                    //- Return RTE emission source term [W/m3]
                    inline tmp<volScalarField> Ep() const;


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
