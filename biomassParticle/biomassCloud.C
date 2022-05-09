/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "biomassCloud.H"
#include "fvMesh.H"
#include "volFields.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::biomassCloud::biomassCloud
(
    const fvMesh& mesh,
    const word& cloudName,
    const volScalarField& rho,
    const volScalarField& YO2,
    const volScalarField& G,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& thermo,    
    bool readFields
)
:
    Cloud<biomassParticle>(mesh, cloudName, false),

    mesh_(mesh),

    particleProperties_
    (
        IOobject
        (
            "biomassProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    setParticlesDict
    (
        IOobject
        (
            "setParticlesDict",
            mesh_.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    
    rho_(rho),
    YO2_(YO2),
    G_(G),
    U_(U),
    thermo_(thermo),
    g_(g),

    rhoTrans_
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + "_rhoTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimMass/dimVol, 0.0)
        )
    ),
    rhoYTrans_
    (
        thermo.carrier().species().size()
    ),
    UTrans_
    (
        new DimensionedField<vector, volMesh>
        (
            IOobject
            (
                this->name() + "_UTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimForce/dimVol*dimTime, Zero)
        )
    ),
    QconvTrans_
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + "_QconvTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimVol, 0.0)
        )
    ),
    absorptionCoeff_
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + "_absorptionCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimLength , 0.0)
        )
    ),        
    emissionTrans_
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + "_emissionTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimEnergy/dimVol , 0.0)
        )
    ),

    wetSolid(),
    drySolid(),
    Char(),
    ash(),
    R1(),
    R2(),
    R3(),
    R4(),

    nParticles(readScalar(particleProperties_.lookup("nParticlesPerCell"))),
    fuel(particleProperties_.lookup("fuelComposition")),

    firebrands(readBool(particleProperties_.lookup("firebrands"))),
    particleElasticity(readScalar(particleProperties_.lookup("particleElasticity"))),
    particleViscosity(readScalar(particleProperties_.lookup("particleViscosity"))),
    dragModel(particleProperties_.lookup("dragModel")),
    dragCoeff(readScalar(particleProperties_.lookup("dragCoeff"))),
    
    maxIter(readLabel(particleProperties_.lookup("maxIter"))),
    meshResolution(readScalar(particleProperties_.lookup("meshResolution"))),
    timeStepThreshold(readScalar(particleProperties_.lookup("timeStepThreshold"))),
    temperatureThreshold(readScalar(particleProperties_.lookup("temperatureThreshold"))),
    solidSpeciesThreshold(readScalar(particleProperties_.lookup("solidSpeciesThreshold"))),
    O2Threshold(readScalar(particleProperties_.lookup("O2Threshold"))),
    pressureThreshold(readScalar(particleProperties_.lookup("pressureThreshold"))),
    temperatureURF(readScalar(particleProperties_.lookup("temperatureURF"))),
    solidSpeciesURF(readScalar(particleProperties_.lookup("solidSpeciesURF"))),
    O2URF(readScalar(particleProperties_.lookup("O2URF"))),
    pressureURF(readScalar(particleProperties_.lookup("pressureURF"))),
    flagRemeshing(readBool(particleProperties_.lookup("flagRemeshing"))),

    shapeName(particleProperties_.lookup("shapeName")),
    length(readScalar(particleProperties_.lookup("length"))),
    width(readScalar(particleProperties_.lookup("width"))),
    
    wetSolidDensity(readScalar(particleProperties_.lookup("wetSolidDensity"))),
    drySolidDensity(readScalar(particleProperties_.lookup("drySolidDensity"))),
    charDensity(readScalar(particleProperties_.lookup("charDensity"))),
    ashDensity(readScalar(particleProperties_.lookup("ashDensity"))),
   
    k0_ws(readScalar(particleProperties_.lookup("k0_ws"))),
    k0_ds(readScalar(particleProperties_.lookup("k0_ds"))),
    k0_c(readScalar(particleProperties_.lookup("k0_c"))),
    k0_a(readScalar(particleProperties_.lookup("k0_a"))),
    nk_ws(readScalar(particleProperties_.lookup("nk_ws"))),
    nk_ds(readScalar(particleProperties_.lookup("nk_ds"))),
    nk_c(readScalar(particleProperties_.lookup("nk_c"))),
    nk_a(readScalar(particleProperties_.lookup("nk_a"))),

    gamma_ws(readScalar(particleProperties_.lookup("gamma_ws"))),
    gamma_ds(readScalar(particleProperties_.lookup("gamma_ds"))),
    gamma_c(readScalar(particleProperties_.lookup("gamma_c"))),
    gamma_a(readScalar(particleProperties_.lookup("gamma_a"))),

    c0_ws(readScalar(particleProperties_.lookup("c0_ws"))),
    c0_ds(readScalar(particleProperties_.lookup("c0_ds"))),
    c0_c(readScalar(particleProperties_.lookup("c0_c"))),
    c0_a(readScalar(particleProperties_.lookup("c0_a"))),
    nc_ws(readScalar(particleProperties_.lookup("nc_ws"))),
    nc_ds(readScalar(particleProperties_.lookup("nc_ds"))),
    nc_c(readScalar(particleProperties_.lookup("nc_c"))),
    nc_a(readScalar(particleProperties_.lookup("nc_a"))),

    eps_ws(readScalar(particleProperties_.lookup("eps_ws"))),
    eps_ds(readScalar(particleProperties_.lookup("eps_ds"))),
    eps_c(readScalar(particleProperties_.lookup("eps_c"))),
    eps_a(readScalar(particleProperties_.lookup("eps_a"))),

    psi_ws(readScalar(particleProperties_.lookup("psi_ws"))),
    psi_ds(readScalar(particleProperties_.lookup("psi_ds"))),
    psi_c(readScalar(particleProperties_.lookup("psi_c"))),
    psi_a(readScalar(particleProperties_.lookup("psi_a"))),

    Kperm_ws(readScalar(particleProperties_.lookup("Kperm_ws"))),
    Kperm_ds(readScalar(particleProperties_.lookup("Kperm_ds"))),
    Kperm_c(readScalar(particleProperties_.lookup("Kperm_c"))),
    Kperm_a(readScalar(particleProperties_.lookup("Kperm_a"))),

    A_R1(readScalar(particleProperties_.lookup("A_R1"))),
    Ea_R1(readScalar(particleProperties_.lookup("Ea_R1"))),
    n_R1(readScalar(particleProperties_.lookup("n_R1"))),
    DeltaH_R1(readScalar(particleProperties_.lookup("DeltaH_R1"))),
    drySolidYield(readScalar(particleProperties_.lookup("drySolidYield"))),

    A_R2(readScalar(particleProperties_.lookup("A_R2"))),
    Ea_R2(readScalar(particleProperties_.lookup("Ea_R2"))),
    n_R2(readScalar(particleProperties_.lookup("n_R2"))),
    DeltaH_R2(readScalar(particleProperties_.lookup("DeltaH_R2"))),
    thermalCharYield(readScalar(particleProperties_.lookup("thermalCharYield"))),
    
    A_R3(readScalar(particleProperties_.lookup("A_R3"))),
    Ea_R3(readScalar(particleProperties_.lookup("Ea_R3"))),
    n_R3(readScalar(particleProperties_.lookup("n_R3"))),
    nO2_R3(readScalar(particleProperties_.lookup("nO2_R3"))),
    DeltaH_R3(readScalar(particleProperties_.lookup("DeltaH_R3"))),
    oxidativeCharYield(readScalar(particleProperties_.lookup("oxidativeCharYield"))),
    
    A_R4(readScalar(particleProperties_.lookup("A_R4"))),
    Ea_R4(readScalar(particleProperties_.lookup("Ea_R4"))),
    n_R4(readScalar(particleProperties_.lookup("n_R4"))),
    nO2_R4(readScalar(particleProperties_.lookup("nO2_R4"))),
    DeltaH_R4(readScalar(particleProperties_.lookup("DeltaH_R4"))),
    ashYield(readScalar(particleProperties_.lookup("ashYield"))),

    packingRatio
    (
        IOobject
        (
            this->name() + "_packingRatio",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0, 0 ,0) , 0.0)
    ),  
    state
    (
        IOobject
        (
            this->name() + "_state",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0, 0 ,0) , 0.0)
    ),    
    massLossRatePUVbed
    (
        IOobject
        (
            this->name() + "_massLossRatePUVbed",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMass/dimTime , 0.0)
    ), 
    gasFuelReleaseRatePUVbed
    (
        IOobject
        (
            this->name() + "_gasFuelReleaseRatePUVbed",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMass/dimTime , 0.0)
    ),
    dryingRatePUVbed
    (
        IOobject
        (
            this->name() + "_dryingRatePUVbed",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMass/dimTime , 0.0)
    ),
    pyrolysisRatePUVbed
    (
        IOobject
        (
            this->name() + "_pyrolysisRatePUVbed",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMass/dimTime , 0.0)
    ),
    oxidPyrolysisRatePUVbed
    (
        IOobject
        (
            this->name() + "_oxidPyrolysisRatePUVbed",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMass/dimTime , 0.0)
    ), 
    charOxidRatePUVbed
    (
        IOobject
        (
            this->name() + "_charOxidRatePUVbed",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMass/dimTime , 0.0)
    ),                                     
    surfaceTemp
    (
        IOobject
        (
            this->name() + "_surfaceTemp",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 1, 0, 0 ,0) , 0.0)
    ), 
    surfaceHeatFluxConv
    (
        IOobject
        (
            this->name() + "_surfaceHeatFluxConv",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimEnergy/dimTime/dimArea , 0.0)
    ),
    surfaceHeatFluxRad
    (
        IOobject
        (
            this->name() + "_surfaceHeatFluxRad",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimEnergy/dimTime/dimArea , 0.0)
    ),                  
    surfaceO2MassFrac
    (
        IOobject
        (
            this->name() + "_surfaceO2MassFrac",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0, 0 ,0) , 0.0)
    ),
    momentum
    (
        IOobject
        (
            this->name() + "_momentum",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dimensionSet(1, 1, -1, 0, 0, 0 ,0) , Zero)
    )
{

    Info << "Reading biomass particle properties" << endl;

    if (readFields)
    {
        biomassParticle::readFields(*this);
    }
    
    Info << "   particle geometry is: " << shapeName << endl;
    Info << "   motion is: " << firebrands<< endl;

    Info << "   setting solid material and reaction properties " << endl;

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

    R1.set_reaction(A_R1, Ea_R1, n_R1, 0.0, DeltaH_R1, drySolidYield, 0.0);
    R2.set_reaction(A_R2, Ea_R2, n_R2, 0.0, DeltaH_R2, thermalCharYield, 0.0);
    R3.set_reaction(A_R3, Ea_R3, n_R3, nO2_R3, DeltaH_R3, oxidativeCharYield, 0.1*(1.0-oxidativeCharYield));
    R4.set_reaction(A_R4, Ea_R4, n_R4, nO2_R4, DeltaH_R4, ashYield, 2.0*(1.0-ashYield));

    superParticle.air = &air;
    superParticle.wetSolid = &wetSolid;
    superParticle.drySolid = &drySolid;
    superParticle.Char = &Char;
    superParticle.ash = &ash;

    superParticle.R1 = &R1;
    superParticle.R2 = &R2;
    superParticle.R3 = &R3;
    superParticle.R4 = &R4;

    Info << "   setting particle geometry " << endl;
    superParticle.setGeometry(
                                shapeName, 
                                readScalar(setParticlesDict.lookup("initSize")), 
                                length, 
                                width
                            );

    Info << "   setting particle solution control " << endl;
    superParticle.setSolutionControl(
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
                                    pressureURF,
                                    flagRemeshing
                                );

    Info << "   setting particle information at time Zero  " << endl;
    superParticle.initialize(
                                readScalar(setParticlesDict.lookup("initTemp")),
                                readScalar(setParticlesDict.lookup("initWetSolid")),
                                readScalar(setParticlesDict.lookup("initDrySolid")),
                                readScalar(setParticlesDict.lookup("initChar")),
                                readScalar(setParticlesDict.lookup("initAsh")),
                                readScalar(setParticlesDict.lookup("initYO2")),
                                readScalar(setParticlesDict.lookup("initGaugeP"))
                            );

    // Set storage for mass source fields and initialise to zero
    forAll(rhoYTrans_, i)
    {
        const word& specieName = thermo.carrier().species()[i];
        rhoYTrans_.set
        (
            i,
            new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    this->name() + "_rhoYTrans_" + specieName,
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh(),
                dimensionedScalar("zero", dimMass/dimVol, 0.0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::biomassCloud::hasWallImpactDistance() const
{
    return true;
}


void Foam::biomassCloud::evolve()
{

    interpolationCellPoint<scalar> rhoInterp(rho_);
    interpolationCellPoint<scalar> YO2Interp(YO2_);
    interpolationCellPoint<scalar> GInterp(G_);     
    interpolationCellPoint<vector> UInterp(U_);
    interpolationCellPoint<scalar> TInterp(thermo_.thermo().T());
    interpolationCellPoint<scalar> pInterp(thermo_.thermo().p());

    biomassParticle::trackingData
        td(*this, rhoInterp, YO2Interp, GInterp, UInterp, TInterp, pInterp, g_.value());

    // reset source terms
    biomassCloud::resetSourceTerms();    

    // evolve the cloud
    Cloud<biomassParticle>::move(td, mesh_.time().deltaTValue());

    Info << "packing ratio: " << "min = " << min(packingRatio()).value() 
    << " , " << "max = " << max(packingRatio()).value() << endl;

	Info << "surface temperature: " << "min = " << min(surfaceTemp()).value()
	<< " , " << "max = " << max(surfaceTemp()).value() << endl;
}


void Foam::biomassCloud::resetSourceTerms()
{
    rhoTrans().field() = 0.0;

    forAll(rhoYTrans_, i)
    {
        rhoYTrans_[i].field() = 0.0;
    }

    UTrans().field() = Zero;

    QconvTrans().field() = 0.0;

    absorptionCoeff().field() = 0.0;
    emissionTrans().field() = 0.0;
}



// ************************************************************************* //
