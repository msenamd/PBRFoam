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
    const volScalarField& rho_g,
    const volVectorField& U_g,
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
    
    nParticles_(readScalar(particleProperties_.lookup("nParticles"))),

    fuel_(particleProperties_.lookup("fuelComposition")),

    firebrands_(particleProperties_.lookup("firebrands")),
    e_(readScalar(particleProperties_.lookup("e"))),
    mu_(readScalar(particleProperties_.lookup("mu"))),

    dragModel_(particleProperties_.lookup("dragModel")),
    dragCoeff_(readScalar(particleProperties_.lookup("dragCoeff"))),
    
    geometry_(particleProperties_.lookup("particleGeometry")),
    resolution_(readScalar(particleProperties_.lookup("resolution"))),
    cylinderLength_(readScalar(particleProperties_.lookup("cylinderLength"))),
    rectangleLength_(readScalar(particleProperties_.lookup("rectangleLength"))),
    rectangleWidth_(readScalar(particleProperties_.lookup("rectangleWidth"))),
    
    eta_c_(readScalar(particleProperties_.lookup("charYield"))),
    
    rho_m_(readScalar(particleProperties_.lookup("moistureDensity"))),
    rho_vs_(readScalar(particleProperties_.lookup("solidDensity"))),
    rho_c_(readScalar(particleProperties_.lookup("charDensity"))),

    c_m_(readScalar(particleProperties_.lookup("moistureHeatCapacity"))),
    c_vs_(readScalar(particleProperties_.lookup("solidHeatCapacity"))),
    c_c_(readScalar(particleProperties_.lookup("charHeatCapacity"))),

    k_m_(readScalar(particleProperties_.lookup("moistureConductivity"))),
    k_vs_(readScalar(particleProperties_.lookup("solidConductivity"))),
    k_c_(readScalar(particleProperties_.lookup("charConductivity"))),

    e_m_(readScalar(particleProperties_.lookup("moistureEmissivity"))),
    e_vs_(readScalar(particleProperties_.lookup("solidEmissivity"))),
    e_c_(readScalar(particleProperties_.lookup("charEmissivity"))),

    rho_g_(rho_g),
    U_g_(U_g),
    mu_g_(thermo.thermo().mu()),
    T_g_(thermo.thermo().T()),
    G_g_
    (
            IOobject
            (
                this->name() + "_G",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(1, 0, -3, 0, 0, 0 ,0) , 0.0)
    ),
    g_(g),
    thermo_(thermo),
    
    packingRatio_
    (
        new volScalarField
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
        )
    ),
    surfaceTemp_
    (
        new volScalarField
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
        )
    ),    
    surfaceToVolumeRatio_
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "_surfaceToVolumeRatio",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(0, -1, 0, 0, 0, 0 ,0) , 0.0)
        )
    ),
    momentum_
    (
        new volVectorField
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
    ),
    rhoTrans_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + "_rhoTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
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
        new volVectorField::Internal
        (
            IOobject
            (
                this->name() + "_UTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimForce/dimVol*dimTime, Zero)
        )
    ),
    Qr3Trans_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + "_Qr3Trans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimVol, 0.0)
        )
    ),
    QconvTrans_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + "_QconvTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimVol, 0.0)
        )
    ),
    absorption_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + "_absorptionCoeff",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimLength , 0.0)
        )
    ),        
    emissionTrans_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + "_emissionTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimMass/dimLength/pow4(dimTime)  , 0.0)
        )
    )
{

    if (readFields)
    {
        biomassParticle::readFields(*this);
    }
    Info << "Reading biomass properties" << endl;

    if (mesh_.objectRegistry::foundObject<volScalarField>("G"))
	{
		G_g_=mesh_.objectRegistry::template lookupObject<volScalarField>("G");
		Info << "	Found radiation field G" << endl;
	}
	else
	{
		Info << "	Warning: No radiation field G found, ignoring radiation transfer" << endl;
	}
    
    Info << "   particle geometry is: " << geometry_ << endl;
    Info << "   motion is: " << firebrands_<< endl;
    Info << "   particle grid resolution is: " << resolution_ << endl;
    Info << "   volatile fuel composition is: " << fuel_ << endl;
    Info << "   virgin solid density = " << rho_vs_ << endl;
    Info << "   virgin solid heat capacity = " << c_vs_ << endl;
    Info << "   virgin solid conductivity = " << k_vs_ << endl;
    Info << "   virgin solid emissivity = " << e_vs_ << endl;


    // Set storage for mass source fields and initialise to zero
    forAll(rhoYTrans_, i)
    {
        const word& specieName = thermo.carrier().species()[i];
        rhoYTrans_.set
        (
            i,
            new volScalarField::Internal
            (
                IOobject
                (
                    this->name() + "_rhoYTrans_" + specieName,
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
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

    interpolationCellPoint<scalar> rhoInterp(rho_g_);
    interpolationCellPoint<vector> UInterp(U_g_);
    interpolationCellPoint<scalar> muInterp(mu_g_);
    interpolationCellPoint<scalar> TInterp(T_g_);

    if (mesh().objectRegistry::foundObject<volScalarField>("G"))
	{
		G_g_=mesh().objectRegistry::template lookupObject<volScalarField>("G");
	}

    //must be divided by 4 because G = int(I*omega)
    interpolationCellPoint<scalar> GInterp(G_g_); 

    biomassParticle::trackingData
        td(*this, rhoInterp, UInterp, muInterp, TInterp, GInterp, g_.value());


    // reset source terms
    biomassCloud::resetSourceTerms();    

    //evolveCloud
    Cloud<biomassParticle>::move(td, mesh_.time().deltaTValue());

	Info << "packing ratio: " << "min= " << min(packingRatio()).value() << endl;
	Info << "packing ratio: " << "max= " << max(packingRatio()).value() << endl;
	Info << "surface temperature: " << "min= " << min(surfaceTemp()).value() << endl;
	Info << "surface temperature: " << "max= " << max(surfaceTemp()).value() << endl;
}


void Foam::biomassCloud::resetSourceTerms()
{
    packingRatio().field() = 0.0;
    surfaceTemp().field() = 0.0;
    surfaceToVolumeRatio().field() = 0.0;
    momentum().field() = Zero;
    
    rhoTrans().field() = 0.0;

    UTrans().field() = Zero;

    forAll(rhoYTrans_, i)
    {
        rhoYTrans_[i].field() = 0.0;
    }

    Qr3Trans().field() = 0.0;
    QconvTrans().field() = 0.0;

    absorption().field() = 0.0;
    emissionTrans().field() = 0.0;
}



// ************************************************************************* //
