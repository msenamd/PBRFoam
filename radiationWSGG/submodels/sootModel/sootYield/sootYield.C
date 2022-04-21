/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "sootYield.H"
#include "singleStepReactingMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::radiation::sootYield<ThermoType>::sootYield
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& modelType
)
:
    sootModel(dict, mesh, modelType),
    fv_
    (
        IOobject
        (
            "fv",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("fv", dimless, scalar(0.0))
    ),
    Ysoot_
    (
        IOobject
        (
            "Ysoot",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Ysoot", dimless, scalar(0.0))
    ),

    coeffsDict_(dict.subOrEmptyDict(modelType + "Coeffs")),
    sootYield_(readScalar(coeffsDict_.lookup("sootYield"))),
    thermo_(mesh.lookupObject<fluidThermo>(basicThermo::dictName)),
    thermo(mesh_.lookupObject<psiReactionThermo>("thermophysicalProperties"))
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::radiation::sootYield<ThermoType>::~sootYield()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::radiation::sootYield<ThermoType>::correct()
{

    forAll(fv_, cellI)
	{
       	fv_[cellI] = 0.0;
  	}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
