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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::biomassParticle::trackingData::trackingData
(
    biomassCloud& bC,
    const interpolationCell<scalar>& rhoInterp,
    const interpolationCell<scalar>& YO2Interp,
    const interpolationCell<scalar>& GInterp,
    const interpolationCellPoint<vector>& UInterp,
    const interpolationCell<scalar>& TInterp,
    const interpolationCell<scalar>& pInterp,
    const vector& g
)
:
    particle::TrackingData<biomassCloud>(bC),
    rhoInterp_(rhoInterp),
    YO2Interp_(YO2Interp),
    GInterp_(GInterp),
    UInterp_(UInterp),
    TInterp_(TInterp),
    pInterp_(pInterp),
    g_(g)
{

}


inline Foam::biomassParticle::biomassParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    const label tetFacei,
    const label tetPtI,
    const label particleID,
    const label particleState,
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
)
:
    particle(mesh, position, celli, tetFacei, tetPtI),

    particleID_(particleID),    
    particleState_(particleState),
    particledt_(particledt),
    particleSize_(particleSize),
    particleVelo_(particleVelo),

    particleTemp_(particleTemp),
    particlePressure_(particlePressure),
    particleO2MassFraction_(particleO2MassFraction),    
    wetSolidVolFraction_(wetSolidVolFraction),
    drySolidVolFraction_(drySolidVolFraction),
    charVolFraction_(charVolFraction),
    ashVolFraction_(ashVolFraction),

    integralWetSolidMass_(integralWetSolidMass),
    integralDrySolidMass_(integralDrySolidMass),
    integralCharMass_(integralCharMass),
    
    particleMass_(0.0),
    particleVol_(0.0), 
    surfaceTemp_(0.0),
    coreTemp_(0.0),
    convFlux_(0.0),
    radFlux_(0.0),
    massFlux_(0.0),
    surfaceO2MassFrac_(0.0),
    hConv_(0.0),
    CD_(0.0),
    dryingRate_(0.0),
    pyrolysisRate_(0.0),
    oxidPyrolysisRate_(0.0),
    charOxidRate_(0.0),
    massLossRate_(0.0),
    heatReleaseRate_(0.0),
    outputPath("")
{

    // Setting the output path and writing headers
    if (Pstream::parRun())
    {
        outputPath = mesh.time().path()/".."/"particlePostProcessing";
    }
    else
    {
        outputPath = mesh.time().path()/"particlePostProcessing";
    } 
    mkDir(outputPath);

    if(particleID_ > 0)
    {
        std::ofstream outParticle(outputPath/ "particle" + name(particleID_) + ".csv", ios::app);

        outParticle << "time(s)" << "," << "x(m)" << "," << "y(m)" << "," << "z(m)" << "," 
                    << "state" << "," << "dt(s)" << "," << "size(m)" << "," << "surfToVol(1/m)" << ","  
                    << "Ux(m/s)" << "," << "Uy(m/s)" << "," << "Uz(m/s)" << ","
                    << "mass(kg)" << "," << "surfaceTemp(K)"  << "," << "coreTemp(K)"  << ","
                    << "convFlux(W/m2)" << "," << "radFlux(W/m2)" << "," << "massFlux(kg/s/m2)" << ","
                    << "surfaceO2MassFrac(-)" << "," << "hConv(W/m2/K)" << "," << "CD(-)" << ","
                    << "dryingRate(kg/s)" << "," << "pyrolysisRate(kg/s)" << "," << "oxidPyrolysisRate(kg/s)" << ","
                    << "charOxidRate(kg/s)" << "," << "massLossRate(kg/s)" << "," << "heatReleaseRate(W)" << "\n";        
    }

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::interpolationCell<Foam::scalar>&
Foam::biomassParticle::trackingData::rhoInterp() const
{
    return rhoInterp_;
}

inline const Foam::interpolationCell<Foam::scalar>&
Foam::biomassParticle::trackingData::YO2Interp() const
{
    return YO2Interp_;
}

inline const Foam::interpolationCell<Foam::scalar>&
Foam::biomassParticle::trackingData::GInterp() const
{
    return GInterp_;
}

inline const Foam::interpolationCellPoint<Foam::vector>&
Foam::biomassParticle::trackingData::UInterp() const
{
    return UInterp_;
}

inline const Foam::interpolationCell<Foam::scalar>&
Foam::biomassParticle::trackingData::TInterp() const
{
    return TInterp_;
}

inline const Foam::interpolationCell<Foam::scalar>&
Foam::biomassParticle::trackingData::pInterp() const
{
    return pInterp_;
}

inline const Foam::vector& Foam::biomassParticle::trackingData::g() const
{
    return g_;
}

// ************************************************************************* //
