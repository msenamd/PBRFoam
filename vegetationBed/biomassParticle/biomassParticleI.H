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
    biomassCloud& spc,
    const interpolationCellPoint<scalar>& rhoInterp,
    const interpolationCellPoint<vector>& UInterp,
    const interpolationCellPoint<scalar>& muInterp,
    const interpolationCellPoint<scalar>& TInterp,
    const interpolationCellPoint<scalar>& GInterp,
    const vector& g
)
:
    particle::TrackingData<biomassCloud>(spc),
    rhoInterp_(rhoInterp),
    UInterp_(UInterp),
    muInterp_(muInterp),
    TInterp_(TInterp),
    GInterp_(GInterp),
    g_(g)
{

}


inline Foam::biomassParticle::biomassParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const scalar delta,
    const vector& U,
    const bool status,
    const scalarField T,
    const scalarField x_m,
    const scalarField x_vs
)
:
    particle(mesh, position, celli, tetFacei, tetPti),
    delta_(delta),
    U_(U),
    status_(status),
    T_(T),
    x_m_(x_m),
    x_vs_(x_vs),
    Tsurf_(0.0),
    volProdRate_tot_(0.0),
    volProdRate_fuel_(0.0),
    volProdRate_H2O_(0.0),
    volHRR_(0.0),
    surfToVol_(0.0),
    mass_(0.0),
    vol_(0.0),
    hConv_(0.0),
    CD_(0.0)
{

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::interpolationCellPoint<Foam::scalar>&
Foam::biomassParticle::trackingData::rhoInterp() const
{
    return rhoInterp_;
}


inline const Foam::interpolationCellPoint<Foam::vector>&
Foam::biomassParticle::trackingData::UInterp() const
{
    return UInterp_;
}


inline const Foam::interpolationCellPoint<Foam::scalar>&
Foam::biomassParticle::trackingData::muInterp() const
{
    return muInterp_;
}

inline const Foam::interpolationCellPoint<Foam::scalar>&
Foam::biomassParticle::trackingData::TInterp() const
{
    return TInterp_;
}

inline const Foam::interpolationCellPoint<Foam::scalar>&
Foam::biomassParticle::trackingData::GInterp() const
{
    return GInterp_;
}

inline const Foam::vector& Foam::biomassParticle::trackingData::g() const
{
    return g_;
}

// Access
inline Foam::scalar Foam::biomassParticle::delta() const
{
    return delta_;
}

inline const Foam::vector& Foam::biomassParticle::U() const
{
    return U_;
}

inline bool Foam::biomassParticle::status() const
{
    return status_;
}

inline Foam::scalarField Foam::biomassParticle::T() const
{
    return T_;
}

inline Foam::scalarField Foam::biomassParticle::x_m() const
{
    return x_m_;
}

inline Foam::scalarField Foam::biomassParticle::x_vs() const
{
    return x_vs_;
}


inline Foam::scalar Foam::biomassParticle::Tsurf() const
{
    return Tsurf_;
}

inline Foam::scalar Foam::biomassParticle::volProdRate_tot() const
{
    return volProdRate_tot_;
}

inline Foam::scalar Foam::biomassParticle::volProdRate_fuel() const
{
    return volProdRate_fuel_;
}

inline Foam::scalar Foam::biomassParticle::volProdRate_H2O() const
{
    return volProdRate_H2O_;
}

inline Foam::scalar Foam::biomassParticle::volHRR() const
{
    return volHRR_;
}

// ************************************************************************* //
