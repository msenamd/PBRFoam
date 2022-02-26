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

#include "biomassCloud.H"
#include "constants.H"
#include "physicoChemicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<biomassParticle>, 0);
}

// * * * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

bool Foam::biomassParticle::move
(
	trackingData& td,
	const scalar trackTime
)
{

    td.switchProcessor = false;
    td.keepParticle = true;

    const polyBoundaryMesh& pbMesh = mesh_.boundaryMesh();

    scalar tEnd = (1.0 - stepFraction())*trackTime;
    scalar dtMax = tEnd;

    while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
   	{
        // set the lagrangian time-step
        scalar dt = min(dtMax, tEnd);

        // remember which cell the parcel is in
        // since this will change if a face is hit
        label celli = cell();

        // Track particle to a given position and returns 1.0 if the
        // trajectory is completed without hitting a face otherwise
        // stops at the face and returns the fraction of the trajectory
        // completed.
        dt *= trackToFace(position() + dt*particleVelo_, td);

        tEnd -= dt;
        stepFraction() = 1.0 - tEnd/trackTime;        

        // update particle variables
        updateParticle(td, dt, celli);

        // check if a patch is hit 
        if (onBoundary() && td.keepParticle)
        {
            if (isA<processorPolyPatch>(pbMesh[patch(face())]))
            {
                td.switchProcessor = true;
            }
        }
    }
    
    return td.keepParticle;
}


bool Foam::biomassParticle::hitPatch
(
    const polyPatch&,
    trackingData&,
    const label,
    const scalar,
    const tetIndices&
)
{
    return false;
}


void Foam::biomassParticle::hitProcessorPatch
(
    const processorPolyPatch&,
    trackingData& td
)
{
    td.switchProcessor = true;
}


void Foam::biomassParticle::hitWallPatch
(
    const wallPolyPatch& wpp,
    trackingData& td,
    const tetIndices& tetIs
)
{
    vector nw = tetIs.faceTri(mesh_).normal();
    nw /= mag(nw);

    scalar Un = particleVelo_ & nw;
    vector Ut = particleVelo_ - Un*nw;

    if (Un > 0)
    {
        particleVelo_ -= (1.0 + td.cloud().particleElasticity)*Un*nw;
    }

    particleVelo_ -= td.cloud().particleViscosity*Ut;
}


void Foam::biomassParticle::hitPatch
(
    const polyPatch&,
    trackingData& td
)
{
    td.keepParticle = false;
}


void Foam::biomassParticle::transformProperties (const tensor& T)
{
    particle::transformProperties(T);
    particleVelo_ = transform(T, particleVelo_);
}


void Foam::biomassParticle::transformProperties(const vector& separation)
{
    particle::transformProperties(separation);
}


Foam::scalar Foam::biomassParticle::wallImpactDistance(const vector&) const
{
    return particleSize_;
}


// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

void Foam::biomassParticle::updateParticle
(
    trackingData& td, 
    const scalar dt_,
    const label celli
)
{
    // get the external gas phase properties of this cell
    // Note: density and pressure are of gas only (not gas and solid)
    // solution obtained from continuum solver is for bulk density

    cellPointWeight cpw(mesh_, position(), celli, face());

    label indexFuel = (td.cloud().thermo()).carrier().species()[td.cloud().fuel];
    label indexO2   = (td.cloud().thermo()).carrier().species()["O2"];
    label indexH2O  = (td.cloud().thermo()).carrier().species()["H2O"];
    label indexCO2  = (td.cloud().thermo()).carrier().species()["CO2"];

    scalar externalGasDensity     = td.rhoInterp().interpolate(cpw) / (1.0 - td.cloud().packingRatio[celli]);
    scalar externalGasPressure    = (td.pInterp().interpolate(cpw) - 101325) / (1.0 - td.cloud().packingRatio[celli]);
    scalar externalGasTemp        = td.TInterp().interpolate(cpw);
    vector externalGasVelo        = td.UInterp().interpolate(cpw);
    scalar externalO2MassFraction = td.YO2Interp().interpolate(cpw);
    scalar externalIrradiation    = td.GInterp().interpolate(cpw);

    // relative velocity
    vector relativeVelo = externalGasVelo - particleVelo_;

    // copy material, reaction and solution data
    p1D = td.cloud().testParticle;

    // update geometrical properties
    p1D.setGeometry(td.cloud().shapeName, particleSize_, td.cloud().length, td.cloud().width);

    // initialize fields to particle model std vector format
    particleTemp_std.resize(particleTemp_.size());
    wetSolidVolFraction_std.resize(wetSolidVolFraction_.size());
    drySolidVolFraction_std.resize(drySolidVolFraction_.size());
    charVolFraction_std.resize(charVolFraction_.size());
    ashVolFraction_std.resize(ashVolFraction_.size());
    particleO2MassFraction_std.resize(particleO2MassFraction_.size());
    particlePressure_std.resize(particlePressure_.size());

    for (int j=0 ; j<particleTemp_.size() ; j++)
    {
        particleTemp_std[j]           = particleTemp_[j];
        wetSolidVolFraction_std[j]    = wetSolidVolFraction_[j];
        drySolidVolFraction_std[j]    = drySolidVolFraction_[j];
        charVolFraction_std[j]        = charVolFraction_[j];
        ashVolFraction_std[j]         = ashVolFraction_[j];
        particleO2MassFraction_std[j] = particleO2MassFraction_[j];
        particlePressure_std[j]       = particlePressure_[j];
    }

    p1D.initialize(
                    particleTemp_std, 
                    wetSolidVolFraction_std, 
                    drySolidVolFraction_std, 
                    charVolFraction_std, 
                    ashVolFraction_std,
                    particleO2MassFraction_std,
                    particlePressure_std
                );

    // compute thermo-chemical degradation
    p1D.stepForward(
                    dt_,
                    externalGasTemp,
                    mag(relativeVelo),
                    externalO2MassFraction,
                    externalGasPressure,
                    externalIrradiation
                );

    particleTemp_.resize(p1D.numCells);
    particlePressure_.resize(p1D.numCells);
    particleO2MassFraction_.resize(p1D.numCells);
    wetSolidVolFraction_.resize(p1D.numCells);
    drySolidVolFraction_.resize(p1D.numCells);
    charVolFraction_.resize(p1D.numCells);
    ashVolFraction_.resize(p1D.numCells);

    for (int j=0 ; j<particleTemp_.size() ; j++)
    {
        particleTemp_[j]           = p1D.Temp[j];
        particlePressure_[j]       = p1D.particlePressure[j];
        particleO2MassFraction_[j] = p1D.particleO2MassFraction[j];
        wetSolidVolFraction_[j]    = p1D.wetSolidVolFraction[j];
        drySolidVolFraction_[j]    = p1D.drySolidVolFraction[j];
        charVolFraction_[j]        = p1D.charVolFraction[j];
        ashVolFraction_[j]         = p1D.ashVolFraction[j];
    }

    particleState_     = p1D.state;
    particleSize_      = p1D.particleSize;
    surfaceTemp_       = p1D.surfaceTemp;
    surfaceO2MassFrac_ = p1D.surfaceO2MassFrac;
    hConv_             = p1D.h_conv;

    if(td.cloud().dragModel == "constant")
    {
        CD_ = td.cloud().dragCoeff;
    }
    else
    {
        CD_ = p1D.dragCoeff;
    } 

    if (td.cloud().firebrands)
    {
        scalar B = 0.5 * externalGasDensity * CD_
                 * p1D.projectedAreaRatio * p1D.particleSurfToVolRatio * p1D.particleVol 
                 / p1D.particleMass * mag(relativeVelo);

        particleVelo_ = (particleVelo_ + dt_ * (B * externalGasVelo + td.g())) / (1.0 + B * dt_);
    }
    else
    {
        particleVelo_ = Zero;
    }


    //  Accumulate carrier phase source terms for this dt
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Packing Ratio (-)
    td.cloud().packingRatio[celli] = td.cloud().nParticles * p1D.particleVol / mesh_.cellVolumes()[celli];
    td.cloud().gasVolRatio[celli] = 1.0 - td.cloud().packingRatio[celli];

    // Mass (kg/s/m3 *s = kg/m3)
    td.cloud().rhoTrans()[celli] +=  p1D.globalMassLossRate / p1D.particleVol * td.cloud().packingRatio()[celli] * dt_;

    // Species (kg/s/m3 *s = kg/m3)
    td.cloud().rhoYTrans(indexH2O)[celli]  += p1D.globalMoistureReleaseRate / p1D.particleVol * td.cloud().packingRatio()[celli] * dt_;
    td.cloud().rhoYTrans(indexFuel)[celli] += p1D.globalGasFuelReleaseRate  / p1D.particleVol * td.cloud().packingRatio()[celli] * dt_;
    td.cloud().rhoYTrans(indexCO2)[celli]  += p1D.globalCO2ReleaseRate      / p1D.particleVol * td.cloud().packingRatio()[celli] * dt_;

    td.cloud().rhoYTrans(indexO2)[celli]   -= p1D.h_mass * (externalO2MassFraction - p1D.surfaceO2MassFrac)
                                            * p1D.particleSurfToVolRatio * td.cloud().packingRatio()[celli] * dt_;

    // Momentum (N/m3 *s)
    td.cloud().UTrans()[celli] -= 0.5 * CD_ * p1D.projectedAreaRatio * externalGasDensity * mag(relativeVelo) * relativeVelo
                                * p1D.particleSurfToVolRatio * td.cloud().packingRatio()[celli] * dt_;

    // Energy (W/m3 *s)
    td.cloud().QconvTrans()[celli] -= p1D.h_conv * (externalGasTemp - surfaceTemp_)
                                    * p1D.particleSurfToVolRatio * td.cloud().packingRatio()[celli] * dt_;

    // RTE absorption coefficient (1/m)
    td.cloud().absorptionCoeff()[celli] =  p1D.particleSurfToVolRatio * td.cloud().packingRatio()[celli] * p1D.surfaceEmissivity / 4.0;

    // RTE emission (W/m3 *s)
    td.cloud().emissionTrans()[celli] += 4.0 * td.cloud().absorptionCoeff()[celli] * physicoChemical::sigma.value()
                                        * Foam::pow(surfaceTemp_, 4.0) * dt_;


    //  Additional diagnostic fields
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //Surface temperature (K)
    td.cloud().surfaceTemp[celli] = surfaceTemp_;

    //Surface O2 mass fraction (-)
    td.cloud().surfaceO2MassFrac[celli] = surfaceO2MassFrac_;

    //Momentum (kg m/s)
    td.cloud().momentum[celli] = td.cloud().nParticles * p1D.particleMass * particleVelo_;


}

// ************************************************************************* //
