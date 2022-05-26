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
        // Set the Lagrangian time-step
        scalar dt = min(dtMax, tEnd);

        // Cache the parcel current cell as this will change if a face is hit
        const label celli = cell();

        // Track particle to a given position and returns 1.0 if the
        // trajectory is completed without hitting a face otherwise
        // stops at the face and returns the fraction of the trajectory
        // completed.
        if (mag(particleVelo_) > ROOTVSMALL)
        {
            dt *= trackToFace(position() + dt*particleVelo_, td);                         
        }
        tEnd -= dt;
        stepFraction() = 1.0 - tEnd/trackTime;  

        // Calculate thermal degradation, avoid problems with extremely small timesteps
        if (dt > ROOTVSMALL)
        {
            updateParticle(td, dt, celli);
        }

        // Check if particle is consumed
        if (particleState_ == Particle::eState::consumed)
        {
            // remove the particle
            td.keepParticle = false;
        }
          
        // Check if a patch is hit 
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
    // Get the external gas phase properties of this cell
    // Note: assume very small volume fraction of solid in a cell
    label indexFuel = (td.cloud().thermo()).carrier().species()[td.cloud().fuel];
    label indexO2   = (td.cloud().thermo()).carrier().species()["O2"];
    label indexH2O  = (td.cloud().thermo()).carrier().species()["H2O"];
    label indexCO2  = (td.cloud().thermo()).carrier().species()["CO2"];

    cellPointWeight cpw(mesh_, position(), celli, face());

    scalar externalGasDensity     = td.rhoInterp().interpolate(position(), celli);
    scalar externalGasPressure    = td.pInterp().interpolate(position(), celli) - 101325;
    scalar externalGasTemp        = td.TInterp().interpolate(position(), celli);
    scalar externalO2MassFraction = td.YO2Interp().interpolate(position(), celli);
    scalar externalIrradiation    = td.GInterp().interpolate(position(), celli);
    vector externalGasVelo        = td.UInterp().interpolate(cpw);

    // Pre-computing step
    p1D = td.cloud().superParticle;

    particleTemp_std.resize(particleTemp_.size(), 0.0);
    wetSolidVolFraction_std.resize(wetSolidVolFraction_.size(), 0.0);
    drySolidVolFraction_std.resize(drySolidVolFraction_.size(), 0.0);
    charVolFraction_std.resize(charVolFraction_.size(), 0.0);
    ashVolFraction_std.resize(ashVolFraction_.size(), 0.0);
    particleO2MassFraction_std.resize(particleO2MassFraction_.size(), 0.0);
    particlePressure_std.resize(particlePressure_.size(), 0.0);
    integralWetSolidMass_std.resize(integralWetSolidMass_.size(), 0.0); 
    integralDrySolidMass_std.resize(integralDrySolidMass_.size(), 0.0); 
    integralCharMass_std.resize(integralCharMass_.size(), 0.0); 

    for (int j=0 ; j<particleTemp_.size() ; j++)
    {
        particleTemp_std[j]           = particleTemp_[j];
        wetSolidVolFraction_std[j]    = wetSolidVolFraction_[j];
        drySolidVolFraction_std[j]    = drySolidVolFraction_[j];
        charVolFraction_std[j]        = charVolFraction_[j];
        ashVolFraction_std[j]         = ashVolFraction_[j];
        particleO2MassFraction_std[j] = particleO2MassFraction_[j];
        particlePressure_std[j]       = particlePressure_[j];
        integralWetSolidMass_std[j]   = integralWetSolidMass_[j];
        integralDrySolidMass_std[j]   = integralDrySolidMass_[j];
        integralCharMass_std[j]       = integralCharMass_[j];
    }

    p1D.preStepForward(
                            (Particle::eState)particleState_,
                            particledt_,
                            particleSize_,
                            particleTemp_std, 
                            wetSolidVolFraction_std, 
                            drySolidVolFraction_std, 
                            charVolFraction_std, 
                            ashVolFraction_std,
                            particleO2MassFraction_std,
                            particlePressure_std,
                            integralWetSolidMass_std,
                            integralDrySolidMass_std,
                            integralCharMass_std
                        );

    // Compute thermo-chemical degradation for time step dt_
    p1D.stepForward(
                        dt_,
                        externalGasTemp,
                        mag(externalGasVelo - particleVelo_),
                        externalO2MassFraction,
                        externalGasPressure,
                        externalIrradiation
                    );

    // Update current particle data from model outputs
    particleTemp_.resize(p1D.numCells, 0.0);
    particlePressure_.resize(p1D.numCells, 0.0);
    particleO2MassFraction_.resize(p1D.numCells, 0.0);
    wetSolidVolFraction_.resize(p1D.numCells, 0.0);
    drySolidVolFraction_.resize(p1D.numCells, 0.0);
    charVolFraction_.resize(p1D.numCells, 0.0);
    ashVolFraction_.resize(p1D.numCells, 0.0);
    integralWetSolidMass_.resize(p1D.numCells, 0.0);
    integralDrySolidMass_.resize(p1D.numCells, 0.0);
    integralCharMass_.resize(p1D.numCells, 0.0);

    for (int j=0 ; j<particleTemp_.size() ; j++)
    {
        particleTemp_[j]           = p1D.Temp[j];
        particlePressure_[j]       = p1D.particlePressure[j];
        particleO2MassFraction_[j] = p1D.particleO2MassFraction[j];
        wetSolidVolFraction_[j]    = p1D.wetSolidVolFraction[j];
        drySolidVolFraction_[j]    = p1D.drySolidVolFraction[j];
        charVolFraction_[j]        = p1D.charVolFraction[j];
        ashVolFraction_[j]         = p1D.ashVolFraction[j];
        integralWetSolidMass_[j]   = p1D.integralWetSolidMass[j];
        integralDrySolidMass_[j]   = p1D.integralDrySolidMass[j];
        integralCharMass_[j]       = p1D.integralCharMass[j];
    }

    particleState_     = p1D.state;
    particledt_        = p1D.localTimeStepSize;
    particleSize_      = p1D.particleSize;
    particleMass_      = p1D.particleMass;
    particleVol_       = p1D.particleVol;
    surfaceTemp_       = p1D.surfaceTemp;
    coreTemp_          = p1D.coreTemp;
    convFlux_          = p1D.surfaceHeatFluxConv;
    radFlux_           = p1D.surfaceHeatFluxRad; 
    massFlux_          = p1D.surfaceMassFlux;   
    surfaceO2MassFrac_ = p1D.surfaceO2MassFrac;
    hConv_             = p1D.h_conv;
    dryingRate_        = p1D.globalR1reactionRate;
    pyrolysisRate_     = p1D.globalR2reactionRate;
    oxidPyrolysisRate_ = p1D.globalR3reactionRate;
    charOxidRate_      = p1D.globalR4reactionRate;
    massLossRate_      = p1D.globalMassLossRate;
    heatReleaseRate_   = p1D.globalHeatReleaseRate;

    // update particle drag
    if(td.cloud().dragModel == "constant")
    {
        CD_ = td.cloud().dragCoeff;
    }
    else if (particleState_ == Particle::eState::ashed)
    {
        CD_ = 1.0;
    }
    else
    {
        CD_ = p1D.dragCoeff;
    }

    // Update particle velocuty
    if (td.cloud().firebrands)
    {
        scalar B = 0.5 * externalGasDensity * CD_
                 * p1D.projectedAreaRatio * p1D.particleSurfToVolRatio * p1D.particleVol 
                 / p1D.particleMass * mag(externalGasVelo - particleVelo_);

        particleVelo_ = (particleVelo_ + dt_ * (B * externalGasVelo + td.g())) / (1.0 + B * dt_);
    }
    else if (particleState_ == Particle::eState::ashed && !onBoundary())
    {
        particleVelo_ = particleVelo_ + dt_ *  td.g();
    }      
    else
    {
        particleVelo_ = Zero;
    }

    //  Accumulate carrier phase source terms for this dt
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Packing Ratio (-)
    td.cloud().packingRatio[celli] = td.cloud().nParticles * p1D.particleVol / mesh_.cellVolumes()[celli];

    // Species (kg/s *1/m3 *s = kg/m3)
    td.cloud().rhoTrans(indexH2O)[celli]  += p1D.globalMoistureReleaseRate / p1D.particleVol * td.cloud().packingRatio()[celli] * dt_;
    td.cloud().rhoTrans(indexFuel)[celli] += p1D.globalGasFuelReleaseRate  / p1D.particleVol * td.cloud().packingRatio()[celli] * dt_;
    td.cloud().rhoTrans(indexCO2)[celli]  += p1D.globalCO2ReleaseRate      / p1D.particleVol * td.cloud().packingRatio()[celli] * dt_;

    td.cloud().rhoTrans(indexO2)[celli]   -= p1D.h_mass * (externalO2MassFraction - p1D.surfaceO2MassFrac)
                                            * p1D.particleSurfToVolRatio * td.cloud().packingRatio()[celli] * dt_;

    // Momentum (N/m3 *s)
    td.cloud().UTrans()[celli] -= 0.5 * externalGasDensity * CD_ * p1D.projectedAreaRatio 
                                * mag(externalGasVelo - particleVelo_) * (externalGasVelo - particleVelo_)
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

    // vegetaion bed state (-)
    td.cloud().state[celli] = particleState_;

    // vegetaion bed mass (kg)
    td.cloud().mass[celli] = td.cloud().nParticles * p1D.particleMass;

    // Mass rates (kg/s/m3)
    td.cloud().massLossRatePUVbed[celli]        = td.cloud().nParticles * p1D.globalMassLossRate / mesh_.cellVolumes()[celli];
    td.cloud().gasFuelReleaseRatePUVbed[celli]  = td.cloud().nParticles * p1D.globalGasFuelReleaseRate / mesh_.cellVolumes()[celli];
    td.cloud().dryingRatePUVbed[celli]          = td.cloud().nParticles * p1D.globalR1reactionRate / mesh_.cellVolumes()[celli];
    td.cloud().pyrolysisRatePUVbed[celli]       = td.cloud().nParticles * p1D.globalR2reactionRate / mesh_.cellVolumes()[celli];
    td.cloud().oxidPyrolysisRatePUVbed[celli]   = td.cloud().nParticles * p1D.globalR3reactionRate / mesh_.cellVolumes()[celli];
    td.cloud().charOxidRatePUVbed[celli]        = td.cloud().nParticles * p1D.globalR4reactionRate / mesh_.cellVolumes()[celli];

    // Surface temperature (K)
    td.cloud().surfaceTemp[celli] = surfaceTemp_;

    // Surface heat fluxes per unit area of vegetation bed (W/m2)
    td.cloud().surfaceHeatFluxConv[celli] = p1D.surfaceHeatFluxConv * td.cloud().nParticles;
    td.cloud().surfaceHeatFluxRad[celli]  = p1D.surfaceHeatFluxRad * td.cloud().nParticles;

    // Surface O2 mass fraction (-)
    td.cloud().surfaceO2MassFrac[celli] = surfaceO2MassFrac_;

    // Momentum (kg m/s)
    td.cloud().momentum[celli] = td.cloud().nParticles * p1D.particleMass * particleVelo_;


    //  Wrtie diagnostic files
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    p1D.pCurrentTime = &mesh_.time().value();


    if(particleID_ > 0)
    {
        std::ofstream outParticle(outputPath/ "particle" + name(particleID_) + ".csv", ios::app);

        outParticle << mesh_.time().value() << "," << this->position()[0] << "," << this->position()[1] << "," << this->position()[2] << "," 
                    << particleState_ << "," << particledt_ << "," << particleSize_ << ","  << p1D.particleSurfToVolRatio << "," 
                    << particleVelo_[0] << "," << particleVelo_[1] << "," << particleVelo_[2] << ","
                    << particleMass_ << "," << surfaceTemp_  << "," << coreTemp_  << ","
                    << convFlux_ << "," << radFlux_ << "," << massFlux_ << ","
                    << surfaceO2MassFrac_ << "," << hConv_ << "," << CD_ << ","
                    << dryingRate_ << "," << pyrolysisRate_ << "," << oxidPyrolysisRate_ << ","
                    << charOxidRate_ << "," << massLossRate_ << "," << heatReleaseRate_ << "\n";


        std::ofstream TempFile(outputPath/ "particle" + name(particleID_) + "_temperature.csv", ios::app);
        p1D.writeTempLine(TempFile);

        std::ofstream wetSolidFile(outputPath/ "particle" + name(particleID_) + "_wetSolid.csv", ios::app);
        p1D.writeWetSolidLine(wetSolidFile);

        std::ofstream drySolidFile(outputPath/ "particle" + name(particleID_) + "_drySolid.csv", ios::app);
        p1D.writeDrySolidLine(drySolidFile);

        std::ofstream charFile(outputPath/ "particle" + name(particleID_) + "_char.csv", ios::app);
        p1D.writeCharLine(charFile);

        std::ofstream ashFile(outputPath/ "particle" + name(particleID_) + "_ash.csv", ios::app);
        p1D.writeAshLine(ashFile);

        std::ofstream O2File(outputPath/ "particle" + name(particleID_) + "_O2.csv", ios::app);
        p1D.writeO2Line(O2File);

        std::ofstream pressureFile(outputPath/ "particle" + name(particleID_) + "_pressure.csv", ios::app);
        p1D.writePressureLine(pressureFile); 

        std::ofstream coordFile(outputPath/ "particle" + name(particleID_) + "_cellCenter.csv", ios::app);
        p1D.writeCoordLine(coordFile);                                                       
    }
}

// ************************************************************************* //
