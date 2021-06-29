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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<biomassParticle>, 0);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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

    if(td.cloud().geometry() == "rectangle")
    {
		vol_=2.0*delta_*td.cloud().rectangleLength()*td.cloud().rectangleWidth();
    }
    else if(td.cloud().geometry() == "cylinder")
    {
		vol_= constant::mathematical::pi*Foam::pow(delta_,2.0)*td.cloud().cylinderLength();
    }
    else if(td.cloud().geometry() == "sphere")
    {
		vol_= 4.0/3.0*constant::mathematical::pi*Foam::pow(delta_,3.0);
    } 

    // note: dt changes if the particle is moving
    // trackToFace function has its own step fractions
    while (td.keepParticle && !td.switchProcessor && tEnd > SMALL)
   	{

        // set the lagrangian time-step
        scalar dt = min(dtMax, tEnd);

        // remember which cell the parcel is in
        // since this will change if a face is hit
        label celli = cell();
    
        // check relative size of bed to cell
        if (vol_ > mesh_.cellVolumes()[celli])
        {
            FatalErrorInFunction
                << "Particle size is greater than the mesh cell size" << nl
                << "cell volume is: " << mesh_.cellVolumes()[celli] << " , particle volume is " << vol_
                << exit(FatalError);
        }

        if (mag(U_) > SMALL)
        {

            dt *= trackToFace(position() + dt*U_, td);
        }
        
        tEnd -= dt;
        stepFraction() = 1.0 - tEnd/trackTime;        

        // Main calculation step
        if (dt > ROOTVSMALL)
        {
        	calcAll(td, dt, celli);
        }
        

        // Remove particle if status is burnt
        if(status_ == false)
        {
            td.keepParticle = false;
        }
        

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

    scalar Un = U_ & nw;
    vector Ut = U_ - Un*nw;

    if (Un > 0)
    {
        U_ -= (1.0 + td.cloud().e())*Un*nw;
    }

    U_ -= td.cloud().mu()*Ut;
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
    U_ = transform(T, U_);
}


void Foam::biomassParticle::transformProperties(const vector& separation)
{
    particle::transformProperties(separation);
}


Foam::scalar Foam::biomassParticle::wallImpactDistance(const vector&) const
{
    return delta_;
}


//------------ Calculation Funcitons ---------

void Foam::biomassParticle::calcAll
(
    trackingData& td, 
    const scalar dt_,
    const label celli
)
{

    // carrier phase properties
    cellPointWeight cpw(mesh_, position(), celli, face());
    scalar rhoc = td.rhoInterp().interpolate(cpw);
    scalar Tc = td.TInterp().interpolate(cpw);
    scalar Gc = td.GInterp().interpolate(cpw) / 4.0;
    vector Uc = td.UInterp().interpolate(cpw);
    scalar muc = td.muInterp().interpolate(cpw);
    word fuelComposition = td.cloud().fuel();
    label indexH2O= (td.cloud().thermo()).carrier().species()["H2O"];
    label indexFuel= (td.cloud().thermo()).carrier().species()[fuelComposition];
    label indexO2= (td.cloud().thermo()).carrier().species()["O2"];
    scalar YO2c = 0.0;  //needs update to get info from carrier phase


    // Sources transfer auxiliary variaables
    //~~~~~~~~
    vector dUTrans = Zero;


    // Thermo-chemical Degradatoin
    // ~~~~~~
    calcDegradation
        (
            td, 
            dt_,
            Tc,
            Uc,
            Gc,
            YO2c
        );


    // Motion
    // ~~~~~~

    word moveParticle = td.cloud().firebrands();

    U_ = calcVelocity(
                            td,
                            moveParticle, 
                            dt_, 
                            mass_,
                            delta_,
                            rhoc,
                            muc,
                            Uc,
                            dUTrans
                    );

    //Update the Packing Ratio and Surface to Volume Ratio
    td.cloud().packingRatio()[celli] = td.cloud().nParticles()*vol_ / mesh_.cellVolumes()[celli];
    td.cloud().surfaceToVolumeRatio()[celli] = surfToVol_;


    //  Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Note: multipy by dt_ here, then integrate over the cloud global dt in biomassCloudI.H
    // "Because dt_local  NOT Eq. dt_cloud" due to trackToFace function

        //mass (kg/s/m3 *s = kg/m3)
        td.cloud().rhoTrans()[celli] +=  volProdRate_tot_* td.cloud().packingRatio()[celli] * dt_;

        //Species (kg/s/m3 *s = kg/m3)
        td.cloud().rhoYTrans(indexH2O)[celli] += volProdRate_H2O_* td.cloud().packingRatio()[celli] * dt_;
        td.cloud().rhoYTrans(indexFuel)[celli] += volProdRate_fuel_* td.cloud().packingRatio()[celli] * dt_;

        //momentum (force/m3 *s)
        td.cloud().UTrans()[celli] -= dUTrans * surfToVol_ * td.cloud().packingRatio()[celli] * dt_;

        //energy (W/m3 *s)
        td.cloud().Qr3Trans()[celli] += volHRR_* (0.5) * td.cloud().packingRatio()[celli] * dt_; //assuming 50% of heat goes to gas phase
        td.cloud().QconvTrans()[celli] -= surfToVol_ * hConv_ * (Tc-Tsurf_) * td.cloud().packingRatio()[celli] * dt_; 

        //radiation
        td.cloud().emissionTrans()[celli] += (
                                            surfToVol_*td.cloud().packingRatio()[celli]/(4.0*(1.0-td.cloud().packingRatio()[celli])) 
                                            * 5.670367e-8 *Foam::pow(Tsurf_,4.0) / constant::mathematical::pi
                                        ) * dt_;

        td.cloud().absorption()[celli] =  surfToVol_ * td.cloud().packingRatio()[celli] / (4.0*(1.0-td.cloud().packingRatio()[celli]));


    // Update extra diagnostics
        //Momentum (kg m/s)
        td.cloud().momentum()[celli] = td.cloud().nParticles() * mass_ * U_;
        //Surface Temperature (K)
        td.cloud().surfaceTemp()[celli] = Tsurf_;
}



void Foam::biomassParticle::calcDegradation
(
    trackingData& td, 
    const scalar dt_,
    const scalar Tc,
    const vector Uc,
    const scalar Gc,
    const scalar YO2c
)
{

    //- Defining particle 1-D object
    particle_1D thisParticle;

	//passing inputs to PBM (converting to std::vector type)
	std::vector<double> T_std (T_.size());
	std::vector<double> x_m_std (x_m_.size());
	std::vector<double> x_vs_std (x_vs_.size());

	for (int j=0 ; j<T_.size() ; j++)
	{
		T_std[j] = T_[j];
		x_m_std[j] = x_m_[j];
		x_vs_std[j] = x_vs_[j];
	}


	//setting 1-D particle conditions
	thisParticle.set(
						td.cloud().geometry(),
                        td.cloud().resolution(),
						delta_,
						td.cloud().cylinderLength(),
						td.cloud().rectangleWidth()*td.cloud().rectangleLength(),
						td.cloud().eta_c(),
                        td.cloud().rho_m(),
                        td.cloud().rho_vs(),
                        td.cloud().rho_c(),
                        td.cloud().c_m(),
                        td.cloud().c_vs(),
                        td.cloud().c_c(),
                        td.cloud().k_m(),
                        td.cloud().k_vs(),
                        td.cloud().k_c(),
                        td.cloud().e_m(),
                        td.cloud().e_vs(),
                        td.cloud().e_c(),                        
						T_std,
						x_m_std,
						x_vs_std
					);
    
    thisParticle.stepForward(
    							dt_, 
    							Tc, 
    							mag(Uc-U_), 
    							Gc, 
    							YO2c,
    							status_
    						);

	// updating particle condition from model outputs
    status_ = thisParticle.getState();
    delta_ = thisParticle.getDelta();
    mass_ = thisParticle.getMass();
    vol_ = thisParticle.getVol();

    T_.resize(thisParticle.getNcells());
    x_m_.resize(thisParticle.getNcells());
    x_vs_.resize(thisParticle.getNcells());

	for (unsigned int j=0 ; j<thisParticle.getNcells() ; j++)
	{
		T_[j] = thisParticle.getT()[j];
		x_m_[j] = thisParticle.getXm()[j];
		x_vs_[j] = thisParticle.getXvs()[j];
	}

	Tsurf_ = thisParticle.getT().back();
	volProdRate_tot_ = thisParticle.getMLR();
	volProdRate_fuel_ = thisParticle.getGFRR(); 
	volProdRate_H2O_ = volProdRate_tot_ - volProdRate_fuel_;
	volHRR_ = thisParticle.getCharHRR();
    hConv_ = thisParticle.getHconv(); 
}
 


Foam::vector Foam::biomassParticle::calcVelocity
(
    trackingData& td,
    const word moveParticle,
    const scalar dt_,
    const scalar mass_p,
    const scalar delta_p,
    const scalar rhoc,
    const scalar muc,
    const vector Uc,
    vector& dUTrans
)
{

    // relative velocity
    vector Urel = Uc - U_;

    // Reynolds number
    scalar Re = rhoc*mag(Urel)*(2.0*delta_)/muc;

    scalar Astar = 0.0;
    scalar C_s = 0.0;

    if(td.cloud().geometry() == "rectangle")
    {
        C_s = 1.0;
        Astar = td.cloud().rectangleWidth() * td.cloud().rectangleLength();
        surfToVol_ = 1.0/delta_p; 
        scalar AR = td.cloud().rectangleWidth() / (2.0*delta_p);
        CD_ = 1.98 - 0.8 * (1-Foam::exp(-20.0/AR));
    }
    else if(td.cloud().geometry() == "cylinder")
    {

        C_s = 1.0/constant::mathematical::pi;

        Astar = 2.0*delta_p * td.cloud().cylinderLength();

        surfToVol_ = 2.0/delta_p;

        //safety capping of the drag coefficient
        if(Re <= 0.001)
        {
            CD_ = 1000; 
        }
        else if(Re <= 1.0)
        {
            CD_ = 10.0 / Foam::pow(Re,0.8);
        }
        else if(Re >1.0 && Re <=1000)
        {
        	CD_ = 10*(0.6+0.4*Foam::pow(Re,0.8))/Re;
        }
        else
        {
        	CD_ = 1.0;
        }
    }
    else if(td.cloud().geometry() == "sphere")
    {
        C_s = 1.0/4.0;

        Astar = constant::mathematical::pi * Foam::pow(delta_p,2.0);

		surfToVol_ = 3.0/delta_p;

        if(Re <= 0.001)
        {
            CD_ = 100; 
        }
        else if(Re <= 1.0)
        {
            CD_ = 24.0/Re;
        }
        else if(Re >1.0 && Re <=1000)
        {
            CD_ = 24.0/Re*(0.85+0.15*Foam::pow(Re,0.687));
        }
        else
        {
            CD_ = 0.44;
        }
    }

    if(td.cloud().dragModel() == "constant")
    {
        CD_ = td.cloud().dragCoeff();

    }  

    vector Unew = Zero;
    
    if (moveParticle == "on")
    {
    	scalar B = 0.5*rhoc*CD_*Astar/mass_p * mag(Urel);

        Unew = (U_ + dt_*(B*Uc+ td.g()))/ (1.0+B*dt_);
    }

    dUTrans = 0.5 * rhoc * CD_ * C_s * mag(Uc-Unew)*(Uc-Unew);

    return Unew;
}


// ************************************************************************* //
