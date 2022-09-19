/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "biomassParticle.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const std::size_t Foam::biomassParticle::sizeofFields_
(
    offsetof(biomassParticle, particleTemp_)
    - offsetof(biomassParticle, particleID_)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::biomassParticle::biomassParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    particle(mesh, is, readFields),

    particleID_(0),
    particleState_(0),
    nParticlesPerSuperParticle_(0),
    particledt_(0.0),
    particleSize_(0.0),
    particleVelo_(Zero),    

    particleTemp_(0.0),
    particlePressure_(0.0),
    particleO2MassFraction_(0.0),    
    wetSolidVolFraction_(0.0),
    drySolidVolFraction_(0.0),
    charVolFraction_(0.0),
    ashVolFraction_(0.0),
    integralWetSolidMass_(0.0),
    integralDrySolidMass_(0.0),
    integralCharMass_(0.0),

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

    DynamicList<scalar> T;
    DynamicList<scalar> p;
    DynamicList<scalar> YO2;
    DynamicList<scalar> x_ws;
    DynamicList<scalar> x_ds;
    DynamicList<scalar> x_c;
    DynamicList<scalar> x_a;
    DynamicList<scalar> m_ws;
    DynamicList<scalar> m_ds;
    DynamicList<scalar> m_c;


    if (readFields)
    {

       if (is.format() == IOstream::ASCII)
        {
            particleID_ = readLabel(is);
            particledt_ = readScalar(is);
            particleSize_ = readScalar(is);
            is >> particleVelo_;

            is >> T;
            is >> p;
            is >> YO2;
			is >> x_ws;
			is >> x_ds;
            is >> x_c;
            is >> x_a;
            is >> m_ws;
            is >> m_ds;
            is >> m_c;

        }
        else
        {
            is.read(reinterpret_cast<char*>(&particleID_), sizeofFields_);
            
            is >> T;
            is >> p;
            is >> YO2;            
            is >> x_ws;
            is >> x_ds;
            is >> x_c;
            is >> x_a;
            is >> m_ws;
            is >> m_ds;
            is >> m_c;
        }

        particleTemp_.transfer(T);
        particlePressure_.transfer(p);
        particleO2MassFraction_.transfer(YO2);
        wetSolidVolFraction_.transfer(x_ws);
        drySolidVolFraction_.transfer(x_ds);
        charVolFraction_.transfer(x_c);
        ashVolFraction_.transfer(x_a);
        integralWetSolidMass_.transfer(m_ws);
        integralDrySolidMass_.transfer(m_ds);
        integralCharMass_.transfer(m_c);
    }

    // Check state of Istream
    is.check("biomassParticle::biomassParticle(Istream&)");

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
}

void Foam::biomassParticle::readFields(Cloud<biomassParticle>& c)
{
    if (!c.size())
    {
        return;
    }

    particle::readFields(c);

    IOField<label> particleID
    (
        c.fieldIOobject("particleID", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, particleID);


    IOField<label> particleState
    (
        c.fieldIOobject("particleState", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, particleState);

    IOField<label> nParticlesPerSuperParticle
    (
        c.fieldIOobject("nParticlesPerSuperParticle", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, nParticlesPerSuperParticle);

    IOField<scalar> particledt
        (
            c.fieldIOobject("particledt", IOobject::MUST_READ)
        );
    c.checkFieldIOobject(c, particledt);

    IOField<scalar> particleSize
    	(
    		c.fieldIOobject("particleSize", IOobject::MUST_READ)
    	);
    c.checkFieldIOobject(c, particleSize);


    IOField<vector> particleVelo
    	(
    		c.fieldIOobject("particleVelo", IOobject::MUST_READ)
    	);
    c.checkFieldIOobject(c, particleVelo);


    IOField<scalarField> particleTemp
    	(
    		c.fieldIOobject("particleTemp", IOobject::MUST_READ)
    	);
    c.checkFieldIOobject(c, particleTemp);


    IOField<scalarField> particlePressure
        (
            c.fieldIOobject("particlePressure", IOobject::MUST_READ)
        );
    c.checkFieldIOobject(c, particlePressure);


    IOField<scalarField> particleO2MassFraction
        (
            c.fieldIOobject("particleO2MassFraction", IOobject::MUST_READ)
        );
    c.checkFieldIOobject(c, particleO2MassFraction);


    IOField<scalarField> wetSolidVolFraction
    	(
    		c.fieldIOobject("wetSolidVolFraction", IOobject::MUST_READ)
    	);
    c.checkFieldIOobject(c, wetSolidVolFraction);


    IOField<scalarField> drySolidVolFraction
        (
            c.fieldIOobject("drySolidVolFraction", IOobject::MUST_READ)
        );
    c.checkFieldIOobject(c, drySolidVolFraction);


    IOField<scalarField> charVolFraction
    (
        c.fieldIOobject("charVolFraction", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, charVolFraction);


    IOField<scalarField> ashVolFraction
    (
        c.fieldIOobject("ashVolFraction", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, ashVolFraction);    

    IOField<scalarField> integralWetSolidMass
    (
        c.fieldIOobject("integralWetSolidMass", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, integralWetSolidMass);   

    IOField<scalarField> integralDrySolidMass
    (
        c.fieldIOobject("integralDrySolidMass", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, integralDrySolidMass);  

    IOField<scalarField> integralCharMass
    (
        c.fieldIOobject("integralCharMass", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, integralCharMass);  


    label i = 0;
    forAllIter(Cloud<biomassParticle>, c, iter)
    {
        biomassParticle& bp = iter();

        bp.particleID_              = particleID[i];
        bp.particleState_           = particleState[i];
        bp.nParticlesPerSuperParticle_ = nParticlesPerSuperParticle[i];
        bp.particledt_              = particledt[i];
        bp.particleSize_            = particleSize[i];
        bp.particleVelo_            = particleVelo[i];

        bp.particleTemp_            = particleTemp[i];
        bp.particlePressure_        = particlePressure[i];
        bp.particleO2MassFraction_  = particleO2MassFraction[i];
        bp.wetSolidVolFraction_     = wetSolidVolFraction[i];
        bp.drySolidVolFraction_     = drySolidVolFraction[i];
        bp.charVolFraction_         = charVolFraction[i];
        bp.ashVolFraction_          = ashVolFraction[i];
        bp.integralWetSolidMass_    = integralWetSolidMass[i];
        bp.integralDrySolidMass_    = integralDrySolidMass[i];
        bp.integralCharMass_        = integralCharMass[i];

        if(bp.particleID_ > 0)
        {
            std::ofstream outParticle(bp.outputPath/ "particle" + name(bp.particleID_) + ".csv", ios::app);

            outParticle << "time(s)" << "," << "x(m)" << "," << "y(m)" << "," << "z(m)" << "," 
                        << "state" << "," << "dt(s)" << "," << "size(m)" << "," << "surfToVol(1/m)" << ","  
                        << "Ux(m/s)" << "," << "Uy(m/s)" << "," << "Uz(m/s)" << ","
                        << "mass(kg)" << "," << "surfaceTemp(K)"  << "," << "coreTemp(K)"  << ","
                        << "convFlux(W/m2)" << "," << "radFlux(W/m2)" << "," << "massFlux(kg/s/m2)" << ","
                        << "surfaceO2MassFrac(-)" << "," << "hConv(W/m2/K)" << "," << "CD(-)" << ","
                        << "dryingRate(kg/s)" << "," << "pyrolysisRate(kg/s)" << "," << "oxidPyrolysisRate(kg/s)" << ","
                        << "charOxidRate(kg/s)" << "," << "massLossRate(kg/s)" << "," << "heatReleaseRate(kg/s)" << "\n";        
        }
        i++;
    }


}


void Foam::biomassParticle::writeFields(const Cloud<biomassParticle>& c)
{
    particle::writeFields(c);

    label np = c.size();

    IOField<label> particleID(c.fieldIOobject("particleID", IOobject::NO_READ), np);    
    IOField<label> particleState(c.fieldIOobject("particleState", IOobject::NO_READ), np);    
    IOField<label> nParticlesPerSuperParticle(c.fieldIOobject("nParticlesPerSuperParticle", IOobject::NO_READ), np);    
    IOField<scalar> particledt(c.fieldIOobject("particledt", IOobject::NO_READ), np);
    IOField<scalar> particleSize(c.fieldIOobject("particleSize", IOobject::NO_READ), np);
    IOField<vector> particleVelo(c.fieldIOobject("particleVelo", IOobject::NO_READ), np);
    
    IOField<scalarField> particleTemp(c.fieldIOobject("particleTemp", IOobject::NO_READ), np);
    IOField<scalarField> particlePressure(c.fieldIOobject("particlePressure", IOobject::NO_READ), np);
    IOField<scalarField> particleO2MassFraction(c.fieldIOobject("particleO2MassFraction", IOobject::NO_READ), np);
    IOField<scalarField> wetSolidVolFraction(c.fieldIOobject("wetSolidVolFraction", IOobject::NO_READ), np);
    IOField<scalarField> drySolidVolFraction(c.fieldIOobject("drySolidVolFraction", IOobject::NO_READ), np);
    IOField<scalarField> charVolFraction(c.fieldIOobject("charVolFraction", IOobject::NO_READ), np);
    IOField<scalarField> ashVolFraction(c.fieldIOobject("ashVolFraction", IOobject::NO_READ), np);
    IOField<scalarField> integralWetSolidMass(c.fieldIOobject("integralWetSolidMass", IOobject::NO_READ), np);
    IOField<scalarField> integralDrySolidMass(c.fieldIOobject("integralDrySolidMass", IOobject::NO_READ), np);
    IOField<scalarField> integralCharMass(c.fieldIOobject("integralCharMass", IOobject::NO_READ), np);

    IOField<scalar> particleMass(c.fieldIOobject("particleMass", IOobject::NO_READ), np);
    IOField<scalar> particleVol(c.fieldIOobject("particleVol", IOobject::NO_READ), np);
    IOField<scalar> surfaceTemp(c.fieldIOobject("surfaceTemp", IOobject::NO_READ), np);
    IOField<scalar> coreTemp(c.fieldIOobject("coreTemp", IOobject::NO_READ), np);    
    IOField<scalar> convFlux(c.fieldIOobject("convFlux", IOobject::NO_READ), np);
    IOField<scalar> radFlux(c.fieldIOobject("radFlux", IOobject::NO_READ), np);
    IOField<scalar> massFlux(c.fieldIOobject("massFlux", IOobject::NO_READ), np);
    IOField<scalar> surfaceO2MassFrac(c.fieldIOobject("surfaceO2MassFrac", IOobject::NO_READ), np);
    IOField<scalar> hConv(c.fieldIOobject("hConv", IOobject::NO_READ), np);
    IOField<scalar> CD(c.fieldIOobject("CD", IOobject::NO_READ), np);

    IOField<scalar> dryingRate(c.fieldIOobject("dryingRate", IOobject::NO_READ), np);
    IOField<scalar> pyrolysisRate(c.fieldIOobject("pyrolysisRate", IOobject::NO_READ), np);
    IOField<scalar> oxidPyrolysisRate(c.fieldIOobject("oxidPyrolysisRate", IOobject::NO_READ), np);
    IOField<scalar> charOxidRate(c.fieldIOobject("charOxidRate", IOobject::NO_READ), np);
    IOField<scalar> massLossRate(c.fieldIOobject("massLossRate", IOobject::NO_READ), np);
    IOField<scalar> heatReleaseRate(c.fieldIOobject("heatReleaseRate", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(Cloud<biomassParticle>, c, iter)
    {
        const biomassParticle& bp = iter();

        particleID[i]             = bp.particleID_;
        particleState[i]          = bp.particleState_;
        nParticlesPerSuperParticle[i] = bp.nParticlesPerSuperParticle_;
        particledt[i]             = bp.particledt_;
        particleSize[i]           = bp.particleSize_;
        particleVelo[i]           = bp.particleVelo_;

        particleTemp[i]           = bp.particleTemp_;
        particlePressure[i]       = bp.particlePressure_;
        particleO2MassFraction[i] = bp.particleO2MassFraction_;
        wetSolidVolFraction[i]    = bp.wetSolidVolFraction_;
        drySolidVolFraction[i]    = bp.drySolidVolFraction_;
        charVolFraction[i]        = bp.charVolFraction_;
        ashVolFraction[i]         = bp.ashVolFraction_;
        integralWetSolidMass[i]   = bp.integralWetSolidMass_;
        integralDrySolidMass[i]   = bp.integralDrySolidMass_;
        integralCharMass[i]       = bp.integralCharMass_;

        particleMass[i]           = bp.particleMass_;
        particleVol[i]            = bp.particleVol_;
        surfaceTemp[i]            = bp.surfaceTemp_;
        coreTemp[i]               = bp.coreTemp_;
        convFlux[i]               = bp.convFlux_;
        radFlux[i]                = bp.radFlux_;
        massFlux[i]               = bp.massFlux_;
        surfaceO2MassFrac[i]      = bp.surfaceO2MassFrac_;
        hConv[i]                  = bp.hConv_;
        CD[i]                     = bp.CD_;

        dryingRate[i]             = bp.dryingRate_;
        pyrolysisRate[i]          = bp.pyrolysisRate_;
        oxidPyrolysisRate[i]      = bp.oxidPyrolysisRate_;
        charOxidRate[i]           = bp.charOxidRate_;
        massLossRate[i]           = bp.massLossRate_;
        heatReleaseRate[i]        = bp.heatReleaseRate_;

        i++;
    }

    particleID.write();
    particleState.write();
    nParticlesPerSuperParticle.write();
    particledt.write();
    particleSize.write();
    particleVelo.write();
    
    particleTemp.write();
    particlePressure.write();
    particleO2MassFraction.write();
    wetSolidVolFraction.write();
    drySolidVolFraction.write();
    charVolFraction.write();
    ashVolFraction.write();
    integralWetSolidMass.write();
    integralDrySolidMass.write();
    integralCharMass.write();

    particleMass.write();
    particleVol.write();
    surfaceTemp.write();
    coreTemp.write();
    convFlux.write();
    radFlux.write();
    massFlux.write();
    surfaceO2MassFrac.write();
    hConv.write();
    CD.write();

    dryingRate.write();
    pyrolysisRate.write();
    oxidPyrolysisRate.write();
    charOxidRate.write();
    massLossRate.write();
    heatReleaseRate.write();

}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const biomassParticle& bp)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(bp)
            << token::SPACE << bp.particleID_
            << token::SPACE << bp.particleState_
            << token::SPACE << bp.nParticlesPerSuperParticle_
            << token::SPACE << bp.particledt_
            << token::SPACE << bp.particleSize_
            << token::SPACE << bp.particleVelo_
            << token::SPACE << bp.particleTemp_
            << token::SPACE << bp.particlePressure_
            << token::SPACE << bp.particleO2MassFraction_
            << token::SPACE << bp.wetSolidVolFraction_
            << token::SPACE << bp.drySolidVolFraction_
            << token::SPACE << bp.charVolFraction_
            << token::SPACE << bp.ashVolFraction_
            << token::SPACE << bp.integralWetSolidMass_
            << token::SPACE << bp.integralDrySolidMass_
            << token::SPACE << bp.integralCharMass_;
    }
    else
    {
        os  << static_cast<const particle&>(bp);
        os.write
        (
            reinterpret_cast<const char*>(&bp.particleID_),
            biomassParticle::sizeofFields_
        );
        os  << bp.particleTemp_;
        os  << bp.particlePressure_;
        os  << bp.particleO2MassFraction_;
        os  << bp.wetSolidVolFraction_;
        os  << bp.drySolidVolFraction_;
        os  << bp.charVolFraction_;
        os  << bp.ashVolFraction_;
        os  << bp.integralWetSolidMass_;
        os  << bp.integralDrySolidMass_;
        os  << bp.integralCharMass_;
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const biomassParticle&)");

    return os;
}


// ************************************************************************* //
