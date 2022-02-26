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
    - offsetof(biomassParticle, particleState_)
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

    particleState_(true),
    particleSize_(0.0),
    particleVelo_(Zero),    

    particleTemp_(0.0),
    particlePressure_(0.0),
    particleO2MassFraction_(0.0),    
    wetSolidVolFraction_(0.0),
    drySolidVolFraction_(0.0),
    charVolFraction_(0.0),
    ashVolFraction_(0.0),

    surfaceTemp_(0.0),
    surfaceO2MassFrac_(0.0),
    hConv_(0.0),
    CD_(0.0)
{

    DynamicList<scalar> T;
    DynamicList<scalar> p;
    DynamicList<scalar> YO2;
    DynamicList<scalar> x_ws;
    DynamicList<scalar> x_ds;
    DynamicList<scalar> x_c;
    DynamicList<scalar> x_a;

    if (readFields)
    {

       if (is.format() == IOstream::ASCII)
        {
            particleState_ = readLabel(is);
            particleSize_ = readScalar(is);
            is >> particleVelo_;

            is >> T;
            is >> p;
            is >> YO2;
			is >> x_ws;
			is >> x_ds;
            is >> x_c;
            is >> x_a;

        }
        else
        {
            is.read(reinterpret_cast<char*>(&particleState_), sizeofFields_);
            
            is >> T;
            is >> p;
            is >> YO2;            
            is >> x_ws;
            is >> x_ds;
            is >> x_c;
            is >> x_a;

        }

        particleTemp_.transfer(T);
        particlePressure_.transfer(p);
        particleO2MassFraction_.transfer(YO2);
        wetSolidVolFraction_.transfer(x_ws);
        drySolidVolFraction_.transfer(x_ds);
        charVolFraction_.transfer(x_c);
        ashVolFraction_.transfer(x_a);

    }

    // Check state of Istream
    is.check("biomassParticle::biomassParticle(Istream&)");
}

void Foam::biomassParticle::readFields(Cloud<biomassParticle>& c)
{
    if (!c.size())
    {
        return;
    }

    particle::readFields(c);


    IOField<label> particleState
    (
        c.fieldIOobject("particleState", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, particleState);


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


    label i = 0;
    forAllIter(Cloud<biomassParticle>, c, iter)
    {
        biomassParticle& bp = iter();

        bp.particleState_    = particleState[i];
        bp.particleSize_     = particleSize[i];
        bp.particleVelo_     = particleVelo[i];

        bp.particleTemp_             = particleTemp[i];
        bp.particlePressure_         = particlePressure[i];
        bp.particleO2MassFraction_   = particleO2MassFraction[i];
        bp.wetSolidVolFraction_      = wetSolidVolFraction[i];
        bp.drySolidVolFraction_      = drySolidVolFraction[i];
        bp.charVolFraction_          = charVolFraction[i];
        bp.ashVolFraction_           = ashVolFraction[i];

        i++;
    }
}


void Foam::biomassParticle::writeFields(const Cloud<biomassParticle>& c)
{
    particle::writeFields(c);

    label np = c.size();

    IOField<label> particleState(c.fieldIOobject("particleState", IOobject::NO_READ), np);    
    IOField<scalar> particleSize(c.fieldIOobject("particleSize", IOobject::NO_READ), np);
    IOField<vector> particleVelo(c.fieldIOobject("particleVelo", IOobject::NO_READ), np);
    
    IOField<scalarField> particleTemp(c.fieldIOobject("particleTemp", IOobject::NO_READ), np);
    IOField<scalarField> particlePressure(c.fieldIOobject("particlePressure", IOobject::NO_READ), np);
    IOField<scalarField> particleO2MassFraction(c.fieldIOobject("particleO2MassFraction", IOobject::NO_READ), np);
    IOField<scalarField> wetSolidVolFraction(c.fieldIOobject("wetSolidVolFraction", IOobject::NO_READ), np);
    IOField<scalarField> drySolidVolFraction(c.fieldIOobject("drySolidVolFraction", IOobject::NO_READ), np);
    IOField<scalarField> charVolFraction(c.fieldIOobject("charVolFraction", IOobject::NO_READ), np);
    IOField<scalarField> ashVolFraction(c.fieldIOobject("ashVolFraction", IOobject::NO_READ), np);

    IOField<scalar> surfaceTemp(c.fieldIOobject("surfaceTemp", IOobject::NO_READ), np);
    IOField<scalar> surfaceO2MassFrac(c.fieldIOobject("surfaceO2MassFrac", IOobject::NO_READ), np);
    IOField<scalar> hConv(c.fieldIOobject("hConv", IOobject::NO_READ), np);
    IOField<scalar> CD(c.fieldIOobject("CD", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(Cloud<biomassParticle>, c, iter)
    {
        const biomassParticle& bp = iter();

        particleState[i] = bp.particleState_;
        particleSize[i]  = bp.particleSize_;
        particleVelo[i]  = bp.particleVelo_;

        particleTemp[i]           = bp.particleTemp_;
        particlePressure[i]       = bp.particlePressure_;
        particleO2MassFraction[i] = bp.particleO2MassFraction_;
        wetSolidVolFraction[i]    = bp.wetSolidVolFraction_;
        drySolidVolFraction[i]    = bp.drySolidVolFraction_;
        charVolFraction[i]        = bp.charVolFraction_;
        ashVolFraction[i]         = bp.ashVolFraction_;

        surfaceTemp[i]       = bp.surfaceTemp_;
        surfaceO2MassFrac[i] = bp.surfaceO2MassFrac_;
        hConv[i]             = bp.hConv_;
        CD[i]                = bp.CD_;

        i++;
    }

    particleState.write();
    particleSize.write();
    particleVelo.write();
    
    particleTemp.write();
    particlePressure.write();
    particleO2MassFraction.write();
    wetSolidVolFraction.write();
    drySolidVolFraction.write();
    charVolFraction.write();
    ashVolFraction.write();

    surfaceTemp.write();
    surfaceO2MassFrac.write();
    hConv.write();
    CD.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const biomassParticle& bp)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(bp)
            << token::SPACE << bp.particleState_
            << token::SPACE << bp.particleSize_
            << token::SPACE << bp.particleVelo_
            << token::SPACE << bp.particleTemp_
            << token::SPACE << bp.particlePressure_
            << token::SPACE << bp.particleO2MassFraction_
            << token::SPACE << bp.wetSolidVolFraction_
            << token::SPACE << bp.drySolidVolFraction_
            << token::SPACE << bp.charVolFraction_
            << token::SPACE << bp.ashVolFraction_;
    }
    else
    {
        os  << static_cast<const particle&>(bp);
        os.write
        (
            reinterpret_cast<const char*>(&bp.particleState_),
            biomassParticle::sizeofFields_
        );
        os  << bp.particleTemp_;
        os  << bp.particlePressure_;
        os  << bp.particleO2MassFraction_;
        os  << bp.wetSolidVolFraction_;
        os  << bp.drySolidVolFraction_;
        os  << bp.charVolFraction_;
        os  << bp.ashVolFraction_;
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const biomassParticle&)");

    return os;
}


// ************************************************************************* //
