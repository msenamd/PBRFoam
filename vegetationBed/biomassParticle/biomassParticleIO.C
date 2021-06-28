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
    offsetof(biomassParticle, T_)
    - offsetof(biomassParticle, delta_)
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
    delta_(0),
    U_(Zero),
    status_(false),
    T_(0),
    x_m_(0),
    x_vs_(0),
    
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

    DynamicList<scalar> TT;
    DynamicList<scalar> XM;
    DynamicList<scalar> XVS;

    if (readFields)
    {

       if (is.format() == IOstream::ASCII)
        {
            delta_ = readScalar(is);
            is >> U_;
            status_ = readBool(is);
            
            is >> TT;
			is >> XM;
			is >> XVS;

        }
        else
        {
            is.read(reinterpret_cast<char*>(&delta_), sizeofFields_);
            is >> TT;
            is >> XM;
            is >> XVS;
        }

        T_.transfer(TT);
        x_m_.transfer(XM);
        x_vs_.transfer(XVS);
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

    IOField<scalar> delta
    	(
    		c.fieldIOobject("delta", IOobject::MUST_READ)
    	);
    c.checkFieldIOobject(c, delta);

    IOField<vector> U
    	(
    		c.fieldIOobject("U", IOobject::MUST_READ)
    	);
    c.checkFieldIOobject(c, U);

    IOField<label> status
    (
        c.fieldIOobject("status", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, status);

    IOField<scalarField> T
    	(
    		c.fieldIOobject("T", IOobject::MUST_READ)
    	);
    c.checkFieldIOobject(c, T);

    IOField<scalarField> x_m
    	(
    		c.fieldIOobject("x_m", IOobject::MUST_READ)
    	);
    c.checkFieldIOobject(c, x_m);


    IOField<scalarField> x_vs
    	(
    		c.fieldIOobject("x_vs", IOobject::MUST_READ)
    	);
    c.checkFieldIOobject(c, x_vs);


    label i = 0;
    forAllIter(Cloud<biomassParticle>, c, iter)
    {
        biomassParticle& p = iter();

        p.delta_ = delta[i];
        p.U_ = U[i];
        p.status_ = status[i];
        p.T_ = T[i];
        p.x_m_ = x_m[i];
        p.x_vs_ = x_vs[i];
        i++;
    }
}


void Foam::biomassParticle::writeFields(const Cloud<biomassParticle>& c)
{
    particle::writeFields(c);

    label np = c.size();

    IOField<scalar> delta(c.fieldIOobject("delta", IOobject::NO_READ), np);
    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    
    IOField<label> status(c.fieldIOobject("status", IOobject::NO_READ), np);    
    IOField<scalarField> T(c.fieldIOobject("T", IOobject::NO_READ), np);
    IOField<scalarField> x_m(c.fieldIOobject("x_m", IOobject::NO_READ), np);
    IOField<scalarField> x_vs(c.fieldIOobject("x_vs", IOobject::NO_READ), np);

    IOField<scalar> mass(c.fieldIOobject("mass", IOobject::NO_READ), np);
    IOField<scalar> volume(c.fieldIOobject("volume", IOobject::NO_READ), np);

    IOField<scalar> Tsurf(c.fieldIOobject("Tsurf", IOobject::NO_READ), np);
    IOField<scalar> volProdRate_tot(c.fieldIOobject("volProdRate_tot", IOobject::NO_READ), np);
    IOField<scalar> volProdRate_fuel(c.fieldIOobject("volProdRate_fuel", IOobject::NO_READ), np);
    IOField<scalar> volProdRate_H2O(c.fieldIOobject("volProdRate_H2O", IOobject::NO_READ), np);
    IOField<scalar> volHRR(c.fieldIOobject("volHRR_", IOobject::NO_READ), np);
    IOField<scalar> hConv(c.fieldIOobject("hConv", IOobject::NO_READ), np);
    IOField<scalar> CD(c.fieldIOobject("CD", IOobject::NO_READ), np);


    label i = 0;
    forAllConstIter(Cloud<biomassParticle>, c, iter)
    {
        const biomassParticle& p = iter();

        delta[i] = p.delta_;
        U[i] = p.U_;

        status[i] = p.status_;
        T[i] = p.T_;
        x_m[i] = p.x_m_;
        x_vs[i] = p.x_vs_;

        mass[i] = p.mass_;
        volume[i] = p.vol_;
        Tsurf[i] = p.Tsurf_;
        volProdRate_tot[i] = p.volProdRate_tot_;
        volProdRate_fuel[i] = p.volProdRate_fuel_;
        volProdRate_H2O[i] = p.volProdRate_H2O_;
        volHRR[i] = p.volHRR_;
        hConv[i] = p.hConv_;
        CD[i] = p.CD_;

        i++;
    }

    
    delta.write();
    U.write();
    
    status.write();
    T.write();
    x_m.write();
    x_vs.write();

    mass.write();
    volume.write();
    Tsurf.write();
    volProdRate_tot.write();
    volProdRate_fuel.write();
    volProdRate_H2O.write();
    volHRR.write();
    hConv.write();
    CD.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const biomassParticle& p)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.delta()
            << token::SPACE << p.U()
            << token::SPACE << p.status()
            << token::SPACE << p.T()
            << token::SPACE << p.x_m()
            << token::SPACE << p.x_vs();
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.delta_),
            biomassParticle::sizeofFields_
        );
        os  << p.T();
        os  << p.x_m();
        os  << p.x_vs();
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const biomassParticle&)");

    return os;
}


// ************************************************************************* //
