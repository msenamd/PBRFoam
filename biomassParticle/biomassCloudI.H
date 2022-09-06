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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::fvMesh& Foam::biomassCloud::mesh() const
{
    return mesh_;
}


inline const Foam::SLGThermo& Foam::biomassCloud::thermo() const
{
    return thermo_;
}


//---------------------------Species Source-----------------------------

inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::biomassCloud::rhoTrans(const label i)
{
    return rhoTrans_[i];
}

inline const Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::biomassCloud::rhoTrans() const
{
    return rhoTrans_;
}

inline Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::biomassCloud::rhoTrans()
{
    return rhoTrans_;
}

inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::biomassCloud::SYi(const label i) const
{
    tmp<DimensionedField<scalar, volMesh>> trhoTrans
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + "_rhoTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar
            (
                "zero",
                dimMass/dimVolume/dimTime,
                0.0
            )
        )
    );

    scalarField& sourceField = trhoTrans.ref();
    
    sourceField = rhoTrans_[i] / this->db().time().deltaTValue();

    return trhoTrans;
}


//---------------------------Mass Source-----------------------------

inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::biomassCloud::Srho() const
{
    tmp<DimensionedField<scalar, volMesh>> trhoTrans
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + "_rhoTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar
            (
                "zero",
                dimMass/dimVolume/dimTime,
                0.0
            )
        )
    );

    scalarField& sourceField = trhoTrans.ref();
    
    forAll(rhoTrans_, i)
    {
        sourceField += rhoTrans_[i] / this->db().time().deltaTValue();
    }

    return trhoTrans;
}

//---------------------------Momentum Source-----------------------------

inline Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::biomassCloud::UTrans()
{
    return UTrans_();
}

inline const Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::biomassCloud::UTrans() const
{
    return UTrans_();
}

inline Foam::tmp<Foam::DimensionedField<Foam::vector, Foam::volMesh>>
Foam::biomassCloud::SU() const
{
    tmp<DimensionedField<vector, volMesh>> tUTrans
    (
        new DimensionedField<vector, volMesh>
        (
            IOobject
            (
                this->name() + "_UTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedVector
            (
                "zero",
                dimForce/dimVolume,
                Zero
            )
        )
    );

    vectorField& sourceField = tUTrans.ref();
    
    sourceField = UTrans() / (this->db().time().deltaTValue());

    return tUTrans;
}


//---------------------------Energy Sources-----------------------------

inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::biomassCloud::QconvTrans()
{
    return QconvTrans_();
}

inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::biomassCloud::QconvTrans() const
{
    return QconvTrans_();
}

inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::biomassCloud::Qconv() const
{
    tmp<DimensionedField<scalar, volMesh>> tQconvTrans
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + "_QconvTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar
            (
                "zero",
                dimEnergy/dimVolume/dimTime,
                0.0
            )
        )
    );


    scalarField& sourceField = tQconvTrans.ref();
    
    sourceField = QconvTrans() / (this->db().time().deltaTValue());

    return tQconvTrans;
}



//-------------------- Radiation ---------------------------------
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::biomassCloud::absorptionCoeff()
{
    return absorptionCoeff_();
}

inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::biomassCloud::absorptionCoeff() const
{
    return absorptionCoeff_();
}

inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::biomassCloud::emissionTrans()
{
    return emissionTrans_();
}

inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::biomassCloud::emissionTrans() const
{
    return emissionTrans_();
}

inline Foam::tmp<Foam::volScalarField> Foam::biomassCloud::ap() const
{
    tmp<volScalarField> tap
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "_ap",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );

        scalarField& ap = tap.ref().primitiveFieldRef();

        ap = absorptionCoeff() / (this->db().time().deltaTValue());

    return tap;
}

inline Foam::tmp<Foam::volScalarField> Foam::biomassCloud::Ep() const
{
    tmp<volScalarField> tEp
    (
        new volScalarField
        (
            IOobject
            (
                this->name() + "_Ep",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );

    scalarField& sourceField = tEp.ref().primitiveFieldRef();
    
    sourceField = emissionTrans() / (this->db().time().deltaTValue());

    return tEp;
}


// ************************************************************************* //
