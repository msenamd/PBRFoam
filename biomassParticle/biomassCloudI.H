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


//---------------------------Mass Source-----------------------------

inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::biomassCloud::massTrans()
{
    return massTrans_();
}

inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::biomassCloud::massTrans() const
{
    return massTrans_();
}

inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::biomassCloud::Srho() const
{
    tmp<DimensionedField<scalar, volMesh>> tmassTrans
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + "_massTrans",
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

    scalarField& sourceField = tmassTrans.ref();
    
    sourceField = massTrans() / (this->db().time().deltaT());

    return tmassTrans;
}


//---------------------------Species Source-----------------------------

inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::biomassCloud::speciesMassTrans(const label i)
{
    return speciesMassTrans_[i];
}

inline const Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::biomassCloud::speciesMassTrans() const
{
    return speciesMassTrans_;
}

inline Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::biomassCloud::speciesMassTrans()
{
    return speciesMassTrans_;
}

inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::biomassCloud::SYi(const label i) const
{
    tmp<DimensionedField<scalar, volMesh>> tspeciesMassTrans
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                this->name() + "_speciesMassTrans",
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

    scalarField& sourceField = tspeciesMassTrans.ref();
    
    sourceField = speciesMassTrans()[i] / (this->db().time().deltaT());

    return tspeciesMassTrans;
}


//---------------------------Momentum Source-----------------------------

inline Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::biomassCloud::momentumTrans()
{
    return momentumTrans_();
}

inline const Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::biomassCloud::momentumTrans() const
{
    return momentumTrans_();
}

inline Foam::tmp<Foam::DimensionedField<Foam::vector, Foam::volMesh>>
Foam::biomassCloud::SU() const
{
    tmp<DimensionedField<vector, volMesh>> tmomentumTrans
    (
        new DimensionedField<vector, volMesh>
        (
            IOobject
            (
                this->name() + "_momentumTrans",
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

    vectorField& sourceField = tmomentumTrans.ref();
    
    sourceField = momentumTrans() / (this->db().time().deltaT());

    return tmomentumTrans;
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
    
    sourceField = QconvTrans() / (this->db().time().deltaT());

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

        ap = absorptionCoeff();

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
    
    sourceField = emissionTrans() / (this->db().time().deltaT());

    return tEp;
}


// ************************************************************************* //
