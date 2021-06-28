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


inline Foam::scalar Foam::biomassCloud::nParticles() const
{
    return nParticles_;
}

inline Foam::word Foam::biomassCloud::fuel() const
{
    return fuel_;
}


inline Foam::word Foam::biomassCloud::firebrands() const
{
    return firebrands_;
}


inline Foam::scalar Foam::biomassCloud::e() const
{
    return e_;
}


inline Foam::scalar Foam::biomassCloud::mu() const
{
    return mu_;
}

inline Foam::word Foam::biomassCloud::geometry() const
{
    return geometry_;
}

inline Foam::scalar Foam::biomassCloud::resolution() const
{
    return resolution_;
}

inline Foam::scalar Foam::biomassCloud::cylinderLength() const
{
    return cylinderLength_;
}

inline Foam::scalar Foam::biomassCloud::rectangleLength() const
{
    return rectangleLength_;
}

inline Foam::scalar Foam::biomassCloud::rectangleWidth() const
{
    return rectangleWidth_;
}

inline Foam::scalar Foam::biomassCloud::eta_c() const
{
    return eta_c_;
}


inline Foam::scalar Foam::biomassCloud::rho_m() const
{
    return rho_m_;
}
inline Foam::scalar Foam::biomassCloud::rho_vs() const
{
    return rho_vs_;
}
inline Foam::scalar Foam::biomassCloud::rho_c() const
{
    return rho_c_;
}


inline Foam::scalar Foam::biomassCloud::c_m() const
{
    return c_m_;
}
inline Foam::scalar Foam::biomassCloud::c_vs() const
{
    return c_vs_;
}
inline Foam::scalar Foam::biomassCloud::c_c() const
{
    return c_c_;
}


inline Foam::scalar Foam::biomassCloud::k_m() const
{
    return k_m_;
}
inline Foam::scalar Foam::biomassCloud::k_vs() const
{
    return k_vs_;
}
inline Foam::scalar Foam::biomassCloud::k_c() const
{
    return k_c_;
}

inline Foam::scalar Foam::biomassCloud::e_m() const
{
    return e_m_;
}
inline Foam::scalar Foam::biomassCloud::e_vs() const
{
    return e_vs_;
}
inline Foam::scalar Foam::biomassCloud::e_c() const
{
    return e_c_;
}

inline const Foam::SLGThermo& Foam::biomassCloud::thermo() const
{
    return thermo_;
}

//Packing Ratio
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::biomassCloud::packingRatio()
{
    return packingRatio_();
}

inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::biomassCloud::packingRatio() const
{
    return packingRatio_();
}

//Surface Temperature
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::biomassCloud::surfaceTemp()
{
    return surfaceTemp_();
}

inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::biomassCloud::surfaceTemp() const
{
    return surfaceTemp_();
}

//Surface to volume ratio
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::biomassCloud::surfaceToVolumeRatio()
{
    return surfaceToVolumeRatio_();
}

inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::biomassCloud::surfaceToVolumeRatio() const
{
    return surfaceToVolumeRatio_();
}

//Momentum
inline Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::biomassCloud::momentum()
{
    return momentum_();
}

inline const Foam::DimensionedField<Foam::vector, Foam::volMesh>&
Foam::biomassCloud::momentum() const
{
    return momentum_();
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
    tmp<volVectorField::Internal> tUTrans
    (
        new volVectorField::Internal
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
    
    //integrating over cloud dt    
    sourceField = UTrans() / (this->db().time().deltaT());

    return tUTrans;
}


//---------------------------Mass Source-----------------------------
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::biomassCloud::rhoTrans()
{
    return rhoTrans_();
}

inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::biomassCloud::rhoTrans() const
{
    return rhoTrans_();
}

inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::biomassCloud::Srho() const
{
    tmp<volScalarField::Internal> trhoTrans
    (
        new volScalarField::Internal
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
    
    //integrating over cloud dt    
    sourceField = rhoTrans() / (this->db().time().deltaT());

    return trhoTrans;
}


//---------------------------Species Source-----------------------------
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::biomassCloud::rhoYTrans(const label i)
{
    return rhoYTrans_[i];
}


inline const Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::biomassCloud::rhoYTrans() const
{
    return rhoYTrans_;
}


inline Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::biomassCloud::rhoYTrans()
{
    return rhoYTrans_;
}

inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::biomassCloud::SYi(const label i) const
{
    tmp<volScalarField::Internal> trhoYTrans
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + "_rhoYTrans",
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


    scalarField& sourceField = trhoYTrans.ref();
    
    //integrating over cloud dt    
    sourceField = rhoYTrans_[i] / (this->db().time().deltaT());

    return trhoYTrans;
}

//---------------------------Energy Sources-----------------------------
//------- Char oxidation reaction Heat Release (R3)
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::biomassCloud::Qr3Trans()
{
    return Qr3Trans_();
}

inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::biomassCloud::Qr3Trans() const
{
    return Qr3Trans_();
}

inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::biomassCloud::Qr3() const
{
    tmp<volScalarField::Internal> tQr3Trans
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + "_Qr3Trans",
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


    scalarField& sourceField = tQr3Trans.ref();
    
    //integrating over cloud dt    
    sourceField = Qr3Trans() / (this->db().time().deltaT());

    return tQr3Trans;
}

//------- Soild-Gas convectice heat transfer ------------------
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
    tmp<volScalarField::Internal> tQconvTrans
    (
        new volScalarField::Internal
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
    
    //integrating over cloud dt    
    sourceField = QconvTrans() / (this->db().time().deltaT());

    return tQconvTrans;
}



//-------------------- Radiation ---------------------------------
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::biomassCloud::absorption()
{
    return absorption_();
}

inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::biomassCloud::absorption() const
{
    return absorption_();
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
            dimensionedScalar("zero", dimMass/dimLength/pow3(dimTime), 0.0)
        )
    );


    scalarField& sourceField = tEp.ref();
    
    //integrating over cloud dt    
    sourceField = emissionTrans() / (this->db().time().deltaT());

    return tEp;
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

        ap = absorption();


    return tap;
}




// ************************************************************************* //
