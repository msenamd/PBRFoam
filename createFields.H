Info<< "Creating combustion model\n" << endl;

autoPtr<combustionModels::psiCombustionModel> combustion
(
    combustionModels::psiCombustionModel::New
    (
        mesh
    )
);

Info<< "Reading thermophysical properties\n" << endl;

psiReactionThermo& thermo = combustion->thermo();
thermo.validate(args.executable(), "h", "e");

SLGThermo slgThermo(mesh, thermo);

basicMultiComponentMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

scalar  YO2Inf = readScalar(thermo.lookup("YO2Inf"));
scalar  YFInf = readScalar(thermo.lookup("YFInf"));

const word inertSpecie(thermo.lookup("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species()
        << exit(FatalIOError);
}

#include "readAdditionalThermo.H"

Info<< "Creating field rho\n" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volScalarField& p = thermo.p();

Info<< "\nReading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "compressibleCreatePhi.H"

#include "createMRF.H"


Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

// Set the turbulence into the combustion model
combustion->setTurbulence(turbulence());


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"
#include "readpRef.H"

volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

thermo.correct();
rho = thermo.rho();

mesh.setFluxRequired(p_rgh.name());

#include "phrghEqn.H"


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

volScalarField dQ
(
    IOobject
    (
        "dQ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dQ", dimEnergy/dimTime/dimVolume, 0.0)
);

volScalarField dQsoot
(
    IOobject
    (
        "dQsoot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dQsoot", dimEnergy/dimTime/dimVolume, 0.0)  
);


Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));


volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.Cp()*thermo.alpha()
);

volScalarField cp
(
    IOobject
    (
        "cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.Cp()
);

singleStepReactingMixture<gasHThermoPhysics>& singleMixture
(
    dynamic_cast<singleStepReactingMixture<gasHThermoPhysics>&>
    (thermo)
);

// stoichiometric O2 to fuel ratio
scalar s(singleMixture.s().value());

Info <<"s = " << s << endl;

label fuelIndex(singleMixture.fuelIndex());

const label O2Index(composition.species()["O2"]);
const label CO2Index(composition.species()["CO2"]);
const label inertIndex(composition.species()[inertSpecie]);

const volScalarField& O2 = Y[O2Index];  

const volScalarField& fu = Y[fuelIndex];

// Fuel molecular weight
scalar Wu_ = singleMixture.speciesData()[fuelIndex].W();  

// Fuel enthalpy of formation
scalar Hu_ = singleMixture.speciesData()[fuelIndex].hc();  

// Enthalpy of formation of fuel in [J/kg]
dimensionedScalar h0_F("h0_F", dimEnergy/dimMass, Hu_/Wu_);

// Enthalpy of formation of CO in [J/kg]
dimensionedScalar h0_CO("h0_CO", dimEnergy/dimMass, -3945.0e3);

Info << "fuel enthalphy of formation = " << h0_F.value() << endl;
Info << "CO enthalphy of formation = " << h0_CO.value() << endl;

// compute stoichiometric mixture fraction
scalar ftSt = (YO2Inf/s) / ( YFInf + YO2Inf/s );  
Info << "stoichiometric mixture fraction = " << ftSt << nl;  

Info << "specie stoich. coeffs. = " << singleMixture.specieStoichCoeffs() << endl;


// create fileds for mixture fraction
volScalarField ft
(
    IOobject
    (
        "ft",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
(fu-O2/s+YO2Inf/s)/(YFInf+YO2Inf/s)
);

#include "createRadiationModel.H"

volScalarField dQrad
(
    IOobject
    (
        "dQrad",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dQrad", dimEnergy/dimTime/dimVolume, 0.0)
);

volScalarField dQradAbs
(
     IOobject
     (
      "dQradAbs",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedScalar("dQradAbs", dimEnergy/dimTime/dimVolume,0.0)
 );

volScalarField dQradEm
(
     IOobject
     (
      "dQradEm",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedScalar("dQradEm", dimEnergy/dimTime/dimVolume,0.0)
 );

volScalarField G_projected
(
     IOobject
     (
      "G_projected",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedScalar("G_projected", dimEnergy/dimTime/dimArea,0.0)
 );