
// for outputing fire front
OFstream outFireFront("outFireFront_"+runTime.timeName());

outFireFront << "time(s)" << "," << "flameHeight(m)" << "," << "flameFront(m)" 
<< "," << "pyrolysisBack(m)" << "," << "pyrolysisFront(m)" << "," << "smolderingFront(m)" << endl;

volScalarField flameHeight
(
    IOobject
    (
        "flameHeight",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimLength,0.0)
);

volScalarField flameFront
(
    IOobject
    (
        "flameFront",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimLength,0.0)
);

volScalarField pyrolysisBack
(
    IOobject
    (
        "pyrolysisBack",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimLength,0.0)
);

volScalarField pyrolysisFront
(
    IOobject
    (
        "pyrolysisFront",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimLength,0.0)
);

volScalarField smolderingFront
(
    IOobject
    (
        "smolderingFront",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimLength,0.0)
);

// Momentum of the flow [kg m/s]
volVectorField fluidMomentum
(
    IOobject
    (
        "fluidMomentum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimMass*dimLength/dimTime, Zero)
);


// Convective heat flux at boundary surfaces [W/m2]
surfaceScalarField convectiveHeatFlux
(
    IOobject
    (
        "convectiveHeatFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(1, 0, -3, 0 , 0, 0, 0), 0.0)
);


// Radiative heat to the particles per unit bed volume [W/m3]
volScalarField QconvPUVbed
(
    IOobject
    (
        "QconvPUVbed",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
);

// Convective heat flux to the particles per unit bed volume [W/m3]
volScalarField QradPUVbed
(
    IOobject
    (
        "QradPUVbed",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
);


// Radiation absorption/emission by the dispersed particles
volScalarField radAbsPUVbed
(
     IOobject
     (
      "radAbsPUVbed",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedScalar("radAbsPUVbed", dimEnergy/dimTime/dimVolume,0.0)
 );

volScalarField radEmPUVbed
(
     IOobject
     (
      "radEmPUVbed",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedScalar("radEmPUVbed", dimEnergy/dimTime/dimVolume,0.0)
 );