
// for outputing fire front
    OFstream outFireFront("outFireFront_"+runTime.timeName());

    outFireFront << "Time(s)" << " " << "FlameHeight(m)" << " " << "FlameFront(m)" 
    << " " << "PyrolysisFront(m)" << " " << "SmolderingFront(m)" << endl;

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

    volScalarField PyrolysisFront
    (
        IOobject
        (
            "PyrolysisFront",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimLength,0.0)
    );

    volScalarField SmolderingFront
    (
        IOobject
        (
            "SmolderingFront",
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

// Vegetation bed mass rates per unit bed volume [W/m2/m3]
    volScalarField mFuelPUVbed
    (
        IOobject
        (
            "mFuelPUVbed",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/dimTime/dimVolume, 0.0)
    );

    volScalarField RR1PUVbed
    (
        IOobject
        (
            "RR1PUVbed",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/dimTime/dimVolume, 0.0)
    );

    volScalarField RR2PUVbed
    (
        IOobject
        (
            "RR2PUVbed",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/dimTime/dimVolume, 0.0)
    );
    volScalarField RR3PUVbed
    (
        IOobject
        (
            "RR3PUVbed",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/dimTime/dimVolume, 0.0)
    );
    volScalarField RR4PUVbed
    (
        IOobject
        (
            "RR4PUVbed",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/dimTime/dimVolume, 0.0)
    );
