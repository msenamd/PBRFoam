// for outputing flame height
    OFstream outFlameHeightFront("outFlameHeightFront_"+runTime.timeName());
    outFlameHeightFront << "Time(s)" << " " << "FlameHeight(m)" << " " << "FlameFront(m)" << endl;

    volScalarField flameHeight
    (
        IOobject
        (
            "flameHeight",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ
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
            IOobject::NO_READ
        ),
        mesh,
        dimensionedScalar("zero",dimLength,0.0)
    );


// for outputing flame front
    OFstream outPyrolysisFront("outPyrolysisFront_"+runTime.timeName());
    outPyrolysisFront << "Time(s)" << " " << "PyrolysisFront(m)" << endl;

    volScalarField PyrolysisFront
    (
        IOobject
        (
            "PyrolysisFront",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ
        ),
        mesh,
        dimensionedScalar("zero",dimLength,0.0)
    );
    

// Radiative surface heat flux to the particles [W/m2]
    volScalarField QradToFreshParticles
    (
        IOobject
        (
            "QradToFreshParticles",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimEnergy/dimTime/dimArea, 0.0)
    );

// Convective surface heat flux to the particles [W/m2]
    volScalarField QconvToFreshParticles
    (
        IOobject
        (
            "QconvToFreshParticles",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimEnergy/dimTime/dimArea, 0.0)
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
    surfaceScalarField convectiveHeatFlux_L
    (
        IOobject
        (
            "convectiveHeatFlux_L",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(1, 0, -3, 0 , 0, 0, 0), 0.0)
    );

    surfaceScalarField convectiveHeatFlux_T
    (
        IOobject
        (
            "convectiveHeatFlux_T",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimensionSet(1, 0, -3, 0 , 0, 0, 0), 0.0)
    );