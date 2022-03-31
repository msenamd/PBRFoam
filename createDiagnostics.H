
// for outputing fire front
    OFstream outFireFront("outFireFront_"+runTime.timeName());

    outFireFront << "time(s)" << "," << "flameHeight(m)" << "," << "flameFront(m)" 
    << "," << "dryingFront(m)" << "," << "pyrolysisFront(m)" << "," << "smolderingFront(m)" << endl;

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

    volScalarField dryingFront
    (
        IOobject
        (
            "dryingFront",
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


// Write particle properties if their IDs are > 0
forAllConstIter(biomassCloud, particles, iter)
{
    if(iter().particleID() > 0)
    {
        std::ofstream outParticle("particle" + name(iter().particleID()) + ".csv", ios::app);

        outParticle << "time(s)" << "," << "x(m)" << "," << "y(m)" << "," << "z(m)" << "," 
                    << "state" << "," << "dt(s)" << "," << "size(m)" << ","  
                    << "Ux(m/s)" << "," << "Uy(m/s)" << "," << "Uz(m/s)" << ","
                    << "mass(kg)" << "," << "surfaceTemp(K)"  << "," << "coreTemp(K)"  << ","
                    << "convFlux(W/m2)" << "," << "radFlux(W/m2)" << "," << "massFlux(kg/s/m2)" << ","
                    << "surfaceO2MassFrac(-)" << "," << "hConv(W/m2/K)" << "," << "CD(-)" << ","
                    << "dryingRate(kg/s)" << "," << "pyrolysisRate(kg/s)" << "," << "oxidPyrolysisRate(kg/s)" << ","
                    << "charOxidRate(kg/s)" << "," << "massLossRate(kg/s)" << endl;
    }
}