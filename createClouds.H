Info<< "\nConstructing solid particle cloud" << endl;
biomassCloud particles
(
    mesh,
    "vegetationBed",
    rho,
    U,
    g,
    slgThermo,
    true
);