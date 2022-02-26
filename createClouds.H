Info<< "\nConstructing solid particle cloud" << endl;
biomassCloud particles
(
    mesh,
    "vegetationBed",
    rho,
    Y[O2Index],
    G_projected,
    U,
    g,
    slgThermo,
    true
);