/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.
    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
Application
    setParticles -- written by Mohamed Ahmed
Description
    outputs locations of lagrangian particles
    Reads in initLOVDict.
\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
#include "IOstreams.H"
#include "fvMesh.H"

void writeHeader(OFstream& os, word className, word objectName)
{
    os << "FoamFile" << nl;
    os << "{" << nl;
    os << "     version     2.0;" << nl;
    os << "     format      ascii;" << nl;
    os << "     class       "<< className << ";" << nl;
    os << "     location    0;" << nl;
    os << "     object      "<< objectName << ";" << nl;
    os << "}" << nl;
    os << nl;
}

template <typename objectType>
void writeData(
                OFstream& os, fvMesh& mesh, objectType init, int numSuperParticles, 
                scalar xMin, scalar yMin, scalar zMin,
                scalar xMax, scalar yMax, scalar zMax
                )
{
    os  << numSuperParticles << nl;
    os  << '(' << nl;
    forAll(mesh.C(), celli)
    {
        if( mesh.C()[celli][0] >= xMin && mesh.C()[celli][1] >= yMin && mesh.C()[celli][2] >= zMin
            && mesh.C()[celli][0] <= xMax && mesh.C()[celli][1] <= yMax && mesh.C()[celli][2] <= zMax
            )
        {

            os  << init << nl;
            
        }
    }
    os  << ')' << nl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

    OFstream os_pos(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"positions");
    OFstream os_state(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"particleState");
    OFstream os_dt(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"particledt");
    OFstream os_size(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"particleSize");
    OFstream os_velo(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"particleVelo");
    OFstream os_T(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"particleTemp");
    OFstream os_p(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"particlePressure");
    OFstream os_O2(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"particleO2MassFraction");
    OFstream os_wetSolidVolFrac(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"wetSolidVolFraction");
    OFstream os_drySolidVolFrac(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"drySolidVolFraction");
    OFstream os_charVolFrac(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"charVolFraction");
    OFstream os_ashVolFrac(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"ashVolFraction");
    OFstream os_wetSolidMass(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"integralWetSolidMass");
    OFstream os_drySolidMass(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"integralDrySolidMass");
    OFstream os_charMass(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"integralCharMass");


    writeHeader(os_pos, "Cloud<solidParticle>", "positions");
    writeHeader(os_state, "labelField", "particleState");
    writeHeader(os_dt, "scalarField", "particledt");
    writeHeader(os_size, "scalarField", "particleSize");
    writeHeader(os_velo, "vectorField", "particleVelo");
    writeHeader(os_T, "scalarFieldField", "particleTemp");
    writeHeader(os_p, "scalarFieldField", "particlePressure");
    writeHeader(os_O2, "scalarFieldField", "particleO2MassFraction");
    writeHeader(os_wetSolidVolFrac, "scalarFieldField", "wetSolidVolFraction");
    writeHeader(os_drySolidVolFrac, "scalarFieldField", "drySolidVolFraction");
    writeHeader(os_charVolFrac, "scalarFieldField", "charVolFraction");
    writeHeader(os_ashVolFrac, "scalarFieldField", "ashVolFraction");
    writeHeader(os_wetSolidMass, "scalarFieldField", "integralWetSolidMass");
    writeHeader(os_drySolidMass, "scalarFieldField", "integralDrySolidMass");
    writeHeader(os_charMass, "scalarFieldField", "integralCharMass");

    IOdictionary setParticlesDict
    (
        IOobject
        (
            "setParticlesDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );    

    IOdictionary biomassProperties
    (
        IOobject
        (
            "biomassProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    Info<< "Reading inputs" << endl;
    const vector boxMin(vector(setParticlesDict.lookup("boxMin")));
    const vector boxMax(vector(setParticlesDict.lookup("boxMax")));
    const vector ignMin(vector(setParticlesDict.lookup("ignMin")));
    const vector ignMax(vector(setParticlesDict.lookup("ignMax")));

    const label initState(readLabel(setParticlesDict.lookup("initState")));
    const scalar initTimeStep(readScalar(setParticlesDict.lookup("initTimeStep")));
    const scalar initSize(readScalar(setParticlesDict.lookup("initSize")));
    const vector initVelo(vector(setParticlesDict.lookup("initVelo")));

    const scalar meshResolution(readScalar(biomassProperties.lookup("meshResolution")));

    int numCells = round(initSize / meshResolution);
    numCells = max(5, numCells);

    Info<< "Setting initial particle fields" << endl;
    const scalarField ignTemp(numCells, readScalar(setParticlesDict.lookup("ignTemp")));
    const scalarField initTemp(numCells, readScalar(setParticlesDict.lookup("initTemp")));
    const scalarField initGaugeP(numCells,readScalar(setParticlesDict.lookup("initGaugeP")));
    const scalarField initYO2(numCells, readScalar(setParticlesDict.lookup("initYO2")));
    const scalarField initWetSolid(numCells, readScalar(setParticlesDict.lookup("initWetSolid")));
    const scalarField initDrySolid(numCells, readScalar(setParticlesDict.lookup("initDrySolid")));
    const scalarField initChar(numCells, readScalar(setParticlesDict.lookup("initChar")));
    const scalarField initAsh(numCells, readScalar(setParticlesDict.lookup("initAsh")));
    const scalarField Zero(numCells, 0.0);   


    Info << "box minimum  [m]          = " << boxMin << nl
         << "box miaximum [m]          = " << boxMax << nl
         << endl;

    scalar xMin=boxMin[0];
    scalar yMin=boxMin[1];
    scalar zMin=boxMin[2];

    scalar xMax=boxMax[0];
    scalar yMax=boxMax[1];
    scalar zMax=boxMax[2];

    scalar xMin_ign=ignMin[0];
    scalar yMin_ign=ignMin[1];
    scalar zMin_ign=ignMin[2];

    scalar xMax_ign=ignMax[0];
    scalar yMax_ign=ignMax[1];
    scalar zMax_ign=ignMax[2];


    int numSuperParticles = 0;
    forAll(mesh.C(), celli)
    {
        if( mesh.C()[celli][0] >= xMin && mesh.C()[celli][1] >= yMin && mesh.C()[celli][2] >= zMin
            && mesh.C()[celli][0] <= xMax && mesh.C()[celli][1] <= yMax && mesh.C()[celli][2] <= zMax
            )
        {
            numSuperParticles += 1;
        }
    }
    Info << "total number of super particles = " << numSuperParticles << endl;
    

    Info << "Writing particles positions " << endl;
    os_pos  << numSuperParticles << nl;
    os_pos  << '(' << nl;
    forAll(mesh.C(), celli)
    {
        if( mesh.C()[celli][0] >= xMin && mesh.C()[celli][1] >= yMin && mesh.C()[celli][2] >= zMin
            && mesh.C()[celli][0] <= xMax && mesh.C()[celli][1] <= yMax && mesh.C()[celli][2] <= zMax
            )
        {
            os_pos  << mesh.C()[celli] << token::SPACE << 0 << nl;
        }
    }
    os_pos  << ')' << nl;


    Info << "Setting ignition temperature" << endl;
    os_T  << numSuperParticles << nl;
    os_T  << '(' << nl;
    forAll(mesh.C(), celli)
    {
        if( mesh.C()[celli][0] >= xMin && mesh.C()[celli][1] >= yMin && mesh.C()[celli][2] >= zMin
            && mesh.C()[celli][0] <= xMax && mesh.C()[celli][1] <= yMax && mesh.C()[celli][2] <= zMax
            )
        {
            if( mesh.C()[celli][0] >= xMin_ign && mesh.C()[celli][1] >= yMin_ign && mesh.C()[celli][2] >= zMin_ign
            && mesh.C()[celli][0] <= xMax_ign && mesh.C()[celli][1] <= yMax_ign && mesh.C()[celli][2] <= zMax_ign
            )
            {
                os_T  << ignTemp << nl;
            }
            else
            {
                os_T   << initTemp << nl;
            }
            
        }
    }
    os_T  << ')' << nl;


    Info << "Writing other particle data" << endl;

    writeData( os_state, mesh, initState, numSuperParticles, 
                xMin, yMin, zMin,
                xMax, yMax, zMax
            );

    writeData( os_dt, mesh, initTimeStep, numSuperParticles, 
                xMin, yMin, zMin,
                xMax, yMax, zMax
            );

    writeData( os_size, mesh, initSize, numSuperParticles, 
                xMin, yMin, zMin,
                xMax, yMax, zMax
            );

    writeData( os_velo, mesh, initVelo, numSuperParticles, 
                xMin, yMin, zMin,
                xMax, yMax, zMax
            );

    writeData( os_p, mesh, initGaugeP, numSuperParticles, 
                xMin, yMin, zMin,
                xMax, yMax, zMax
            );

    writeData( os_O2, mesh, initYO2, numSuperParticles, 
                xMin, yMin, zMin,
                xMax, yMax, zMax
            );

    writeData( os_wetSolidVolFrac, mesh, initWetSolid, numSuperParticles, 
                xMin, yMin, zMin,
                xMax, yMax, zMax
            );

    writeData( os_drySolidVolFrac, mesh, initDrySolid, numSuperParticles, 
                xMin, yMin, zMin,
                xMax, yMax, zMax
            );

    writeData( os_charVolFrac, mesh, initChar, numSuperParticles, 
                xMin, yMin, zMin,
                xMax, yMax, zMax
            );

    writeData( os_ashVolFrac, mesh, initAsh, numSuperParticles, 
                xMin, yMin, zMin,
                xMax, yMax, zMax
            );

    writeData( os_wetSolidMass, mesh, Zero, numSuperParticles, 
                xMin, yMin, zMin,
                xMax, yMax, zMax
            );

    writeData( os_drySolidMass, mesh, Zero, numSuperParticles, 
                xMin, yMin, zMin,
                xMax, yMax, zMax
            );

    writeData( os_charMass, mesh, Zero, numSuperParticles, 
                xMin, yMin, zMin,
                xMax, yMax, zMax
            );

    Info << "Done setting particles" << endl;

    return(0);
}


// ************************************************************************* //