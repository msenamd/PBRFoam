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
    OFstream os_size(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"particleSize");
    OFstream os_velo(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"particleVelo");
    OFstream os_T(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"particleTemp");
    OFstream os_p(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"particlePressure");
    OFstream os_O2(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"particleO2MassFraction");
    OFstream os_wetSolid(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"wetSolidVolFraction");
    OFstream os_drySolid(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"drySolidVolFraction");
    OFstream os_char(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"charVolFraction");
    OFstream os_ash(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"ashVolFraction");

    writeHeader(os_pos, "Cloud<solidParticle>", "positions");
    writeHeader(os_state, "labelField", "particleState");
    writeHeader(os_size, "scalarField", "particleSize");
    writeHeader(os_velo, "vectorField", "particleVelo");
    writeHeader(os_T, "scalarFieldField", "particleTemp");
    writeHeader(os_p, "scalarFieldField", "particlePressure");
    writeHeader(os_O2, "scalarFieldField", "particleO2MassFraction");
    writeHeader(os_wetSolid, "scalarFieldField", "wetSolidVolFraction");
    writeHeader(os_drySolid, "scalarFieldField", "drySolidVolFraction");
    writeHeader(os_char, "scalarFieldField", "charVolFraction");
    writeHeader(os_ash, "scalarFieldField", "ashVolFraction");


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

    Info<< "Reading the parameters from the setParticlesDict file" << endl;
    const vector boxMin(vector(setParticlesDict.lookup("boxMin")));
    const vector boxMax(vector(setParticlesDict.lookup("boxMax")));
    const vector ignMin(vector(setParticlesDict.lookup("ignMin")));
    const vector ignMax(vector(setParticlesDict.lookup("ignMax")));

    const label initState(readLabel(setParticlesDict.lookup("initState")));
    const scalar initSize(readScalar(setParticlesDict.lookup("initSize")));
    const vector initVelo(vector(setParticlesDict.lookup("initVelo")));

    const scalarField ignTemp(1, readScalar(setParticlesDict.lookup("ignTemp")));
    const scalarField initTemp(1, readScalar(setParticlesDict.lookup("initTemp")));
    const scalarField initGaugeP(1,readScalar(setParticlesDict.lookup("initGaugeP")));
    const scalarField initYO2(1, readScalar(setParticlesDict.lookup("initYO2")));
    const scalarField initWetSolid(1, readScalar(setParticlesDict.lookup("initWetSolid")));
    const scalarField initDrySolid(1, readScalar(setParticlesDict.lookup("initDrySolid")));
    const scalarField initChar(1, readScalar(setParticlesDict.lookup("initChar")));
    const scalarField initAsh(1, readScalar(setParticlesDict.lookup("initAsh")));

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

    writeData( os_wetSolid, mesh, initWetSolid, numSuperParticles, 
                xMin, yMin, zMin,
                xMax, yMax, zMax
            );

    writeData( os_drySolid, mesh, initDrySolid, numSuperParticles, 
                xMin, yMin, zMin,
                xMax, yMax, zMax
            );

    writeData( os_char, mesh, initChar, numSuperParticles, 
                xMin, yMin, zMin,
                xMax, yMax, zMax
            );

    writeData( os_ash, mesh, initAsh, numSuperParticles, 
                xMin, yMin, zMin,
                xMax, yMax, zMax
            );


    Info << "Done setting particles" << endl;

    return(0);
}


// ************************************************************************* //