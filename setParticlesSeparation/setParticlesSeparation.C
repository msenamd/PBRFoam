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
    setParticlesSeparation -- written by Mohamed Ahmed
Description
    outputs locations of lagrangian particles
    Reads in setParticlesDict
\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
#include "IOstreams.H"
#include "fvMesh.H"
#include "pointField.H"
#include "treeBoundBoxList.H"

// Auxiliary function to write header to output files
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


// Auxiliary function to write data to output files
template <typename objectType>
void writeData(
                OFstream& os, fvMesh& mesh, objectType init, int nParticles, 
                List<vector> position
                )
{
    os  << nParticles << nl;
    os  << '(' << nl;

    forAll(position, i)
    {
                os  << init << nl;
    }

    os  << ')' << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"


// Open output files

    OFstream os_pos(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"positions");
    OFstream os_ID(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"particleID");
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


// Write headers

    writeHeader(os_pos, "Cloud<solidParticle>", "positions");
    writeHeader(os_ID, "labelField", "particleID");
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


// Read Auxiliary dicts

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


// Read Input data

    Info<< "Reading inputs" << endl;

    const scalar meshResolution(readScalar(biomassProperties.lookup("meshResolution")));
    const scalar width(readScalar(biomassProperties.lookup("width")));
    const scalar length(readScalar(biomassProperties.lookup("length")));

    const label initState(readLabel(setParticlesDict.lookup("initState")));    
    const scalar initTimeStep(readScalar(setParticlesDict.lookup("initTimeStep")));
    const scalar initSize(readScalar(setParticlesDict.lookup("initSize")));
    const vector initVelo(vector(setParticlesDict.lookup("initVelo")));

    int numCells = round(initSize / meshResolution);
    numCells = max(5, numCells);

    const scalarField initTemp(numCells, readScalar(setParticlesDict.lookup("initTemp")));
    const scalarField initGaugeP(numCells,readScalar(setParticlesDict.lookup("initGaugeP")));
    const scalarField initYO2(numCells, readScalar(setParticlesDict.lookup("initYO2")));
    const scalarField initWetSolid(numCells, readScalar(setParticlesDict.lookup("initWetSolid")));
    const scalarField initDrySolid(numCells, readScalar(setParticlesDict.lookup("initDrySolid")));
    const scalarField initChar(numCells, readScalar(setParticlesDict.lookup("initChar")));
    const scalarField initAsh(numCells, readScalar(setParticlesDict.lookup("initAsh")));
    const scalarField Zero(numCells, 0.0);   

    const pointField probeLocations = setParticlesDict.lookup("probeLocations");

    const scalar separationStreamwise(readScalar(setParticlesDict.lookup("separationStreamwise")));
    const scalar separationSpanwise(readScalar(setParticlesDict.lookup("separationSpanwise")));

    const vector boundingBoxMin(vector(setParticlesDict.lookup("boundingBoxMin")));
    const vector boundingBoxMax(vector(setParticlesDict.lookup("boundingBoxMax")));


// Set the number of super particles

    label nParticlesStreamwise = (boundingBoxMax.x() - boundingBoxMin.x()) / (separationStreamwise + 2*initSize);
    label nParticlesSpanwise = (boundingBoxMax.z() - boundingBoxMin.z()) / (separationSpanwise + width);
    label nParticlesHeight = (boundingBoxMax.y() - boundingBoxMin.y()) / (length);

    Info << "nParticlesStreamwise = " << nParticlesStreamwise << endl;
    Info << "nParticlesSpanwise = " << nParticlesSpanwise << endl;
    Info << "nParticlesHeight = " << nParticlesHeight << endl;


    int nParticles = nParticlesStreamwise * nParticlesSpanwise * nParticlesHeight;

    Info << "total number of super particles = " << nParticles << endl;

    List<vector> position(nParticles);

// Write particle positions

    Info << "Writing particles positions " << endl;
    os_pos  << nParticles << nl;
    os_pos  << '(' << nl;

    int n = 0;
    for (int i = 0 ; i < nParticlesStreamwise ; i++ )
    {
        for (int j = 0 ; j < nParticlesHeight ; j++ )
        {
            for (int k = 0 ; k < nParticlesSpanwise ; k++ )
            {
                position[n].x() =  boundingBoxMin.x() + initSize + i*(separationStreamwise + 2*initSize) ;
                position[n].y() =  boundingBoxMin.y() + length/2.0 + j*length ;
                position[n].z() =  boundingBoxMin.z() + width/2.0 + k*(separationSpanwise + width) ;

                os_pos  << position[n]
                        << token::SPACE << 0 
                        << nl;

                n++;        
            }
        }
    }
    os_pos  << ')' << nl;


// Write particle IDs

    Info << "Writing particles IDs" << endl;

    List<label> particleID(nParticles, 0);
    List<scalar> magDistance(nParticles, 1e12);

    forAll(probeLocations, probeID)
    {
        forAll(position, i)
        {
            magDistance[i] = mag(probeLocations[probeID] - position[i]);
        }
        label k = findMin(magDistance, 0); 
        particleID[k] = probeID + 1;
    }


    os_ID  << nParticles << nl;
    os_ID  << '(' << nl;
    forAll(position, i)
    {
            os_ID  << particleID[i] << nl;
    }
    os_ID  << ')' << nl;


// Write particle data

    Info << "Writing other particle data" << endl;

    writeData( os_T, mesh, initTemp, nParticles, position);

    writeData( os_state, mesh, initState, nParticles, position); 

    writeData( os_dt, mesh, initTimeStep, nParticles, position);

    writeData( os_size, mesh, initSize, nParticles,  position); 

    writeData( os_velo, mesh, initVelo, nParticles, position); 

    writeData( os_p, mesh, initGaugeP, nParticles, position);

    writeData( os_O2, mesh, initYO2, nParticles, position);

    writeData( os_wetSolidVolFrac, mesh, initWetSolid, nParticles, position);

    writeData( os_drySolidVolFrac, mesh, initDrySolid, nParticles, position);

    writeData( os_charVolFrac, mesh, initChar, nParticles, position);

    writeData( os_ashVolFrac, mesh, initAsh, nParticles, position);

    writeData( os_wetSolidMass, mesh, Zero, nParticles, position);

    writeData( os_drySolidMass, mesh, Zero, nParticles, position);

    writeData( os_charMass, mesh, Zero, nParticles, position);

    Info << "Done setting particles" << endl;


    return(0);
}


// ************************************************************************* //