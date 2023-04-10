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
    Reads in setParticlesDict.
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
                OFstream& os, fvMesh& mesh, objectType initValue, int numSuperParticles, 
                List<vector> positions,
                scalar meshResolution, scalar initSize                
                )
{
    os  << numSuperParticles << nl;
    os  << '(' << nl;


        forAll(positions, posI)
        {
            label numCells = round(initSize / meshResolution);
            numCells = max(5, numCells);

            List<objectType> init(numCells, initValue);

            os  << init << nl;
        }

    os  << ')' << nl;
}


template <typename objectType>
void writeData(
                OFstream& os, fvMesh& mesh, objectType initValue, int numSuperParticles, 
                List<vector> positions            
                )
{
    os  << numSuperParticles << nl;
    os  << '(' << nl;


        forAll(positions, posI)
        {

            objectType init(initValue);

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
    OFstream os_nParticlesPerSuperParticle(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"nParticlesPerSuperParticle");    
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

    OFstream os_size0(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"particleSize0");
    OFstream os_T0(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"particleTemp0");
    OFstream os_p0(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"particlePressure0");
    OFstream os_O20(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"particleO2MassFraction0");
    OFstream os_wetSolidVolFrac0(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"wetSolidVolFraction0");
    OFstream os_drySolidVolFrac0(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"drySolidVolFraction0");
    OFstream os_charVolFrac0(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"charVolFraction0");
    OFstream os_ashVolFrac0(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"ashVolFraction0");

// Write headers

    writeHeader(os_pos, "Cloud<solidParticle>", "positions");
    writeHeader(os_ID, "labelField", "particleID");
    writeHeader(os_state, "labelField", "particleState");
    writeHeader(os_nParticlesPerSuperParticle, "labelField", "nParticlesPerSuperParticle");
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

    writeHeader(os_size0, "scalarField", "particleSize0");
    writeHeader(os_T0, "scalarField", "particleTemp0");
    writeHeader(os_p0, "scalarField", "particlePressure0");
    writeHeader(os_O20, "scalarField", "particleO2MassFraction0");
    writeHeader(os_wetSolidVolFrac0, "scalarField", "wetSolidVolFraction0");
    writeHeader(os_drySolidVolFrac0, "scalarField", "drySolidVolFraction0");
    writeHeader(os_charVolFrac0, "scalarField", "charVolFraction0");
    writeHeader(os_ashVolFrac0, "scalarField", "ashVolFraction0");

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
    const scalar width(readScalar(biomassProperties.lookup("width")));      // needs to be a list
    const scalar length(readScalar(biomassProperties.lookup("length")));    // needs to be a list

    const label initState(readLabel(setParticlesDict.lookup("initState")));
    const scalar initSize(readScalar(setParticlesDict.lookup("initSize")));
    const scalar initTimeStep(readScalar(setParticlesDict.lookup("initTimeStep")));
    const vector initVelo((setParticlesDict.lookup("initVelo")));   
 
    const scalar initTemp(readScalar(setParticlesDict.lookup("initTemp")));
    const scalar initGaugeP(readScalar(setParticlesDict.lookup("initGaugeP")));
    const scalar initYO2(readScalar(setParticlesDict.lookup("initYO2")));
    const scalar initWetSolid(readScalar(setParticlesDict.lookup("initWetSolid")));
    const scalar initDrySolid(readScalar(setParticlesDict.lookup("initDrySolid")));
    const scalar initChar(readScalar(setParticlesDict.lookup("initChar")));
    const scalar initAsh(readScalar(setParticlesDict.lookup("initAsh")));
    const scalar Zero(0.0); 

    const List<point> probeLocations(setParticlesDict.lookup("probeLocations"));     

    const scalar separationStreamwise(readScalar(setParticlesDict.lookup("separationStreamwise")));
    const scalar separationSpanwise(readScalar(setParticlesDict.lookup("separationSpanwise")));

    const vector boundingBoxMin((setParticlesDict.lookup("boundingBoxMin")));
    const vector boundingBoxMax((setParticlesDict.lookup("boundingBoxMax")));


// Assign ids to groups (labeled multiples of million 1000000)

    label groupID( 0);

    groupID = 1000000;


    Info << "groupID = " << groupID << endl;


// Find the number and positions of the original particles

    label numParticles(0);

    label nParticlesStreamwise(0);
    label nParticlesSpanwise(0);
    label nParticlesHeight(0);

    label totNumParticles(0);

    nParticlesStreamwise = (boundingBoxMax.x() - boundingBoxMin.x()) / (separationStreamwise + 2*initSize);
    nParticlesSpanwise = (boundingBoxMax.z() - boundingBoxMin.z()) / (separationSpanwise + width);
    nParticlesHeight = (boundingBoxMax.y() - boundingBoxMin.y()) / (length);

    numParticles = nParticlesStreamwise * nParticlesSpanwise * nParticlesHeight;

    totNumParticles += numParticles;


    List<vector> positions_orig(numParticles);

    Info << "Writing particles positions " << endl;

    label n = 0;

    for (int i = 0 ; i < nParticlesStreamwise ; i++ )
    {
        for (int j = 0 ; j < nParticlesHeight ; j++ )
        {
            for (int k = 0 ; k < nParticlesSpanwise ; k++ )
            {
                positions_orig[n].x() =  boundingBoxMin.x() + initSize + i*(separationStreamwise + 2*initSize) ;
                positions_orig[n].y() =  boundingBoxMin.y() + length/2.0 + j*length ;
                positions_orig[n].z() =  boundingBoxMin.z() + width/2.0 + k*(separationSpanwise + width) ;

                n++;        
            }
        }
    }


// Find cells containing each particle

    List<label> cellIDs(totNumParticles);

    forAll(positions_orig, posI)
    {
        cellIDs[posI] = mesh.findCell(positions_orig[posI]);
    }

// reduce the list if multiple particles exist within same cell

    List<label> uniqueCellIDs;
    List<label> nParticlesPerSuperParticle;

    std::sort(cellIDs.begin(), cellIDs.end());

for (int i = 0; i < cellIDs.size(); ++i)
{
    const int labelValue = cellIDs[i];

    // find the position where the current label should be inserted in the unique list
    int repetitionIndex = 0;
    while (repetitionIndex < uniqueCellIDs.size() && uniqueCellIDs[repetitionIndex] < labelValue)
    {
        ++repetitionIndex;
    }

    if (repetitionIndex == uniqueCellIDs.size() || uniqueCellIDs[repetitionIndex] != labelValue)
    {
        // if the label is not in the unique list, append it and set the count to 1 in the repetition list
        uniqueCellIDs.append(labelValue);
        nParticlesPerSuperParticle.append(1);
    }
    else
    {
        // if the label is already in the unique list, increment its count in the repetition list
        nParticlesPerSuperParticle[repetitionIndex] += 1;
    }
}

    label numSuperParticles = uniqueCellIDs.size();

    Info << "total number of super particles = " << numSuperParticles << endl;


// Write particle positions    

    List<vector> positions(numSuperParticles);

    Info << "Writing particles positions " << endl;
    os_pos  << numSuperParticles << nl;
    os_pos  << '(' << nl;
    forAll(uniqueCellIDs, n)
    {
        positions[n] = mesh.C()[uniqueCellIDs[n]];
        os_pos  << positions[n] << token::SPACE << 0 << nl;
    }
    os_pos  << ')' << nl;


// Write  nParticlesPerSuperParticle;

    os_nParticlesPerSuperParticle  << numSuperParticles << nl;
    os_nParticlesPerSuperParticle  << '(' << nl;
    forAll(nParticlesPerSuperParticle, n)
    {
            os_nParticlesPerSuperParticle  << nParticlesPerSuperParticle[n] << nl;
    }
    os_nParticlesPerSuperParticle  << ')' << nl;



// Write particle IDs

    Info << "Writing particles IDs" << endl;

    List<label> particleID(numSuperParticles);
    forAll(particleID, i)
    {
            particleID[i] = 0;
    }


    List<scalar> magDistance(numSuperParticles, 1e12);
    
    forAll(probeLocations, probeID)
    {
        forAll(positions, n)
        {
            magDistance[n] = mag(probeLocations[probeID] - positions[n]);
        }

        label k = findMin(magDistance, 0); 
        particleID[k] = groupID + probeID + 1;
    }


    os_ID  << numSuperParticles << nl;
    os_ID  << '(' << nl;
    

    forAll(positions, n)
    {
            os_ID  << particleID[n] << nl;
    }

    os_ID  << ')' << nl;


// Write particle data

    Info << "Writing other particle data" << endl;


    // write global quantities

    writeData( os_state, mesh, initState, numSuperParticles, positions); 

    writeData( os_dt, mesh, initTimeStep, numSuperParticles, positions);

    writeData( os_size, mesh, initSize, numSuperParticles,  positions); 

    writeData( os_velo, mesh, initVelo, numSuperParticles, positions); 

    writeData( os_size0, mesh, initSize, numSuperParticles,  positions); 

    writeData( os_T0, mesh, initTemp, numSuperParticles, positions);

    writeData( os_p0, mesh, initGaugeP, numSuperParticles, positions);

    writeData( os_O20, mesh, initYO2, numSuperParticles, positions);

    writeData( os_wetSolidVolFrac0, mesh, initWetSolid, numSuperParticles, positions);

    writeData( os_drySolidVolFrac0, mesh, initDrySolid, numSuperParticles, positions);

    writeData( os_charVolFrac0, mesh, initChar, numSuperParticles, positions);

    writeData( os_ashVolFrac0, mesh, initAsh, numSuperParticles, positions);


    // write particle mesh-based quantities

    writeData( os_T, mesh, initTemp, numSuperParticles, positions, meshResolution, initSize);

    writeData( os_p, mesh, initGaugeP, numSuperParticles, positions, meshResolution, initSize);

    writeData( os_O2, mesh, initYO2, numSuperParticles, positions, meshResolution, initSize);

    writeData( os_wetSolidVolFrac, mesh, initWetSolid, numSuperParticles, positions, meshResolution, initSize);

    writeData( os_drySolidVolFrac, mesh, initDrySolid, numSuperParticles, positions, meshResolution, initSize);

    writeData( os_charVolFrac, mesh, initChar, numSuperParticles, positions, meshResolution, initSize);

    writeData( os_ashVolFrac, mesh, initAsh, numSuperParticles, positions, meshResolution, initSize);

    writeData( os_wetSolidMass, mesh, Zero, numSuperParticles, positions, meshResolution, initSize);

    writeData( os_drySolidMass, mesh, Zero, numSuperParticles, positions, meshResolution, initSize);

    writeData( os_charMass, mesh, Zero, numSuperParticles, positions, meshResolution, initSize);

    Info << "Done setting particles" << endl;

    return(0);

    
}


// ************************************************************************* //