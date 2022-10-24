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
                OFstream& os, fvMesh& mesh, List<objectType> initValue, int totNumSuperParticles, 
                List<List<vector>> positions,
                scalar meshResolution, List<scalar> initSize                
                )
{
    os  << totNumSuperParticles << nl;
    os  << '(' << nl;

    forAll(initValue, i)
    {
        forAll(positions[i], posI)
        {
            label numCells = round(initSize[i] / meshResolution);
            numCells = max(5, numCells);

            List<objectType> init(numCells, initValue[i]);

            os  << init << nl;
        }
    }

    os  << ')' << nl;
}


template <typename objectType>
void writeData(
                OFstream& os, fvMesh& mesh, List<objectType> initValue, int totNumSuperParticles, 
                List<List<vector>> positions            
                )
{
    os  << totNumSuperParticles << nl;
    os  << '(' << nl;

    forAll(initValue, i)
    {
        forAll(positions[i], posI)
        {

            objectType init(initValue[i]);

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

    const List<label> initState(List<label>(setParticlesDict.lookup("initState")));
    const List<label> nParticlesPerSuperParticle(List<label>(setParticlesDict.lookup("nParticlesPerSuperParticle")));   
    const List<scalar> initSize(List<scalar>(setParticlesDict.lookup("initSize")));
    const List<scalar> initTimeStep(List<scalar>(setParticlesDict.lookup("initTimeStep")));
    const List<vector> initVelo(List<vector>(setParticlesDict.lookup("initVelo")));   
 
    const List<scalar> initTemp(List<scalar>(setParticlesDict.lookup("initTemp")));
    const List<scalar> initGaugeP(List<scalar>(setParticlesDict.lookup("initGaugeP")));
    const List<scalar> initYO2(List<scalar>(setParticlesDict.lookup("initYO2")));
    const List<scalar> initWetSolid(List<scalar>(setParticlesDict.lookup("initWetSolid")));
    const List<scalar> initDrySolid(List<scalar>(setParticlesDict.lookup("initDrySolid")));
    const List<scalar> initChar(List<scalar>(setParticlesDict.lookup("initChar")));
    const List<scalar> initAsh(List<scalar>(setParticlesDict.lookup("initAsh")));
    const List<scalar> Zero(initTemp.size(), 0.0); 

    const List<point> probeLocations(setParticlesDict.lookup("probeLocations"));     

    const List<scalar> separationStreamwise(List<scalar>(setParticlesDict.lookup("separationStreamwise")));
    const List<scalar> separationSpanwise(List<scalar>(setParticlesDict.lookup("separationSpanwise")));

    const List<vector> boundingBoxMin(List<vector>(setParticlesDict.lookup("boundingBoxMin")));
    const List<vector> boundingBoxMax(List<vector>(setParticlesDict.lookup("boundingBoxMax")));


// Assign ids to groups (labeled multiples of million 1000000)

    List<label> groupID(boundingBoxMin.size(), 0);

    forAll(boundingBoxMin, boxI)
    {
        groupID[boxI] = (boxI+1)*1000000;
    }

    Info << "groupID = " << groupID << endl;


// Set the number of super particles

    List<label> numSuperParticles(boundingBoxMin.size(), 0);

    List<label> nParticlesStreamwise(boundingBoxMin.size(), 0);
    List<label> nParticlesSpanwise(boundingBoxMin.size(), 0);
    List<label> nParticlesHeight(boundingBoxMin.size(), 0);

    label totNumSuperParticles(0);

    forAll(boundingBoxMin, boxI)
    {

        nParticlesStreamwise[boxI] = (boundingBoxMax[boxI].x() - boundingBoxMin[boxI].x()) / (separationStreamwise[boxI] + 2*initSize[boxI]);
        nParticlesSpanwise[boxI] = (boundingBoxMax[boxI].z() - boundingBoxMin[boxI].z()) / (separationSpanwise[boxI] + width);
        nParticlesHeight[boxI] = (boundingBoxMax[boxI].y() - boundingBoxMin[boxI].y()) / (length);

        numSuperParticles[boxI] = nParticlesStreamwise[boxI] * nParticlesSpanwise[boxI] * nParticlesHeight[boxI];

        totNumSuperParticles += numSuperParticles[boxI];
    }

    Info << "number of super particles per region = " << numSuperParticles << endl;
    Info << "total number of super particles = " << totNumSuperParticles << endl;


// Write particle positions

    List<List<vector>> positions(numSuperParticles);

    Info << "Writing particles positions " << endl;
    os_pos  << totNumSuperParticles << nl;
    os_pos  << '(' << nl;

    forAll(boundingBoxMin, boxI)
    {
        label n = 0;

        for (int i = 0 ; i < nParticlesStreamwise[boxI] ; i++ )
        {
            for (int j = 0 ; j < nParticlesHeight[boxI] ; j++ )
            {
                for (int k = 0 ; k < nParticlesSpanwise[boxI] ; k++ )
                {
                    positions[boxI][n].x() =  boundingBoxMin[boxI].x() + initSize[boxI] + i*(separationStreamwise[boxI] + 2*initSize[boxI]) ;
                    positions[boxI][n].y() =  boundingBoxMin[boxI].y() + length/2.0 + j*length ;
                    positions[boxI][n].z() =  boundingBoxMin[boxI].z() + width/2.0 + k*(separationSpanwise[boxI] + width) ;

                    os_pos  << positions[boxI][n] << token::SPACE << 0 << nl;

                    n++;        
                }
            }
        }
    }

    os_pos  << ')' << nl;


// Write particle IDs

    Info << "Writing particles IDs" << endl;

    List<List<label>> particleID(numSuperParticles);
    forAll(particleID, i)
    {
            particleID[i] = 0;
    }

    forAll(boundingBoxMin, i)
    {
        List<scalar> magDistance(numSuperParticles[i], 1e12);
        
        forAll(probeLocations, probeID)
        {
            forAll(positions[i], n)
            {
                magDistance[n] = mag(probeLocations[probeID] - positions[i][n]);
            }

            label k = findMin(magDistance, 0); 
            particleID[i][k] = groupID[i] + probeID + 1;
        }
    }

    os_ID  << totNumSuperParticles << nl;
    os_ID  << '(' << nl;
    
    forAll(boundingBoxMin, i)
    {
        forAll(positions[i], n)
        {
                os_ID  << particleID[i][n] << nl;
        }
    }

    os_ID  << ')' << nl;


// Write particle data

    Info << "Writing other particle data" << endl;


    // write global quantities

    writeData( os_state, mesh, initState, totNumSuperParticles, positions); 

    writeData( os_nParticlesPerSuperParticle, mesh, nParticlesPerSuperParticle, totNumSuperParticles, positions); 

    writeData( os_dt, mesh, initTimeStep, totNumSuperParticles, positions);

    writeData( os_size, mesh, initSize, totNumSuperParticles,  positions); 

    writeData( os_velo, mesh, initVelo, totNumSuperParticles, positions); 

    writeData( os_size0, mesh, initSize, totNumSuperParticles,  positions); 

    writeData( os_T0, mesh, initTemp, totNumSuperParticles, positions);

    writeData( os_p0, mesh, initGaugeP, totNumSuperParticles, positions);

    writeData( os_O20, mesh, initYO2, totNumSuperParticles, positions);

    writeData( os_wetSolidVolFrac0, mesh, initWetSolid, totNumSuperParticles, positions);

    writeData( os_drySolidVolFrac0, mesh, initDrySolid, totNumSuperParticles, positions);

    writeData( os_charVolFrac0, mesh, initChar, totNumSuperParticles, positions);

    writeData( os_ashVolFrac0, mesh, initAsh, totNumSuperParticles, positions);


    // write particle mesh-based quantities

    writeData( os_T, mesh, initTemp, totNumSuperParticles, positions, meshResolution, initSize);

    writeData( os_p, mesh, initGaugeP, totNumSuperParticles, positions, meshResolution, initSize);

    writeData( os_O2, mesh, initYO2, totNumSuperParticles, positions, meshResolution, initSize);

    writeData( os_wetSolidVolFrac, mesh, initWetSolid, totNumSuperParticles, positions, meshResolution, initSize);

    writeData( os_drySolidVolFrac, mesh, initDrySolid, totNumSuperParticles, positions, meshResolution, initSize);

    writeData( os_charVolFrac, mesh, initChar, totNumSuperParticles, positions, meshResolution, initSize);

    writeData( os_ashVolFrac, mesh, initAsh, totNumSuperParticles, positions, meshResolution, initSize);

    writeData( os_wetSolidMass, mesh, Zero, totNumSuperParticles, positions, meshResolution, initSize);

    writeData( os_drySolidMass, mesh, Zero, totNumSuperParticles, positions, meshResolution, initSize);

    writeData( os_charMass, mesh, Zero, totNumSuperParticles, positions, meshResolution, initSize);

    Info << "Done setting particles" << endl;

    return(0);
}


// ************************************************************************* //