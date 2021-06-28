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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

	OFstream os_pos(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"positions");

	os_pos << "FoamFile" << nl;
	os_pos << "{" << nl;
	os_pos << "		version     2.0;" << nl;
	os_pos << "		format      ascii;" << nl;
	os_pos << "		class       Cloud<solidParticle>;" << nl;
	os_pos << "		location    0;" << nl;
	os_pos << "		object      positions;" << nl;
	os_pos << "}" << nl;
	os_pos << nl;
 

    OFstream os_T(runTime.path()/"0"/"lagrangian"/"vegetationBed"/"T");

    os_T << "FoamFile" << nl;
    os_T << "{" << nl;
    os_T << "     version     2.0;" << nl;
    os_T << "     format      ascii;" << nl;
    os_T << "     class       scalarFieldField;" << nl;
    os_T << "     location    0;" << nl;
    os_T << "     object      T;" << nl;
    os_T << "}" << nl;
    os_T << nl;

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
    const scalar ignTemp(readScalar(setParticlesDict.lookup("ignTemp")));
    const scalar ambientTemp(readScalar(setParticlesDict.lookup("ambientTemp")));

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


    int numParticles = 0;
    forAll(mesh.C(), celli)
    {
    	if( mesh.C()[celli][0] >= xMin && mesh.C()[celli][1] >= yMin && mesh.C()[celli][2] >= zMin
    		&& mesh.C()[celli][0] <= xMax && mesh.C()[celli][1] <= yMax && mesh.C()[celli][2] <= zMax
    		)
    	{
    		numParticles += 1;
    	}
    }

    Info << "numParticles = " << numParticles << endl;
	Info << "Writing particles positions " << endl;

	os_pos  << numParticles << nl;
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



    os_pos  << numParticles << nl;
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

    Info << "Finished setting particles positions" << endl;

    os_T  << numParticles << nl;
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
                os_T  << '(' << ignTemp << ')' << nl;
            }
            else
            {
                os_T  << '(' << ambientTemp << ')' << nl;
            }
            
        }
    }

    os_T  << ')' << nl;


	Info << "Finished setting ignition temperature" << endl;
    return(0);
}


// ************************************************************************* //
