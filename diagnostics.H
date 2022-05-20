#include "constants.H"

// Mixture fraction
if (thermo.composition().contains("O2")) 
{
    ft = (s*fu-O2+YO2Inf)/(s*YFInf+YO2Inf);
}
else
{
    ft = fu;
}

ft.max(0.0);
ft.min(1.0);


// Convective heat flux
convectiveHeatFlux = 
         - fvc::interpolate(turbulence->alphaEff())*fvc::interpolate(thermo.Cp())*fvc::snGrad(T);


// Fire front diagnostics
forAll (mesh.C(),cellI)
{
	flameHeight[cellI] = 0.0;
	flameFront[cellI] = 0.0;
	pyrolysisBack[cellI] = 1e12;
	pyrolysisFront[cellI] = 0.0;
	smolderingFront[cellI] = 0.0;

	// flame front
	if (dQ[cellI] >= 1000)
	{ 
	    flameHeight[cellI] = mag(mesh.C()[cellI] & (g/mag(g)).value() );
	    flameFront[cellI] = mesh.C()[cellI].x();
	}
	// pyrolysis front
	if (particles.state[cellI] == 3) 
	{ 
	    pyrolysisFront[cellI] = mesh.C()[cellI].x();
	}
	// pyrolysis back
	if (particles.state[cellI] == 3) 
	{ 
	    pyrolysisBack[cellI] = mesh.C()[cellI].x();
	}	
	// smoldering front
	if (particles.state[cellI] == 5) 
	{ 
	    smolderingFront[cellI] = mesh.C()[cellI].x();
	}
}

scalar maxFlameHeight = max(flameHeight).value();
scalar maxFlameFront = max(flameFront).value();
scalar maxpyrolysisBack = min(pyrolysisBack).value();
scalar maxPyrolysisFront = max(pyrolysisFront).value();
scalar maxSmolderingFront = max(smolderingFront).value();

outFireFront << runTime.value() << "," << maxFlameHeight << "," << maxFlameFront
<< "," << maxpyrolysisBack << "," << maxPyrolysisFront << "," << maxSmolderingFront << endl;



// Other diagnostics

label nCells = 0;

forAll (mesh.C(),cellI)
{
	// heat flux ahead of the flame (drying particles)
	QconvPUVbed[cellI] = 0.0;
	QradPUVbed[cellI] = 0.0;

	if (particles.state[cellI] == 1) 
	{ 
			nCells = nCells + 1;
	    	QconvPUVbed[cellI] = particles.surfaceHeatFluxConv[cellI] / mesh.V()[cellI];
	    	QradPUVbed[cellI]  = particles.surfaceHeatFluxRad[cellI] / mesh.V()[cellI];				
	}

	//fluid momentum
	fluidMomentum[cellI] = rho[cellI]/(1.0 - particles.packingRatio[cellI]) * U[cellI] * mesh.V()[cellI];
}

QconvPUVbed = QconvPUVbed/max(nCells,1);
QradPUVbed = QradPUVbed/max(nCells,1);


// Radiation absorption/emission by the dispersed particles
radEmPUVbed = particles.Ep();
radAbsPUVbed = particles.ap() * 4*G_projected;

