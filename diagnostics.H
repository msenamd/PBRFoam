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
convectiveHeatFlux_L = 
         - fvc::interpolate(turbulence->alpha())*fvc::interpolate(thermo.Cp())*fvc::snGrad(T);

convectiveHeatFlux_T = 
         - fvc::interpolate(turbulence->alphaEff())*fvc::interpolate(thermo.Cp())*fvc::snGrad(T);


flameHeight = dimensionedScalar("zero",dimLength,0.0);
PyrolysisFront = dimensionedScalar("zero",dimLength,0.0);


forAll (mesh.C(),cellI)
{
	// Flame front is where HRR > threshold
	if (dQ[cellI] >= 1e-9)
	{ 
	    flameHeight[cellI] = mesh.C()[cellI].z();
	    flameFront[cellI] = mesh.C()[cellI].x();
	}

	// Pyrolysis front is where mass production rate of fuel > threshold
	if (particles.gasFuelReleaseRate[cellI] > 1e-9) 
	{ 
	    PyrolysisFront[cellI] = mesh.C()[cellI].x();
	}

	// Smoldering front is where mass production rate of CO2 > threshold
	if (particles.CO2ReleaseRate[cellI] > 1e-9) 
	{ 
	    PyrolysisFront[cellI] = mesh.C()[cellI].x();
	}

	// Momentum of the gas phase
	fluidMomentum[cellI] = rho[cellI]/(1.0 - particles.packingRatio[cellI]) * U[cellI] * mesh.V()[cellI];
}

scalar maxFlameHeight = max(flameHeight).value();
scalar maxFlameFront = max(flameFront).value();
scalar maxPyrolysisFront = max(PyrolysisFront).value();

outFlameHeightFront <<  runTime.value() << " " << maxFlameHeight << " " << maxFlameFront << endl;
outPyrolysisFront <<  runTime.value() << " " << maxPyrolysisFront << endl;


//limit heat transfer diagnostic to cells ahead of the flame only
forAll (mesh.C(),cellI)
{
	convHeatFluxFreshFuel[cellI] = 0.0;
	radHeatFluxFreshFuel[cellI] = 0.0;

	if (mesh.C()[cellI].x() > maxPyrolysisFront) 
	{ 
	    	convHeatFluxFreshFuel[cellI] = particles.surfaceHeatFluxConv[cellI];
	    	radHeatFluxFreshFuel[cellI]  = particles.surfaceHeatFluxRad[cellI];				
	}
}
