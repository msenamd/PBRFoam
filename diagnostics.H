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
	PyrolysisFront[cellI] = 0.0;
	SmolderingFront[cellI] = 0.0;

	// Flame front
	if (ft[cellI] >= ftSt)
	{ 
	    flameHeight[cellI] = mesh.C()[cellI].z();
	    flameFront[cellI] = mesh.C()[cellI].x();
	}

	// Pyrolysis front
	if (particles.pyrolysisRate[cellI] > 1e-12) 
	{ 
	    PyrolysisFront[cellI] = mesh.C()[cellI].x();
	}

	// Smoldering front
	if (particles.charOxidRate[cellI] > 1e-12) 
	{ 
	    SmolderingFront[cellI] = mesh.C()[cellI].x();
	}


}

scalar maxFlameHeight = max(flameHeight).value();
scalar maxFlameFront = max(flameFront).value();
scalar maxPyrolysisFront = max(PyrolysisFront).value();
scalar maxSmolderingFront = max(SmolderingFront).value();

outFireFront << runTime.value() << " " << maxFlameHeight << " " << maxFlameFront 
<< " " << maxPyrolysisFront << " " << maxSmolderingFront << endl;


// Other diagnostics
forAll (mesh.C(),cellI)
{
	// heat flux ahead of the flame
	QconvPUVbed[cellI] = 0.0;
	QradPUVbed[cellI] = 0.0;

	if (mesh.C()[cellI].x() > maxPyrolysisFront) 
	{ 
	    	QconvPUVbed[cellI] = particles.surfaceHeatFluxConv[cellI] / mesh.V()[cellI];
	    	QradPUVbed[cellI]  = particles.surfaceHeatFluxRad[cellI] / mesh.V()[cellI];				
	}

	//mass rates
	mFuelPUVbed[cellI] = particles.gasFuelReleaseRate[cellI] / mesh.V()[cellI];
	RR1PUVbed[cellI] = particles.dryingRate[cellI] / mesh.V()[cellI];
	RR2PUVbed[cellI] = particles.pyrolysisRate[cellI] / mesh.V()[cellI];
	RR3PUVbed[cellI] = particles.oxidPyrolysisRate[cellI] / mesh.V()[cellI];
	RR4PUVbed[cellI] = particles.charOxidRate[cellI] / mesh.V()[cellI];

	//fluid momentum
	fluidMomentum[cellI] = rho[cellI]/(1.0 - particles.packingRatio[cellI]) * U[cellI] * mesh.V()[cellI];
}
