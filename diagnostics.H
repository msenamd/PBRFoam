#include "constants.H"

if (thermo.composition().contains("O2")) 
{
    ft = (fu*s-O2+YO2Inf)/(s+YO2Inf);
}
else
{
    ft = fu;
}

ft.max(0.0);
ft.min(1.0);


flameHeight = dimensionedScalar("zero",dimLength,0.0);
PyrolysisFront = dimensionedScalar("zero",dimLength,0.0);


// Pyrolysis front is where MLR_fuel > 1e-5
// Flame front and height is where dQ >= threshold
forAll (mesh.C(),cellI)
{
	if (dQ[cellI] >= 1.0)
	{ 
	    flameHeight[cellI] = mesh.C()[cellI].z();
	    flameFront[cellI] = mesh.C()[cellI].x();
	}

	if (particles.rhoYTrans()[fuelIndex][cellI] > ROOTVSMALL) 
	{ 
	    PyrolysisFront[cellI] = mesh.C()[cellI].x();
	}
}

scalar maxFlameHeight = max(flameHeight).value();
scalar maxFlameFront = max(flameFront).value();
scalar maxPyrolysisFront = max(PyrolysisFront).value();

outFlameHeightFront <<  runTime.value() << " " << maxFlameHeight << " " << maxFlameFront << endl;
outPyrolysisFront <<  runTime.value() << " " << maxPyrolysisFront << endl;

/*
//limit heat transfer diagnostic to cells ahead of the flame only
forAll (mesh.C(),cellI)
{

	fluidMomentum[cellI] = rho[cellI] * U[cellI] * mesh.V()[cellI];

	MLR[cellI] = particles.rhoTrans()[cellI] / runTime.deltaTValue();

	QconvToFreshParticles[cellI] = 0.0;
	QradToFreshParticles[cellI] = 0.0;

	if (mesh.C()[cellI].x() > maxPyrolysisFront) 
	{ 
	    	QconvToFreshParticles[cellI] = - particles.QconvTrans()[cellI] / runTime.deltaTValue();

			QradToFreshParticles[cellI]  = particles.packingRatio()[cellI] * particles.surfaceToVolumeRatio()[cellI] / 4.0 * 
											(
												G_tot[cellI] - 4*5.670367e-8*Foam::pow(particles.surfaceTemp()[cellI],4.0)
											); 					
	}
}
*/

convectiveHeatFlux_L = 
         - fvc::interpolate(turbulence->alpha())*fvc::interpolate(thermo.Cp())*fvc::snGrad(T);

convectiveHeatFlux_T = 
         - fvc::interpolate(turbulence->alphaEff())*fvc::interpolate(thermo.Cp())*fvc::snGrad(T);