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


// Gas Radiation absorption/emission
const volScalarField::Internal radRu = radiation->Ru();
volScalarField::Internal& radFieldRu = dQradAbs.ref();
radFieldRu = radRu;

const volScalarField::Internal T4 = pow(thermo.T(),4).ref();
const volScalarField::Internal radRp = radiation->Rp().ref()*T4;
volScalarField::Internal& radFieldRp = dQradEm.ref();
radFieldRp = radRp;

dQrad = dQradAbs - dQradEm; 

// Convective heat flux
convectiveHeatFlux = 
         - fvc::interpolate(turbulence->alphaEff())*fvc::interpolate(thermo.Cp())*fvc::snGrad(T);


// Fire front diagnostics
forAll (mesh.C(),cellI)
{
	flameHeight[cellI] = 0.0;
	flameBack[cellI] = 1e12;
	flameFront[cellI] = 0.0;
	pyrolysisBack[cellI] = 1e12;
	pyrolysisFront[cellI] = 0.0;
	smolderingBack[cellI] = 1e12;
	smolderingFront[cellI] = 0.0;

	// flame height
	if (dQ[cellI] >= 1000 && mag(g).value())
	{ 
	    flameHeight[cellI] = mag(mesh.C()[cellI] & (g/mag(g)).value() );
	}
	// flame front
	if (dQ[cellI] >= 1000)
	{ 
	    flameBack[cellI] = mesh.C()[cellI].x();
	    flameFront[cellI] = mesh.C()[cellI].x();	    
	}	
	// pyrolysis front
	if (particles.state[cellI] == 3) 
	{ 
	    pyrolysisBack[cellI] = mesh.C()[cellI].x();		
	    pyrolysisFront[cellI] = mesh.C()[cellI].x();
	}
	// smoldering front
	if (particles.state[cellI] == 5) 
	{ 
	    smolderingBack[cellI] = mesh.C()[cellI].x();		
	    smolderingFront[cellI] = mesh.C()[cellI].x();
	}
}

scalar maxFlameHeight = max(flameHeight.internalField()).value();

scalar maxFlameBack = min(flameBack.internalField()).value();
scalar maxPyrolysisBack = min(pyrolysisBack.internalField()).value();
scalar maxSmolderingBack = min(smolderingBack.internalField()).value();

scalar maxFlameFront = max(flameFront.internalField()).value();
scalar maxPyrolysisFront = max(pyrolysisFront.internalField()).value();
scalar maxSmolderingFront = max(smolderingFront.internalField()).value();

outFireFront << runTime.value() << "," << maxFlameHeight << "," << maxFlameBack << "," << maxFlameFront
<< "," << maxPyrolysisBack << "," << maxPyrolysisFront << "," << maxSmolderingBack << "," << maxSmolderingFront << endl;



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

