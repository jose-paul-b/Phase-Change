    // Momentum predictor
	volScalarField prosityFunc= A_mushy *sqr(1.0-alpha2)/(pow3(alpha2)+q);
    MRF.correctBoundaryVelocity(U);

  //  tmp<fvVectorMatrix> tUEqn
 //   (
 //   fvm::div(phi, U)
//	  + alpha1*fvm::Sp(prosityFunc,U)//Add
  //    + MRF.DDt(U)
 //     + turbulence->divDevReff(U)
 //    ==
 //       fvOptions(U)
 //   
    
	
// create Buoyancy initialized to zero
volVectorField Buoyancy
(
    IOobject
    (
        "phiPositive",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
	dimensionedVector("zero", dimensionSet(1, -2, -2, 0, 0, 0, 0), Zero) 
);

// fill per cell: if T > Tmelt use buoyancy term, else zero
forAll(Buoyancy, cellI)
{
    
    if (T[cellI]>Tmelt.value())
    {   scalar dT=T[cellI]-Tmelt.value();
        Buoyancy[cellI] = -rho[cellI] * g.value() * beta.value() * dT;
    }
    else
    {
        Buoyancy[cellI] = vector::zero;
    }
}


	
	 tmp<fvVectorMatrix> tUEqn
	(
		fvm::ddt(rho, U)+fvm::div(rhoPhi, U)
		+alpha1*fvm::Sp(prosityFunc,U)//Add
		+MRF.DDt(rho,U)
		+turbulence->divDevRhoReff(rho,U)
		==
		fvOptions(rho, U)+Buoyancy
	);
    fvVectorMatrix& UEqn = tUEqn.ref();


    //UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);
    }
