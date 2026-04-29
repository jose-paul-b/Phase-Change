{
  rhoCp = rho * Cp;
  kappaEff = kappa + rho * Cp * turbulence->nut() / Prt;
  const surfaceScalarField rhoCpPhi = fvc::interpolate(rho * Cp) * phi;
  
  Tprev = T;
  

  Tprev = T;
  Info << "Creating TEqn" << endl;
  fvScalarMatrix TEqn(fvm::ddt(rhoCp, T) + fvm::div(rhoCpPhi, T) -
                          fvm::laplacian(kappaEff, T) ==
                      -fvc::ddt(rhoHf, alpha2));
  Info << "Solving TEqn" << endl;

  TEqn.solve();
  T.correctBoundaryConditions();
  Info << "min/max(T)=" << min(T).value() << "," << max(T).value() << ","
       << min(alpha2).value() << "," << max(alpha2).value() << endl;

  alpha2 = alpha2 + (Cp * (T - Tmelt) / Hf);

  
  alpha2 = min(max(alpha2, scalar(0.0)), scalar(1.0));
  // alpha2.correctBoundaryConditions();
  alpha2.correctBoundaryConditions();
  // recompute alpha1
  alpha1 = scalar(1.0) - alpha2;
  alpha1.correctBoundaryConditions();
  // alpha1.correctBoundaryConditions();
  if (runTime.writeTime()) {
    alpha1.write();
    alpha2.write();
  }
}
