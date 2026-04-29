
Info << "ReadingphaseChangeProperties" << endl;  // Add
IOdictionary phaseChangeProperties(IOobject("phaseChangeProperties",
                                            runTime.constant(), mesh,
                                            IOobject::MUST_READ,
                                            IOobject::NO_WRITE));
dimensionedScalar Prt  // prantdl
    ("Prt", dimless, phaseChangeProperties);
dimensionedScalar Hf  // latentheat
    ("Hf", dimensionSet(0, 2, -2, 0, 0, 0, 0), phaseChangeProperties);
dimensionedScalar A_mushy  // mushy zone constant
    ("A_mushy", dimensionSet(1, -3, -1, 0, 0, 0, 0), phaseChangeProperties);
dimensionedScalar q  // constant inDarcyterm
    ("q", dimless, phaseChangeProperties.lookupOrDefault<scalar>("q", 0.001));
dimensionedScalar Tmelt  // melting T
    ("Tmelt", dimensionSet(0, 0, 0, 1, 0, 0, 0), phaseChangeProperties);

dimensionedScalar Tliq  // liquidus T
    ("Tliq", dimensionSet(0, 0, 0, 1, 0, 0, 0), phaseChangeProperties);
dimensionedScalar Tsol  // solidus T
    ("Tsol", dimensionSet(0, 0, 0, 1, 0, 0, 0), phaseChangeProperties);
dimensionedScalar rhoS  // solid rho
    ("rhoS", dimDensity, phaseChangeProperties);

dimensionedVector g  // solid rho
    ("g", dimLength / dimTime / dimTime, phaseChangeProperties);

dimensionedScalar Cp1  // solid rho
    ("Cp1", dimEnergy / dimTemperature / dimMass, phaseChangeProperties);
dimensionedScalar Cp2  // solid rho
    ("Cp2", dimEnergy / dimTemperature / dimMass, phaseChangeProperties);

dimensionedScalar CpS  // solid Cp
    ("CpS", dimEnergy / dimTemperature / dimMass, phaseChangeProperties);
dimensionedScalar beta  // solid rho
    ("beta", dimensionSet(0, 0, 0, -1, 0, 0, 0), phaseChangeProperties);

dimensionedScalar kappaS  // solid kappa
    ("kappaS", dimEnergy / dimTime / dimLength / dimTemperature,
     phaseChangeProperties);
Info << "Reading field p\n" << endl;
volScalarField p(IOobject("p", runTime.timeName(), mesh, IOobject::MUST_READ,
                          IOobject::AUTO_WRITE),
                 mesh);

Info << "Reading field Us\n" << endl;
volVectorField Us(IOobject("Us", runTime.timeName(), mesh, IOobject::MUST_READ,
                           IOobject::AUTO_WRITE),
                  mesh);

Info << "ReadingfieldT\n" << endl;  // Add
volScalarField T(IOobject("T", runTime.timeName(), mesh, IOobject::MUST_READ,
                          IOobject::AUTO_WRITE),
                 mesh);

volScalarField Tprev(IOobject("Tprev", runTime.timeName(), mesh,
                              IOobject::NO_READ, IOobject::NO_WRITE),
                     T);
/*
volScalarField gamma
(
IOobject
(

"gamma",
mesh.time().timeName(),
mesh,
IOobject::NO_READ,
IOobject::AUTO_WRITE
),
mesh,
dimensionedScalar("gamma",dimless,1)
);
*/

// Updategammausingtheerf function
// gamma=0.5 * Foam::erf(4.0* (T-Tmelt)/(Tliq-Tsol))+ 0.5;

Info << "Reading field U\n" << endl;
volVectorField U(IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ,
                          IOobject::AUTO_WRITE),
                 mesh);

Info << "Reading field Alpha1\n" << endl;
volScalarField alpha1(IOobject("alpha1", runTime.timeName(), mesh,
                               IOobject::MUST_READ, IOobject::AUTO_WRITE),
                      mesh);
alpha1.oldTime();
Info << "Reading field Alpha2\n" << endl;
volScalarField alpha2(IOobject("alpha2", runTime.timeName(), mesh,
                               IOobject::MUST_READ, IOobject::AUTO_WRITE),
                      mesh);
alpha2.oldTime();

Info << "Reading field Rho1\n" << endl;
volScalarField rho1(IOobject("rho1", runTime.timeName(), mesh,
                             IOobject::MUST_READ, IOobject::AUTO_WRITE),
                    mesh);

Info << "Reading field Rho2\n" << endl;
volScalarField rho2(IOobject("rho2", runTime.timeName(), mesh,
                             IOobject::MUST_READ, IOobject::AUTO_WRITE),
                    mesh);

Info << "Reading field kappa1\n" << endl;
volScalarField kappa1(IOobject("kappa1", runTime.timeName(), mesh,
                               IOobject::MUST_READ, IOobject::AUTO_WRITE),
                      mesh);

Info << "Reading field kappa2\n" << endl;
volScalarField kappa2(IOobject("kappa2", runTime.timeName(), mesh,
                               IOobject::MUST_READ, IOobject::AUTO_WRITE),
                      mesh);

#include "createPhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence(
    incompressible::turbulenceModel::New(U, phi, laminarTransport));

#include "createFvOptions.H"
#include "createMRF.H"

volScalarField Cp  // Add
    (IOobject("Cp", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT),
     alpha1*(Cp1) + alpha2 * Cp2);
Cp.oldTime();

// Needtostore rho forddt(rho,U)
volScalarField rho(IOobject("rho", runTime.timeName(), mesh,
                            IOobject::READ_IF_PRESENT),
                   alpha1*(rho1) + alpha2 * rho2);
rho.oldTime();
// Massflux
surfaceScalarField rhoPhi(IOobject("rhoPhi", runTime.timeName(), mesh,
                                   IOobject::NO_READ, IOobject::NO_WRITE),
                          fvc::interpolate(rho) * phi);

volScalarField rhoCp(IOobject("rhoCp", runTime.timeName(), mesh,
                              IOobject::NO_READ, IOobject::NO_WRITE),
                     rho* Cp);
rhoCp.oldTime();

volScalarField rhoHf(IOobject("rhoHf", runTime.timeName(), mesh,
                              IOobject::NO_READ, IOobject::NO_WRITE),
                     rho1* Hf);
rhoHf.oldTime();

volScalarField kappa  // Add
    (IOobject("kappa", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT),
     alpha1*(kappa1) + alpha2 * kappa2);
kappa.oldTime();
volScalarField kappaEff  // Add
    (IOobject("kappaEff", runTime.timeName(), mesh, IOobject::NO_READ,
              IOobject::NO_WRITE),
     kappa);
