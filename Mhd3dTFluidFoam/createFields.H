#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;

autoPtr<fluidThermo> pThermo(fluidThermo::New(mesh));
fluidThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField& p = thermo.p();

IOdictionary physicalProperties
(
    IOobject
    (
        "physicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )            
);
IOdictionary electronicPropertiesF
(
    IOobject
    (
        "electronicPropertiesF",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )            
);

/* dimensionedScalar rhoCS
(
    "rhoCS",
    dimDensity,
    physicalProperties 
); */

dimensionedScalar sigmaF
(
    "sigmaF",
    dimensionSet(-1, -3, 3, 0, 0, 2, 0),
    electronicPropertiesF 
);

dimensionedVector B
(
    "B",
    dimensionSet(1, 0, -2, 0, 0, -1, 0),
    electronicPropertiesF
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "compressibleCreatePhi.H"

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::momentumTransportModel> turbulence
(
    compressible::momentumTransportModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<fluidThermophysicalTransportModel> thermophysicalTransport
(
    fluidThermophysicalTransportModel::New(turbulence(), thermo)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"
#include "readpRef.H"


Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField J
(
    IOobject
    (
        "J",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),  
    mesh       
);

volScalarField phie
(
    IOobject
    (
        "phie",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),  
    mesh        
);

volVectorField lorentzF = J^B; //定义洛伦兹力矢量

pressureReference pressureReference
(
    p,
    p_rgh,
    pimple.dict(),
    thermo.incompressible()
);

mesh.schemes().setFluxRequired(p_rgh.name());

//电势参考值
label phieRefCell = 0;
scalar phieRefValue = 0.0;
setRefCell(phie, mesh.solution().dict().subDict("PIMPLE"), phieRefCell, phieRefValue);
mesh.schemes().setFluxRequired(phie.name());

hydrostaticInitialisation
(
    p_rgh,
    p,
    rho,
    U,
    gh,
    ghf,
    pRef,
    thermo,
    pimple.dict()
);

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

dimensionedScalar initialMass = fvc::domainIntegrate(rho);

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"
#include "checkRadiationModel.H"
