surfaceScalarField faceUB = fvc::interpolate(U^B) & mesh.Sf();

fvScalarMatrix phieEqn
(
    fvm::laplacian(phie) == fvc::div(faceUB)
);
  
phieEqn.setReference(phieRefCell, phieRefValue);

phieEqn.relax();
phieEqn.solve();

surfaceScalarField faceJ = sigmaF*(-(fvc::snGrad(phie)*mesh.magSf()) + faceUB);
surfaceVectorField faceJV = faceJ*mesh.Cf();

J = fvc::surfaceIntegrate(faceJV) - (fvc::surfaceIntegrate(faceJ)*mesh.C());

J.correctBoundaryConditions();
fvConstraints.constrain(J);

lorentzF = J^B;