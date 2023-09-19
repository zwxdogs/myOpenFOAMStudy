/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "myPhieBCFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

myPhieBCFvPatchScalarField::
myPhieBCFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    phieNbrName_("undefined-phieNbr"),
    sigmaF_(0),
    sigmaS_(0),
    solidOrFluid_(0) //0为流体，1为固体
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1;
}


myPhieBCFvPatchScalarField::
myPhieBCFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    phieNbrName_(dict.lookup("phieNbr")),
    sigmaF_(0),
    sigmaS_(0),
    solidOrFluid_(0)
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }
    
    if (dict.found("sigmaF"))
    {
        sigmaF_ = readScalar(dict.lookup("sigmaF"));
        sigmaS_ = readScalar(dict.lookup("sigmaS"));
        solidOrFluid_ = readScalar(dict.lookup("solidOrFluid"));
    }
    
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1;
    }
}


myPhieBCFvPatchScalarField::
myPhieBCFvPatchScalarField
(
    const myPhieBCFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    phieNbrName_(ptf.phieNbrName_),
    sigmaF_(ptf.sigmaF_),
    sigmaS_(ptf.sigmaS_),
    solidOrFluid_(ptf.solidOrFluid_)
{}


myPhieBCFvPatchScalarField::
myPhieBCFvPatchScalarField
(
    const myPhieBCFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    phieNbrName_(ptf.phieNbrName_),
    sigmaF_(ptf.sigmaF_),
    sigmaS_(ptf.sigmaS_),
    solidOrFluid_(ptf.solidOrFluid_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void myPhieBCFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const label samplePatchi = mpp.samplePolyPatch().index();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];

    // Calculate 
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    typedef myPhieBCFvPatchScalarField thisType;

    const fvPatchScalarField& nbrTp =
        nbrPatch.lookupPatchField<volScalarField, scalar>(phieNbrName_);

    if (!isA<thisType>(nbrTp))
    {
        FatalErrorInFunction
            << "Patch field for " << internalField().name() << " on "
            << patch().name() << " is of type " << thisType::typeName
            << endl << "The neighbouring patch field " << phieNbrName_ << " on "
            << nbrPatch.name() << " is required to be the same, but is "
            << "currently of type " << nbrTp.type() << exit(FatalError);
    }

    const thisType& nbrField = refCast<const thisType>(nbrTp);

    // Swap to obtain full local values of neighbour internal field
    tmp<scalarField> nbrIntFld(new scalarField(nbrField.size(), 0.0));
    tmp<scalarField> nbrDelta(new scalarField(nbrField.size(), 0.0));

    nbrIntFld.ref() = nbrField.patchInternalField();
    nbrDelta.ref() = nbrPatch.deltaCoeffs();
    
    mpp.distribute(nbrIntFld.ref());
    mpp.distribute(nbrDelta.ref());

    tmp<scalarField> myDelta = patch().deltaCoeffs();

    if (solidOrFluid_ == 0)
    {
        this->refValue() = nbrIntFld();
        this->refGrad() = 0;

        this->valueFraction() = (sigmaS_*nbrDelta())/(sigmaS_*nbrDelta()+sigmaF_*myDelta());
    }
    else
    {
        this->refValue() = nbrIntFld();
        this->refGrad() = 0;

        this->valueFraction() = (sigmaF_*nbrDelta())/(sigmaF_*nbrDelta()+sigmaS_*myDelta());     
    }

    mixedFvPatchScalarField::updateCoeffs();

    // Restore tag
    UPstream::msgType() = oldTag;
}


void myPhieBCFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    writeEntry(os, "phieNbr", phieNbrName_);
    writeEntry(os, "sigmaS", sigmaS_);
    writeEntry(os, "sigmaF", sigmaF_);
    writeEntry(os, "solidOrFluid", solidOrFluid_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    myPhieBCFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
