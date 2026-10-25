/*---------------------------------------------------------------------------*\
  RDEsim acousticWallVelocity
  SPDX-License-Identifier: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "acousticWallVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "volFields.H"

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        acousticWallVelocityFvPatchVectorField
    );
}

Foam::acousticWallVelocityFvPatchVectorField::
acousticWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, fvMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    model_(word(dict.lookup("model"))),
    pName_(word(dict.lookup("p"))),
    pRef_(dict.lookup<scalar>("pRef")),
    ZReal_(0),
    resistiveFloor_(0),
    inertia_(0),
    resistance_(0),
    stiffness_(0),
    porosity_(0),
    initialDisplacement_(0),
    initialVelocity_(0),
    xi_(),
    v_()
{
    if (model_ == "resistive")
    {
        ZReal_ = dict.lookup<scalar>("ZReal");
        resistiveFloor_ = dict.lookup<scalar>("resistiveFloor");
        if (mag(ZReal_) < resistiveFloor_)
        {
            FatalIOErrorInFunction(dict)
                << "patch " << p.name() << " |ZReal| " << mag(ZReal_)
                << " is below resistiveFloor " << resistiveFloor_
                << ". A near-zero real impedance is a pressure-release wall."
                << exit(FatalIOError);
        }
    }
    else if (model_ == "helmholtz")
    {
        inertia_ = dict.lookup<scalar>("inertia");
        resistance_ = dict.lookup<scalar>("resistance");
        stiffness_ = dict.lookup<scalar>("stiffness");
        porosity_ = dict.lookup<scalar>("porosity");
        initialDisplacement_ = dict.lookup<scalar>("initialDisplacement");
        initialVelocity_ = dict.lookup<scalar>("initialVelocity");
        if (inertia_ <= 0 || porosity_ < 0)
        {
            FatalIOErrorInFunction(dict)
                << "patch " << p.name()
                << " helmholtz inertia must be positive and porosity "
                << "must be >= 0. Both come from the case dictionary."
                << exit(FatalIOError);
        }
        xi_.setSize(p.size(), initialDisplacement_);
        v_.setSize(p.size(), initialVelocity_);
    }
    else if (model_ != "slip")
    {
        FatalIOErrorInFunction(dict)
            << "patch " << p.name() << " acoustic model '" << model_
            << "' is not slip, resistive, or helmholtz."
            << exit(FatalIOError);
    }
}


Foam::acousticWallVelocityFvPatchVectorField::
acousticWallVelocityFvPatchVectorField
(
    const acousticWallVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    model_(ptf.model_),
    pName_(ptf.pName_),
    pRef_(ptf.pRef_),
    ZReal_(ptf.ZReal_),
    resistiveFloor_(ptf.resistiveFloor_),
    inertia_(ptf.inertia_),
    resistance_(ptf.resistance_),
    stiffness_(ptf.stiffness_),
    porosity_(ptf.porosity_),
    initialDisplacement_(ptf.initialDisplacement_),
    initialVelocity_(ptf.initialVelocity_),
    xi_(mapper(ptf.xi_)),
    v_(mapper(ptf.v_))
{}


Foam::acousticWallVelocityFvPatchVectorField::
acousticWallVelocityFvPatchVectorField
(
    const acousticWallVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, fvMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    model_(ptf.model_),
    pName_(ptf.pName_),
    pRef_(ptf.pRef_),
    ZReal_(ptf.ZReal_),
    resistiveFloor_(ptf.resistiveFloor_),
    inertia_(ptf.inertia_),
    resistance_(ptf.resistance_),
    stiffness_(ptf.stiffness_),
    porosity_(ptf.porosity_),
    initialDisplacement_(ptf.initialDisplacement_),
    initialVelocity_(ptf.initialVelocity_),
    xi_(ptf.xi_),
    v_(ptf.v_)
{}


void Foam::acousticWallVelocityFvPatchVectorField::map
(
    const fvPatchVectorField& ptf,
    const fieldMapper& mapper
)
{
    fixedValueFvPatchVectorField::map(ptf, mapper);
    const acousticWallVelocityFvPatchVectorField& acoustic =
        refCast<const acousticWallVelocityFvPatchVectorField>(ptf);
    mapper(xi_, acoustic.xi_);
    mapper(v_, acoustic.v_);
}


void Foam::acousticWallVelocityFvPatchVectorField::reset
(
    const fvPatchVectorField& ptf
)
{
    fixedValueFvPatchVectorField::reset(ptf);
    const acousticWallVelocityFvPatchVectorField& acoustic =
        refCast<const acousticWallVelocityFvPatchVectorField>(ptf);
    xi_.reset(acoustic.xi_);
    v_.reset(acoustic.v_);
}


void Foam::acousticWallVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const vectorField normal(patch().nf());

    if (model_ == "slip")
    {
        const vectorField internal(patchInternalField());
        operator==(internal - normal*(normal & internal));
    }
    else if (model_ == "resistive")
    {
        if (mag(ZReal_) < resistiveFloor_)
        {
            FatalErrorInFunction
                << "patch " << patch().name() << " |ZReal| " << mag(ZReal_)
                << " is below resistiveFloor " << resistiveFloor_
                << exit(FatalError);
        }
        const fvPatchField<scalar>& pressure =
            patch().lookupPatchField<volScalarField, scalar>(pName_);
        const scalarField un((pressure - pRef_)/ZReal_);
        operator==(un*normal);
    }
    else if (model_ == "helmholtz")
    {
        if (xi_.size() != size() || v_.size() != size())
        {
            FatalErrorInFunction
                << "patch " << patch().name()
                << " Helmholtz state does not match the patch size."
                << exit(FatalError);
        }
        const fvPatchField<scalar>& pressure =
            patch().lookupPatchField<volScalarField, scalar>(pName_);
        const scalar dt = db().time().deltaTValue();
        const scalarField acceleration
        (
            (pressure - pRef_ - resistance_*v_ - stiffness_*xi_)/inertia_
        );
        v_ += acceleration*dt;
        xi_ += v_*dt;
        operator==((porosity_*v_)*normal);
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::acousticWallVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "model", model_);
    writeEntry(os, "p", pName_);
    writeEntry(os, "pRef", pRef_);
    if (model_ == "resistive")
    {
        writeEntry(os, "ZReal", ZReal_);
        writeEntry(os, "resistiveFloor", resistiveFloor_);
    }
    else if (model_ == "helmholtz")
    {
        writeEntry(os, "inertia", inertia_);
        writeEntry(os, "resistance", resistance_);
        writeEntry(os, "stiffness", stiffness_);
        writeEntry(os, "porosity", porosity_);
        writeEntry(os, "initialDisplacement", initialDisplacement_);
        writeEntry(os, "initialVelocity", initialVelocity_);
    }
    writeEntry(os, "value", *this);
}


// ************************************************************************* //
