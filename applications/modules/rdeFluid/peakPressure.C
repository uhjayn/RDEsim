/*---------------------------------------------------------------------------*\
  RDEsim peak-pressure function object
  SPDX-License-Identifier: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "peakPressure.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(peakPressure, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        peakPressure,
        dictionary
    );
}
}


Foam::functionObjects::peakPressure::peakPressure
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    pMax_
    (
        IOobject
        (
            "pMax",
            mesh_.time().name(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimPressure, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    tAtPMax_
    (
        IOobject
        (
            "tAtPMax",
            mesh_.time().name(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTime, 0),
        zeroGradientFvPatchScalarField::typeName
    )
{}


Foam::functionObjects::peakPressure::~peakPressure()
{}


Foam::wordList Foam::functionObjects::peakPressure::fields() const
{
    return wordList(1, "p");
}


bool Foam::functionObjects::peakPressure::sample()
{
    if (!mesh_.foundObject<volScalarField>("p"))
    {
        return true;
    }

    const volScalarField& p = mesh_.lookupObject<volScalarField>("p");
    const scalar time = mesh_.time().value();
    scalarField& pMax = pMax_.primitiveFieldRef();
    scalarField& tAtPMax = tAtPMax_.primitiveFieldRef();
    const scalarField& pressure = p.primitiveField();

    forAll(pressure, celli)
    {
        if (pressure[celli] >= pMax[celli])
        {
            pMax[celli] = pressure[celli];
            tAtPMax[celli] = time;
        }
    }

    pMax_.correctBoundaryConditions();
    tAtPMax_.correctBoundaryConditions();
    return true;
}


bool Foam::functionObjects::peakPressure::execute()
{
    return sample();
}


bool Foam::functionObjects::peakPressure::write()
{
    pMax_.updateTimeInstance();
    tAtPMax_.updateTimeInstance();
    pMax_.write();
    tAtPMax_.write();
    return true;
}


// ************************************************************************* //
