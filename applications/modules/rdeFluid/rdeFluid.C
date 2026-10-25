/*---------------------------------------------------------------------------*\
  RDEsim rdeFluid
  SPDX-License-Identifier: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "rdeFluid.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(rdeFluid, 0);
    addToRunTimeSelectionTable(solver, rdeFluid, fvMesh);
}
}

Foam::solvers::rdeFluid::rdeFluid(fvMesh& mesh)
:
    shockFluid(mesh),
    properties_
    (
        IOobject
        (
            "rdeProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    lambda_
    (
        IOobject
        (
            "lambda",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    A_
    (
        "A",
        dimless/dimTime,
        properties_.subDict("chemistry").lookup<scalar>("A")
    ),
    beta_(properties_.subDict("chemistry").lookup<scalar>("beta")),
    Ta_
    (
        "Ta",
        dimTemperature,
        properties_.subDict("chemistry").lookup<scalar>("Ta")
    ),
    TRef_
    (
        "TRef",
        dimTemperature,
        properties_.subDict("chemistry").lookup<scalar>("TRef")
    ),
    TFloor_
    (
        "TFloor",
        dimTemperature,
        properties_.subDict("chemistry").lookup<scalar>("TFloor")
    ),
    Tignition_
    (
        "ignitionTemperature",
        dimTemperature,
        properties_.subDict("chemistry").lookup<scalar>("ignitionTemperature")
    ),
    lambdaMin_
    (
        "lambdaMin",
        dimless,
        properties_.subDict("chemistry").lookup<scalar>("lambdaMin")
    ),
    lambdaMax_
    (
        "lambdaMax",
        dimless,
        properties_.subDict("chemistry").lookup<scalar>("lambdaMax")
    ),
    q_
    (
        "q",
        dimEnergy/dimMass,
        properties_.subDict("chemistry").lookup<scalar>("q")
    ),
    sourceIntegration_
    (
        word(properties_.subDict("chemistry").lookup("sourceIntegration"))
    ),
    flowModel_(word(properties_.lookup("flowModel"))),
    equationOfState_(word(properties_.lookup("equationOfState"))),
    ddtScheme_(word(properties_.lookup("ddtScheme")))
{
    if (TRef_.value() <= 0 || TFloor_.value() <= 0)
    {
        FatalIOErrorInFunction(properties_)
            << "chemistry TRef and TFloor must be positive. "
            << "They are read from constant/rdeProperties and have no default."
            << exit(FatalIOError);
    }

    mesh.schemes().setFluxRequired(lambda_.name());
    checkHyperbolicity();
}


Foam::solvers::rdeFluid::~rdeFluid()
{}


void Foam::solvers::rdeFluid::prePredictor()
{
    shockFluid::prePredictor();
    progressEquation();
}


// ************************************************************************* //
