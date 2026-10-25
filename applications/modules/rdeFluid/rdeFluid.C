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

namespace
{

// The one-step constants must be in rdeProperties for oneStepProgress and have
// no default. Other models do not use them.
Foam::scalar oneStepEntry
(
    const Foam::dictionary& chemistry,
    const Foam::word& model,
    const Foam::word& key
)
{
    return model == "oneStepProgress"
      ? chemistry.lookup<Foam::scalar>(key)
      : chemistry.lookupOrDefault<Foam::scalar>(key, 1);
}

Foam::scalar twoStepEntry
(
    const Foam::dictionary& chemistry,
    const Foam::word& model,
    const Foam::word& key
)
{
    return model == "twoStep"
      ? chemistry.lookup<Foam::scalar>(key)
      : chemistry.lookupOrDefault<Foam::scalar>(key, 0);
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
    reactionModel_
    (
        properties_.lookupOrDefault<word>("reactionModel", "oneStepProgress")
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
    xi_(),
    inductionA_
    (
        twoStepEntry(properties_.subDict("chemistry"), reactionModel_, "inductionA")
    ),
    inductionTa_
    (
        twoStepEntry(properties_.subDict("chemistry"), reactionModel_, "inductionTa")
    ),
    releaseA_
    (
        twoStepEntry(properties_.subDict("chemistry"), reactionModel_, "releaseA")
    ),
    releaseTa_
    (
        twoStepEntry(properties_.subDict("chemistry"), reactionModel_, "releaseTa")
    ),
    releaseOrder_
    (
        twoStepEntry(properties_.subDict("chemistry"), reactionModel_, "releaseOrder")
    ),
    inductionMargin_
    (
        properties_.subDict("chemistry").lookupOrDefault<scalar>
        (
            "inductionMargin",
            0.1
        )
    ),
    A_
    (
        "A",
        dimless/dimTime,
        oneStepEntry(properties_.subDict("chemistry"), reactionModel_, "A")
    ),
    beta_(oneStepEntry(properties_.subDict("chemistry"), reactionModel_, "beta")),
    Ta_
    (
        "Ta",
        dimTemperature,
        oneStepEntry(properties_.subDict("chemistry"), reactionModel_, "Ta")
    ),
    TRef_
    (
        "TRef",
        dimTemperature,
        oneStepEntry(properties_.subDict("chemistry"), reactionModel_, "TRef")
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
        oneStepEntry
        (
            properties_.subDict("chemistry"),
            reactionModel_,
            "ignitionTemperature"
        )
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
    ddtScheme_(word(properties_.lookup("ddtScheme"))),
    positivityEnabled_(false),
    positivityTMin_(0),
    nPositivity_(0),
    nEnergyFloor_(0)
{
    if (TRef_.value() <= 0 || TFloor_.value() <= 0)
    {
        FatalIOErrorInFunction(properties_)
            << "chemistry TRef and TFloor must be positive. "
            << "They are read from constant/rdeProperties and have no default."
            << exit(FatalIOError);
    }

    if (properties_.found("positivity"))
    {
        const dictionary& positivity = properties_.subDict("positivity");
        positivityEnabled_ = positivity.lookup<bool>("enabled");
        if (positivityEnabled_)
        {
            positivityTMin_ = positivity.lookup<scalar>("TMin");
            if (positivityTMin_ <= 0)
            {
                FatalIOErrorInFunction(positivity)
                    << "positivity TMin must be positive. "
                    << "It is a floor used when kinetic energy would make "
                    << "temperature non-positive, not a default gas temperature."
                    << exit(FatalIOError);
            }
            Info<< "Positivity guard is on. |U| is scaled when T would fall "
                << "below " << positivityTMin_ << " K, keeping e + 1/2|U|^2. "
                << "Each activation is counted." << nl;
        }
    }

    if (reactionModel_ == "twoStep")
    {
        if
        (
            inductionA_ <= 0 || releaseA_ <= 0
         || releaseOrder_ <= 0 || inductionTa_ < 0 || releaseTa_ < 0
         || inductionMargin_ < 0
        )
        {
            FatalIOErrorInFunction(properties_)
                << "twoStep needs inductionA, releaseA, releaseOrder > 0 and "
                << "inductionTa, releaseTa >= 0."
                << exit(FatalIOError);
        }

        xi_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "xi",
                    runTime.name(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
        mesh.schemes().setFluxRequired(xi_().name());
    }
    else if (reactionModel_ != "oneStepProgress")
    {
        FatalIOErrorInFunction(properties_)
            << "reactionModel '" << reactionModel_ << "' is not implemented. "
            << "Implemented: oneStepProgress, twoStep."
            << exit(FatalIOError);
    }

    Info<< "rdeFluid reaction model: " << reactionModel_ << endl;

    mesh.schemes().setFluxRequired(lambda_.name());
    checkHyperbolicity();
}


Foam::solvers::rdeFluid::~rdeFluid()
{
    if (positivityEnabled_)
    {
        Info<< "Positivity guard total cell-steps: " << nPositivity_
            << ", energy floor: " << nEnergyFloor_ << endl;
    }
}


void Foam::solvers::rdeFluid::prePredictor()
{
    shockFluid::prePredictor();
    if (xi_.valid())
    {
        inductionEquation();
    }
    progressEquation();
}


// ************************************************************************* //
