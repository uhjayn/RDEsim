/*---------------------------------------------------------------------------*\
  RDEsim rdeFluid progress variable
  SPDX-License-Identifier: GPL-3.0-or-later

  lambda = 0 is unburned and lambda = 1 is burned. The rate is the one-step
  Arrhenius expression from constant/rdeProperties. Nothing here is a built-in
  fuel constant.
\*---------------------------------------------------------------------------*/

#include "rdeFluid.H"
#include "fvcDiv.H"
#include "fvcInterpolate.H"
#include "fvmDdt.H"
#include "fvmSup.H"

Foam::tmp<Foam::volScalarField> Foam::solvers::rdeFluid::reactionRate() const
{
    const volScalarField& T = thermo.T();
    const volScalarField safeT(max(T, TFloor_));
    const volScalarField theta(safeT/TRef_);

    // pos0 zeros the rate below ignitionTemperature. TFloor only keeps the
    // exponential argument from hitting zero; it does not ignite the mixture.
    return volScalarField::New
    (
        "reactionRate",
        A_*pow(theta, beta_)*exp(-Ta_/safeT)*pos0(T - Tignition_)
    );
}


void Foam::solvers::rdeFluid::progressEquation()
{
    if (!pos.valid() || !aphiv_pos.valid() || !rho_pos.valid())
    {
        FatalErrorInFunction
            << "shockFluid has not cached the Kurganov-Tadmor face states. "
            << "rdeFluid::prePredictor calls shockFluid::prePredictor first."
            << exit(FatalError);
    }

    const volScalarField rate(reactionRate());
    const surfaceScalarField lambdaPos
    (
        fvc::interpolate(lambda_, pos(), "reconstruct(lambda)")
    );
    const surfaceScalarField lambdaNeg
    (
        fvc::interpolate(lambda_, neg(), "reconstruct(lambda)")
    );
    const surfaceScalarField phiLambda
    (
        "phiLambda",
        aphiv_pos()*rho_pos()*lambdaPos + aphiv_neg()*rho_neg()*lambdaNeg
    );

    const dimensionedScalar zeroLambda("zero", dimless, 0);
    tmp<fvScalarMatrix> tLambdaEqn;

    if (sourceIntegration_ == "pointImplicit")
    {
        tLambdaEqn =
        (
            fvm::ddt(rho_, lambda_) + fvc::div(phiLambda)
         ==
            rho_*rate*lambdaMax_ - fvm::Sp(rho_*rate, lambda_)
        );
    }
    else if (sourceIntegration_ == "explicit")
    {
        tLambdaEqn =
        (
            fvm::ddt(rho_, lambda_) + fvc::div(phiLambda)
         ==
            rho_*rate*max(lambdaMax_ - lambda_, zeroLambda)
        );
    }
    else
    {
        FatalErrorInFunction
            << "chemistry sourceIntegration '" << sourceIntegration_
            << "' is not implemented. Implemented: pointImplicit, explicit."
            << exit(FatalError);
    }

    fvScalarMatrix& lambdaEqn = tLambdaEqn.ref();

    lambdaEqn.relax();
    fvConstraints().constrain(lambdaEqn);
    lambdaEqn.solve();
    fvConstraints().constrain(lambda_);

    lambda_.primitiveFieldRef() = max
    (
        min(lambda_.primitiveField(), lambdaMax_.value()),
        lambdaMin_.value()
    );
    lambda_.correctBoundaryConditions();
}


// ************************************************************************* //
