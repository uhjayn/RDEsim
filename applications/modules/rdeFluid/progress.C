/*---------------------------------------------------------------------------*\
  RDEsim rdeFluid progress variable
  SPDX-License-Identifier: GPL-3.0-or-later

  lambda = 0 is unburned and lambda = 1 is burned. The rate is the one-step
  Arrhenius expression from constant/rdeProperties. Nothing here is a built-in
  fuel constant.
\*---------------------------------------------------------------------------*/

#include "rdeFluid.H"
#include "fvcDiv.H"
#include "surfaceInterpolate.H"
#include "fvmDdt.H"
#include "fvmSup.H"

Foam::tmp<Foam::volScalarField> Foam::solvers::rdeFluid::reactionRate() const
{
    const volScalarField& T = thermo.T();
    const volScalarField safeT(max(T, TFloor_));

    if (reactionModel_ == "detailed")
    {
        // Heat release comes from the species and the energy bookkeeping in
        // chemistryStep, so the progress-variable source is zero.
        return volScalarField::New
        (
            "reactionRate",
            mesh,
            dimensionedScalar(dimless/dimTime, 0)
        );
    }

    if (reactionModel_ == "twoStep")
    {
        // Release starts only where induction is complete (xi >= 1). The rate
        // is the linearised coefficient k in omega = k*(lambdaMax - lambda), so
        // the pointImplicit and explicit integrators below apply unchanged.
        // releaseOrder = 1 gives k = releaseA exp(-releaseTa/T).
        const dimensionedScalar tiny("tiny", dimless, 1e-6);
        const dimensionedScalar one("one", dimless, 1);
        const dimensionedScalar releaseA("releaseA", dimless/dimTime, releaseA_);
        const dimensionedScalar releaseTa("releaseTa", dimTemperature, releaseTa_);
        return volScalarField::New
        (
            "reactionRate",
            releaseA
           *exp(-releaseTa/safeT)
           *pow(max(lambdaMax_ - lambda_, tiny), releaseOrder_ - 1)
           *pos0(xi_() - (one - tiny))
        );
    }

    const volScalarField theta(safeT/TRef_);

    // pos0 zeros the rate below ignitionTemperature. TFloor only keeps the
    // exponential argument from hitting zero; it does not ignite the mixture.
    return volScalarField::New
    (
        "reactionRate",
        A_*pow(theta, beta_)*exp(-Ta_/safeT)*pos0(T - Tignition_)
    );
}


void Foam::solvers::rdeFluid::inductionEquation()
{
    if (!pos.valid() || !aphiv_pos.valid() || !rho_pos.valid())
    {
        FatalErrorInFunction
            << "shockFluid has not cached the Kurganov-Tadmor face states."
            << exit(FatalError);
    }

    volScalarField& xi = xi_();
    const volScalarField safeT(max(thermo.T(), TFloor_));
    const dimensionedScalar inductionA("inductionA", dimless/dimTime, inductionA_);
    const dimensionedScalar inductionTa("inductionTa", dimTemperature, inductionTa_);
    const volScalarField rate(inductionA*exp(-inductionTa/safeT));

    const surfaceScalarField xiPos(fvc::interpolate(xi, pos(), "reconstruct(xi)"));
    const surfaceScalarField xiNeg(fvc::interpolate(xi, neg(), "reconstruct(xi)"));
    const surfaceScalarField phiXi
    (
        "phiXi",
        aphiv_pos()*rho_pos()*xiPos + aphiv_neg()*rho_neg()*xiNeg
    );

    // xi grows at the induction rate and the release gate opens at xi = 1.
    // xi keeps growing up to 1 + inductionMargin and is clipped there. The
    // margin keeps numerical diffusion at the shock and noise in the limited
    // reconstruction from pushing fully induced fluid back below the gate.
    // The gate crossing time, and so the induction time, is unchanged.
    const dimensionedScalar xiTop("xiTop", dimless, 1 + inductionMargin_);
    fvScalarMatrix xiEqn
    (
        fvm::ddt(rho_, xi) + fvc::div(phiXi)
     ==
        rho_*rate*Foam::pos(xiTop - xi)
    );

    xiEqn.relax();
    fvConstraints().constrain(xiEqn);
    xiEqn.solve();
    fvConstraints().constrain(xi);

    xi.primitiveFieldRef() = max
    (
        min(xi.primitiveField(), xiTop.value()),
        scalar(0)
    );
    xi.correctBoundaryConditions();
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
