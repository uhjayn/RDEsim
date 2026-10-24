/*---------------------------------------------------------------------------*\
  RDEsim rdeFluid energy equation
  SPDX-License-Identifier: GPL-3.0-or-later

  Kurganov-Tadmor face flux for sensible internal energy, plus the heat-release
  source rho * omega * q. Face states are the ones shockFluid cached. The
  reconstruction uses the public fvc::interpolate and the user's
  reconstruct(e) scheme. q is constant/rdeProperties chemistry/q.
\*---------------------------------------------------------------------------*/

#include "rdeFluid.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcInterpolate.H"
#include "fvmDdt.H"

void Foam::solvers::rdeFluid::thermophysicalPredictor()
{
    if (!pos.valid() || !aphiv_pos.valid() || !p_pos.valid() || !U_pos.valid())
    {
        FatalErrorInFunction
            << "shockFluid face states are not available for the energy flux."
            << exit(FatalError);
    }

    volScalarField& e = thermo_.he();
    const volScalarField rate(reactionRate());
    const dimensionedScalar zeroLambda("zero", dimless, 0);
    const volScalarField omega(rate*max(lambdaMax_ - lambda_, zeroLambda));

    const surfaceScalarField ePos
    (
        fvc::interpolate(e, pos(), "reconstruct(e)")
    );
    const surfaceScalarField eNeg
    (
        fvc::interpolate(e, neg(), "reconstruct(e)")
    );

    surfaceScalarField phiEp
    (
        "phiEp",
        aphiv_pos()*(rho_pos()*(ePos + 0.5*magSqr(U_pos())) + p_pos())
      + aphiv_neg()*(rho_neg()*(eNeg + 0.5*magSqr(U_neg())) + p_neg())
      + aSf()*(p_pos() - p_neg())
    );

    if (mesh.moving())
    {
        phiEp += mesh.phi()*(a_pos()*p_pos() + a_neg()*p_neg());
    }

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho_, e) + fvc::div(phiEp) + fvc::ddt(rho_, K)
     ==
        fvModels().source(rho_, e) + rho_*omega*q_
    );

    if (!inviscid)
    {
        const surfaceScalarField devTauDotU
        (
            "devTauDotU",
            devTau() & (a_pos()*U_pos() + a_neg()*U_neg())
        );
        EEqn += thermophysicalTransport->divq(e) + fvc::div(devTauDotU);
    }

    EEqn.relax();
    fvConstraints().constrain(EEqn);
    EEqn.solve();
    fvConstraints().constrain(e);
    thermo_.correct();
}


// ************************************************************************* //
