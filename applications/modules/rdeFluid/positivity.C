/*---------------------------------------------------------------------------*\
  RDEsim rdeFluid positivity guard
  SPDX-License-Identifier: GPL-3.0-or-later

  The Kurganov update can give a cell more kinetic energy than its total
  energy can support. Temperature then crosses zero and the next
  shockFluid::fluxPredictor hits a SIGFPE in sqrt(gamma R T). On the long
  annulus that happened at the injector face: pressure had fallen to about
  1 kPa while |U| was about 1100 m/s, and the heat-release source in that
  cell was zero.

  This guard runs after the energy solve. For a cell below TMin it reduces
  |U| and raises sensible energy by the same amount, so e + 1/2|U|^2 stays
  the same. Only when the total energy cannot support TMin even at rest is
  energy added. Both counts are written to the log and to
  postProcessing/positivity/activations. The switch is rdeConfig positivity.
\*---------------------------------------------------------------------------*/

#include "rdeFluid.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "Pstream.H"

void Foam::solvers::rdeFluid::boundPositivity()
{
    if (!positivityEnabled_)
    {
        return;
    }

    volScalarField& e = thermo_.he();
    scalarField& ee = e.primitiveFieldRef();
    vectorField& UU = U_.primitiveFieldRef();
    scalarField& KK = K.primitiveFieldRef();
    const scalarField& TT = thermo_.T().primitiveField();
    const scalarField& Cv = thermo_.Cv().primitiveField();

    label scaled = 0;
    label added = 0;

    forAll(TT, celli)
    {
        if (TT[celli] >= positivityTMin_)
        {
            continue;
        }

        const scalar total = ee[celli] + KK[celli];
        const scalar eTarget =
            ee[celli] + Cv[celli]*(positivityTMin_ - TT[celli]);

        if (total >= eTarget && KK[celli] > SMALL)
        {
            const scalar kTarget = total - eTarget;
            UU[celli] *= sqrt(kTarget/KK[celli]);
            KK[celli] = 0.5*magSqr(UU[celli]);
            ee[celli] = total - KK[celli];
            scaled++;
        }
        else
        {
            UU[celli] = Zero;
            KK[celli] = 0;
            ee[celli] = eTarget;
            added++;
        }
    }

    if (scaled || added)
    {
        U_.correctBoundaryConditions();
        e.correctBoundaryConditions();
        thermo_.correct();
    }

    nPositivity_ += scaled + added;
    nEnergyFloor_ += added;

    if (scaled || added)
    {
        Info<< "Positivity guard at t=" << runTime.value()
            << ": scaled |U| in " << scaled
            << " cells, added energy in " << added
            << " cells, so T stays at or above " << positivityTMin_
            << " K. Running total " << nPositivity_
            << " cell-steps, energy floor " << nEnergyFloor_ << endl;
    }
}


void Foam::solvers::rdeFluid::writePositivityCount() const
{
    if (!positivityEnabled_ || !Pstream::master())
    {
        return;
    }

    const fileName directory = runTime.path()/"postProcessing"/"positivity";
    mkDir(directory);
    OFstream out(directory/"activations.dat");
    out<< "# positivity guard counts. cell-steps sum over time steps" << nl
       << "# and over PIMPLE outer correctors." << nl
       << "enabled " << (positivityEnabled_ ? "yes" : "no") << nl
       << "TMin " << positivityTMin_ << nl
       << "activations " << nPositivity_ << nl
       << "energyFloors " << nEnergyFloor_ << nl;
}


// ************************************************************************* //
