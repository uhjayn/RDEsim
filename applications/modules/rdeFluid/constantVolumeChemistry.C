/*---------------------------------------------------------------------------*\
  RDEsim constant-volume chemistry
  SPDX-License-Identifier: GPL-3.0-or-later

  State y = (Y_0 .. Y_{N-1}, T). Density rho is fixed, so
    c_i = rho Y_i / W_i,   p = rho RR T sum_i(Y_i / W_i)   (perfect gas),
    dY_i/dt = W_i omega_i / rho,
    dT/dt   = - sum_i (dY_i/dt) ea_i(T) / Cv_mix,
  where omega_i is Reaction::dNdtByV and ea_i is the absolute internal energy of
  specie i per unit mass (sensible plus formation). The right-hand side
  conserves the absolute internal energy sum_i Y_i ea_i(T) of the cell.
  The Jacobian is analytic, built from Reaction::ddNdtByVdcTp, which returns the
  derivatives with respect to concentration at fixed T and with respect to T at
  fixed concentration. At fixed density the concentrations do not depend on T,
  so no thermal-expansion correction (as in the constant-pressure model) is
  needed.
\*---------------------------------------------------------------------------*/

#include "constantVolumeChemistry.H"
#include "thermodynamicConstants.H"
#include "fluidMulticomponentThermo.H"

using namespace Foam::constant::thermodynamic;

Foam::constantVolumeChemistry::constantVolumeChemistry
(
    const multicomponentThermo& thermo,
    const dictionary& dict
)
:
    mixture_
    (
        dynamicCast<const multicomponentMixture<thermoType>>(thermo)
    ),
    thermos_(mixture_.specieThermos()),
    nSpecie_(thermos_.size()),
    reactions_(thermo.species(), thermos_, dict),
    solver_(),
    rho_(0),
    c_(nSpecie_, 0),
    work0_(nSpecie_ + 2, 0),
    work1_(nSpecie_ + 2, 0),
    Cv_(nSpecie_, 0),
    ea_(nSpecie_, 0),
    omega_(nSpecie_ + 1, 0),
    dwdc_(nSpecie_ + 1, 0),
    noReduction_()
{
    solver_ = ODESolver::New(*this, dict.subDict("ode"));

    Info<< "rdeFluid detailed chemistry: " << nSpecie_ << " species, "
        << reactions_.size() << " reactions, constant-volume ODE solver "
        << dict.subDict("ode").lookup<word>("solver") << endl;
}


Foam::constantVolumeChemistry::~constantVolumeChemistry()
{}


void Foam::constantVolumeChemistry::state
(
    const scalarField& y,
    scalar& T,
    scalar& p
) const
{
    T = y[nSpecie_];

    scalar sumYByW = 0;
    for (label i=0; i<nSpecie_; i++)
    {
        const scalar Yi = max(y[i], 0);
        sumYByW += Yi/thermos_[i].W();
        c_[i] = rho_*Yi/thermos_[i].W();
    }

    p = rho_*RR*T*sumYByW;
}


void Foam::constantVolumeChemistry::derivatives
(
    const scalar,
    const scalarField& y,
    const label li,
    scalarField& dydt
) const
{
    scalar T, p;
    state(y, T, p);

    omega_ = Zero;
    forAll(reactions_, ri)
    {
        reactions_[ri].dNdtByV(p, T, c_, li, omega_, false, noReduction_, 0);
    }

    scalar Cv = 0;
    scalar dTdt = 0;
    for (label i=0; i<nSpecie_; i++)
    {
        const scalar Yi = max(y[i], 0);
        Cv += Yi*thermos_[i].Cv(p, T);

        const scalar dYidt = omega_[i]*thermos_[i].W()/rho_;
        dydt[i] = dYidt;
        dTdt -= dYidt*thermos_[i].ea(p, T);
    }

    dydt[nSpecie_] = dTdt/Cv;
}


void Foam::constantVolumeChemistry::jacobian
(
    const scalar,
    const scalarField& y,
    const label li,
    scalarField& dfdt,
    scalarSquareMatrix& dfdy
) const
{
    scalar T, p;
    state(y, T, p);

    const label nT = nSpecie_;

    omega_ = Zero;
    for (label i=0; i<=nT; i++)
    {
        for (label j=0; j<=nT; j++)
        {
            dwdc_(i, j) = 0;
        }
    }

    forAll(reactions_, ri)
    {
        reactions_[ri].ddNdtByVdcTp
        (
            p, T, c_, li, omega_, dwdc_, false, noReduction_, 0, nT,
            work0_, work1_
        );
    }

    // Mixture Cv, dCv/dT, and the specie Cv and ea at (p, T)
    scalar Cv = 0, dCvdT = 0;
    for (label i=0; i<nSpecie_; i++)
    {
        const scalar Yi = max(y[i], 0);
        Cv_[i] = thermos_[i].Cv(p, T);
        ea_[i] = thermos_[i].ea(p, T);
        Cv += Yi*Cv_[i];
        // Cv = Cp - R for a perfect gas, so dCv/dT = dCp/dT
        dCvdT += Yi*thermos_[i].dCpdT(p, T);
    }

    // f_i = (W_i/rho) omega_i with c_j = rho Y_j/W_j, so
    // df_i/dY_j = (W_i/W_j) d omega_i/dc_j and df_i/dT = (W_i/rho) d omega_i/dT
    scalar dTdt = 0;
    for (label i=0; i<nSpecie_; i++)
    {
        const scalar Wi = thermos_[i].W();
        const scalar fi = omega_[i]*Wi/rho_;
        dTdt -= fi*ea_[i];

        for (label j=0; j<nSpecie_; j++)
        {
            dfdy(i, j) = (Wi/thermos_[j].W())*dwdc_(i, j);
        }
        dfdy(i, nT) = (Wi/rho_)*dwdc_(i, nT);
    }
    dTdt /= Cv;

    // d(dT/dt)/dY_j and d(dT/dt)/dT from dT/dt = -sum_k f_k ea_k / Cv
    for (label j=0; j<nSpecie_; j++)
    {
        scalar sum = 0;
        for (label k=0; k<nSpecie_; k++)
        {
            sum += dfdy(k, j)*ea_[k];
        }
        dfdy(nT, j) = -(sum + dTdt*Cv_[j])/Cv;
    }

    scalar sumT = 0;
    for (label k=0; k<nSpecie_; k++)
    {
        const scalar fk = omega_[k]*thermos_[k].W()/rho_;
        sumT += dfdy(k, nT)*ea_[k] + fk*Cv_[k];
    }
    dfdy(nT, nT) = -(sumT + dTdt*dCvdT)/Cv;

    // The system does not depend on time explicitly
    dfdt = Zero;
}


void Foam::constantVolumeChemistry::solve
(
    const scalar deltaT,
    const scalar rho,
    scalarField& y,
    scalar& dtTry
) const
{
    rho_ = rho;
    solver_->solve(0, deltaT, y, 0, dtTry);
}


// ************************************************************************* //
