/*---------------------------------------------------------------------------*\
  RDEsim constant-volume chemistry
  SPDX-License-Identifier: GPL-3.0-or-later

  Stiff integration of a CHEMKIN-style mechanism at fixed density and fixed
  absolute internal energy, for use as the chemistry half of a Lie operator
  split with the Kurganov-Tadmor transport step of rdeFluid. The state of one
  cell is (Y_1 .. Y_N, T). The reaction list, the specie thermo and the rate
  laws are OpenFOAM's (ReactionList, Reaction, janafThermo); the ODE solver is
  OpenFOAM's ODESolver (seulex by default). The stock chemistryModel::Standard
  integrates at constant pressure, which does not conserve the density of a
  compressible cell, so this class replaces it. It is not a copy of that class.

  The thermo type is fixed to what rdeFluid asks the case writer for:
  janaf thermo, perfect gas, sensible internal energy, sutherland transport.
\*---------------------------------------------------------------------------*/

#ifndef constantVolumeChemistry_H
#define constantVolumeChemistry_H

#include "boolList.H"
#include "ODESystem.H"
#include "ODESolver.H"
#include "ReactionList.H"
#include "multicomponentMixture.H"
#include "multicomponentThermo.H"
#include "specie.H"
#include "perfectGas.H"
#include "janafThermo.H"
#include "sensibleInternalEnergy.H"
#include "sutherlandTransport.H"
#include "thermo.H"

namespace Foam
{

class constantVolumeChemistry
:
    public ODESystem
{
public:

    typedef
        sutherlandTransport
        <
            species::thermo<janafThermo<perfectGas<specie>>, sensibleInternalEnergy>
        >
        thermoType;


private:

    const multicomponentMixture<thermoType>& mixture_;

    const PtrList<thermoType>& thermos_;

    const label nSpecie_;

    const ReactionList<thermoType> reactions_;

    autoPtr<ODESolver> solver_;

    //- Density of the cell being integrated (kg/m^3)
    mutable scalar rho_;

    mutable scalarField c_;
    mutable scalarField work0_;
    mutable scalarField work1_;
    mutable scalarField Cv_;
    mutable scalarField ea_;
    mutable scalarField omega_;
    mutable scalarSquareMatrix dwdc_;

    const List<label> noReduction_;


    void state
    (
        const scalarField& y,
        scalar& T,
        scalar& p
    ) const;


public:

    //- Construct from the multicomponent thermo and a dictionary that holds
    //  the reaction list ("reactions") and the ODE settings ("ode")
    constantVolumeChemistry
    (
        const multicomponentThermo& thermo,
        const dictionary& dict
    );

    constantVolumeChemistry(const constantVolumeChemistry&) = delete;

    virtual ~constantVolumeChemistry();


    label nEqns() const
    {
        return nSpecie_ + 1;
    }

    label nSpecie() const
    {
        return nSpecie_;
    }

    const thermoType& specieThermo(const label i) const
    {
        return thermos_[i];
    }

    label nReaction() const
    {
        return reactions_.size();
    }

    virtual void derivatives
    (
        const scalar t,
        const scalarField& y,
        const label li,
        scalarField& dydt
    ) const;

    virtual void jacobian
    (
        const scalar t,
        const scalarField& y,
        const label li,
        scalarField& dfdt,
        scalarSquareMatrix& dfdy
    ) const;

    //- Advance y = (Y, T) over deltaT at density rho. dtTry is the first
    //  step guess and returns the solver's last step estimate.
    void solve
    (
        const scalar deltaT,
        const scalar rho,
        scalarField& y,
        scalar& dtTry
    ) const;
};

}

#endif

// ************************************************************************* //
