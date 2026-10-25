# Reaction models for rdeFluid: design note

Status of each model is in the table at the end. This note records what was
read in the OpenFOAM Foundation v14 sources (`/opt/openfoam14`, the
`openfoam14` package from dl.openfoam.org), what was chosen, and why.

## What exists in OpenFOAM 14

`applications/modules/shockFluid`
: Density-based, explicit Kurganov-Tadmor (KT) solver. It builds its thermo
with `psiThermo::New(mesh)` and needs a sensible internal energy `e`. Its
`prePredictor` caches the one-sided face states (`pos`, `neg`, `rho_pos`,
`aphiv_pos`, ...) that `rdeFluid` re-uses for the progress variables and for
the energy flux. `foamRun` calls `preSolve`, `adjustDeltaT`, the PIMPLE loop and
`postSolve` in that order.

`applications/modules/multicomponentFluid`
: Species-carrying fluid, but built on `isothermalFluid`: a pressure-based
PIMPLE solver with `fluidMulticomponentThermo` and a `reactionModel`
(`laminar`, which wraps `chemistryModel`). It has no KT flux. It cannot be
combined with `shockFluid` because the two are separate solver modules with
separate pressure and flux algorithms; a detonation needs the KT flux.

`src/thermophysicalModels/multicomponentThermo`
: `psiMulticomponentThermo` derives from `psiThermo`
(`psiMulticomponentThermo.H`: `virtual public psiThermo`) and is registered with
the thermo selection tables (`psiMulticomponentThermos.C`, using
`makeFluidMulticomponentThermo`). So `shockFluid`, which calls
`psiThermo::New`, can be given a multi-species `physicalProperties`
(`mixture multicomponentMixture`, `janaf`, `perfectGas`, `sensibleInternalEnergy`).
This is what makes a species-carrying `rdeFluid` possible without forking
`shockFluid`.

`src/thermophysicalModels/chemistryModel` and `src/reactionModels/laminar`
: The stock `Standard` chemistry model integrates at constant pressure with
`hs`/`T` as the thermal state. The stock laminar reaction model is only reached
from `multicomponentFluid`. For a compressible, density-based detonation the
operator-split chemistry step has to be at constant volume: after the transport
step a cell has density `rho` and internal energy `e`, and the chemistry step
must keep both.

`src/ODE`
: Stiff solvers (`seulex`, `Rosenbrock34`, `rodas23`, ...) behind `ODESolver`,
which need a derivative and a Jacobian from an `ODESystem`. These are re-used.

## Choice

Same transport for every model, different source terms.

1. **`oneStepProgress`** (default, unchanged). One progress variable `lambda`,
   Arrhenius rate. Verified bit-identical (`p`, `T`, `U`, `lambda`, `rho` after
   5e-6 s on `annulus2D`) against `main` before this change.
2. **`twoStep`**. Two passively advected scalars on the same KT flux and the
   same energy equation. `xi` is the induction progress, `d xi/dt = A_I
   exp(-Ta_I/T)`, thermally neutral. `lambda` is the release progress,
   `d lambda/dt = A_R exp(-Ta_R/T) (1-lambda)^nu`, only where `xi >= 1`, and it
   releases `q d lambda`. This is the Short-Sharpe style structure with the
   parameters left to the user. The `rdeFluid` source integrators
   (`pointImplicit`, `explicit`) are re-used for `lambda` by writing the source
   as `k (lambdaMax - lambda)`.
3. **`detailed`** (next change). `psiMulticomponentThermo` through
   `physicalProperties`, species mass fractions advected with the KT flux, and
   Lie-split constant-volume chemistry in `postSolve` with an OpenFOAM
   `ODESolver`, using the Burke et al. (2012) H2/O2 mechanism converted with
   `chemkinToFoam`.

Why the induction variable has a margin: with a limited reconstruction and a
steep shock, an advected `xi` that has just reached 1 comes back as 0.97 to 0.99
in neighbouring cells. A hard gate `xi == 1` then switches the release off in
most of the induction-complete fluid and the wave decouples (observed in a 1D
test: shock to `lambda = 0.5` distance of 4.8 mm instead of a few cells, wave
speed 0.93 of CJ). `xi` is therefore allowed to grow to `1 + inductionMargin`
(default 0.1, optional entry in `rdeProperties` chemistry) while the gate stays
at `xi >= 1`. The induction time, the time at which `xi` crosses 1, is not
changed.

## What the two-step parameters are

`fuels/hydrogen-air/fuelDict` `twoStep{}`. `inductionA` and `inductionTa` are a
least-squares fit of ln(1/tau) against 1/T to the constant-volume ignition delay
(maximum dT/dt) of stoichiometric H2-air from the Burke 2012 mechanism in
Cantera 3.2 at 20 bar, 1400 to 1900 K; the fit is within 7 % of those delays
over that range. `releaseA` matches the 10 to 90 % duration of the temperature
rise (0.27 to 0.30 microseconds, practically independent of T, hence
`releaseTa 0`). The chain-branching release of the mechanism is not first order,
so only the duration is matched. These are starting values. They are not
a validated model of H2-air detonation structure.

## Limits that come from the existing solver, not from the new models

- `rdeFluid` has one perfect gas with one `gamma` (1.4 for the example) and a
  constant `q`. For the example mixture that gives a CJ speed of 2627 m/s and a
  CJ pressure of 2.48 MPa. The equilibrium CJ state of stoichiometric H2-air at
  1 atm, 300 K for the Burke 2012 species and thermo data is 1975 m/s, 15.8 bar,
  2961 K (Cantera 3.2 equilibrium plus a minimisation of the Hugoniot speed, a
  development-time check). The reduced models are far from it because the
  burned gas has a much lower effective gamma and dissociation limits the
  temperature. Both reduced models inherit this: they
  describe a hotter, faster wave than the real mixture. The detailed model does
  not have this limit.
- The CJ state depends only on `gamma` and `q`, so the one-step and two-step
  models give the same CJ speed and differ in structure only.
- Splitting is first order (Lie). The explicit KT step uses the CFL time step,
  which resolves the transport and not necessarily the fastest chemistry; the
  ODE solvers (detailed model) or the point-implicit source (reduced models)
  handle that.

## Checks

Python: `python3 -m unittest tests.test_rdesim tests.test_runs tests.test_reaction_models`.

1D wave speed against Chapman-Jouguet: `python3 -m rdesim.cjcheck` (see the
README). It is a sanity check of heat release and flux scheme at one resolution,
not a data validation and not a cell-size check.

## Status

| Model | Selectable | Compiled | Run |
|---|---|---|---|
| `oneStepProgress` | default | yes | `annulus2D`, 1D CJ tube |
| `twoStep` | `reaction.model` or `--model` | yes | `annulus2D`, 1D CJ tube |
| `detailed` | see the next change | | |
