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
3. **`detailed`**. A multicomponent thermo (`multicomponentMixture` with
   `janaf`/`perfectGas`/`sensibleInternalEnergy`) through `physicalProperties`;
   the nine species mass fractions of the Burke et al. (2012) H2/O2 mechanism,
   advected with the same KT face states as `rho`, `U` and `e`; and Lie-split
   constant-volume chemistry in `postSolve` with an OpenFOAM `ODESolver`
   (`seulex` by default). Details and what was checked are in the section
   "The detailed model" below.

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

## The detailed model

**Why not `multicomponentFluid` or the stock chemistry.** In OpenFOAM 14
`multicomponentFluid` derives from `isothermalFluid` and is pressure-based, so
it cannot share the explicit KT flux of `shockFluid`/`rdeFluid`. The stock
`chemistryModels::Standard` integrates at constant pressure, which is wrong
for a cell whose density is fixed by the explicit transport step. What is
reused from OpenFOAM: `psiMulticomponentThermo` (and the `multicomponentMixture`
it builds from the `thermo` file), `ReactionList` for the reaction file,
`Reaction::dNdtByV` and `ddNdtByVdcTp` for production rates and Jacobian, and
`ODESolver`. New: `constantVolumeChemistry` (an `ODESystem` for `[Y_i, T]` at
fixed density, `d(Y_i)/dt = W_i omega_i / rho`, `dT/dt = -sum ea_i dY_i/dt /
Cv`, `ea_i` the absolute internal energy) and the species and chemistry-step code in `detailed.C`.

**Algorithm per time step.** (1) `rdeFluid` advances `rho`, `U`, `e` and the
species `Y_i` with the KT flux. (2) In `postSolve` every cell with `T >= TMin`
is advanced over `deltaT` at fixed `rho` and fixed absolute internal energy by
the ODE solver; `Y_i` are renormalised and the transported sensible internal
energy is changed by `-sum hf_i dY_i`, so the absolute energy
`sum Y_i (es_i + hf_i)` is kept. (3) `thermo.correct()` and the pressure
update. `lambda` is kept only as a diagnostic, `1 - Y_fuel/Y_fuel,unburned`, so
the optimizer objectives, `lambdaMax` and the visualisation keep working. The
split is first order (Lie): the explicit KT step is limited by the CFL
condition and does not resolve the chemistry; the ODE solver does.

**Interaction with the positivity guard.** `boundPositivity()` runs in the
thermophysical predictor, before the chemistry step, and acts on `U` and `e`
only; it does not touch species. The chemistry step then conserves the
absolute energy, so the guard and the chemistry do not fight.

**Mechanism data.** `fuels/hydrogen-air/detailed/burke2012/`: `thermo` and
`reactions` generated by `make_mechanism.py` (a developer step, needs PyYAML
and `chemkinToFoam`) from `source/Burke2012.yaml`, the Cantera YAML conversion
of Burke, Chaos, Ju, Dryer, Klippenstein, "Comprehensive H2/O2 kinetic model
for high-pressure combustion", Int. J. Chem. Kinet. 44 (2012) 444-474, as
distributed in the Caltech Explosion Dynamics Laboratory SDToolbox. Only the
H2-air subset is kept (H, H2, O, OH, H2O, O2, HO2, H2O2, N2; 23 of the 27
reactions; the four dropped reactions exist only for argon and helium baths).
Every rate constant, third-body efficiency, falloff parameter and NASA-7
coefficient is copied unchanged. The published file asks the user to choose the bath gas for
H+O2(+M)=HO2(+M); it is shipped set to N2, which is what air needs, and that is
the block that was converted. The source files in `source/` keep the authors'
notice and the reference they ask for. One deliberate change, described in the
script: N2 is published for 300 to 5000 K and the other eight species for 200 K
upward; N2's lower limit is set to 200 K without changing its coefficients, so
that cold gas at the 300 K ambient temperature does not trigger a janaf range
warning on every call (about 110000 warning lines in a 5 microsecond run
before). Between 200 and 300 K that is an extrapolation of the published
polynomial (cp about 1.1 % under tabulated N2 at 200 K); no part of a
detonation is that cold. The upper limit stays 3500 K for most species; a cell
hotter than that makes OpenFOAM print its range warning, which is then real.

**Limits.** Euler only: the mechanism files carry no transport data, so
`flowModel navierStokes` is refused at configuration time rather than silently
staying inviscid. No `progressKernel` ignition (there is no progress
variable to preset). One perfect-gas equation of state; no real-gas or
high-pressure correction, although the mechanism itself is a high-pressure
one. No NOx. The hyperbolicity check evaluates the unburned composition at the
inlet and kernel temperatures and checks the janaf range for them. The
burnt-gas composition is not known before the run, so the solver checks the
initial fields only and does not monitor gamma or c^2 afterwards.

## Checks

Python: `python3 -m unittest tests.test_rdesim tests.test_runs tests.test_reaction_models tests.test_detailed_chemistry`.

1D wave speed against Chapman-Jouguet: `python3 -m rdesim.cjcheck` (see the
README). It is a sanity check of heat release and flux scheme at one resolution,
not a data validation and not a cell-size check.

What was run for `detailed` (OpenFOAM 14, Cantera 3.2 as the reference):

- 0D, `python3 -m rdesim.chemcheck` (`./bin/rdesim-wsl chemcheck`): a closed
  4-cell tube at one state, so only the chemistry step acts, against a Cantera
  constant-volume reactor on `source/Burke2012.yaml`. 1200 K, 1 atm, phi 1:
  ignition delay (maximum dT/dt) 49.5 microseconds for both, final T 2896.3 K
  against 2894.0 K (0.08 %), final p 0.06 %, largest T difference away from the
  ignition 3.9 K. 1000 K, 2 atm: 497 against 499 microseconds (write interval
  2 microseconds), final T 0.04 %. 1500 K, 1 atm: 14.25 for both, final T 0.12 %.
  The T difference of a few kelvin is not explained further here.
- 1D, `python3 -m rdesim.cjcheck --model detailed`: reference is the
  equilibrium CJ state of the same nine species from Cantera (1974.96 m/s,
  15.55 times the initial pressure, 2960 K; published values for stoichiometric H2-air
  at 1 atm and 300 K are close to 1970 m/s). Measured
  lead-shock speed 2005.3 m/s (1.0154 CJ) at 20 micrometre cells and 2004.3
  m/s (1.0149) at 60 micrometre cells. The shock to half-fuel-consumption
  distance was 77 micrometres at 20 and 155 micrometres at 60, so the reaction
  zone is not converged at these resolutions. Only the wave speed is a result.
- Unchanged elsewhere: with `reaction.model` absent the case, the results
  (`p`, `T`, `U`, `lambda`, `rho` of `annulus2D` at 5e-6 s) are byte-identical
  to `main`.
- Cost: `annulus2D-cells` (46800 cells) first 3 microseconds, one core: 174 s
  detailed, 82 s one-step. The 1D tubes above took 24.5 minutes (6000 cells)
  and 3 minutes (2000 cells) for about 56 microseconds of simulated time.

Not validated: any 2D structure or cell width, the burnt-gas hyperbolicity,
anything above the 3500 K janaf limit, and the first-order splitting error in
the presence of a steep front (the 1D result above is the only evidence).

## Status

| Model | Selectable | Compiled | Run |
|---|---|---|---|
| `oneStepProgress` | default | yes | `annulus2D`, 1D CJ tube |
| `twoStep` | `reaction.model` or `--model` | yes | `annulus2D`, 1D CJ tube |
| `detailed` | `reaction.model detailed` or `--model detailed` | yes | `annulus2D` and the fine `annulus2D-cells` channel (3 microseconds), 0D `chemcheck`, 1D CJ tube |
