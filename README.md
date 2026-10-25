# RDEsim

OpenFOAM runs in WSL. ParaView runs on Windows and opens the case through the WSL path.

Edit `cases/annulus2D/rdeConfig` and `fuels/<name>/fuelDict` in Cursor. Those two files are the case. There is no separate GUI. With the folder open in Cursor's WSL window, run a task instead of typing the shell commands: Command Palette → **Tasks: Run Task**.

| Task | What it does |
| --- | --- |
| **RDEsim: Prepare case** | Rewrites `cases/annulus2D` from `rdeConfig` |
| **RDEsim: Build solver** | Sources the OpenFOAM 14 bashrc and runs `Allwmake` |
| **RDEsim: Run annulus2D** | Sources that bashrc and runs `cases/annulus2D/Allrun` (`blockMesh`, then `foamRun`) |
| **RDEsim: Open latest result in Windows ParaView** | Opens `cases/annulus2D/annulus2D.foam` in the Windows ParaView install |

Use that order after you change `rdeConfig` or a fuel file. Build again only when the solver sources change. The same four actions from a WSL terminal are `./bin/rdesim-wsl prepare`, `build`, `run`, and `open`.

The open task finds ParaView at `/mnt/c/Program Files/ParaView*/bin/paraview.exe` and opens:

```text
\\wsl$\<distro>\home\<you>\RDEsim\cases\annulus2D\annulus2D.foam
```

`<distro>` comes from `$WSL_DISTRO_NAME` in the WSL terminal. You do not type the path. Colour **p** first, then **T**, **rho**, **U**, and **lambda**. **e** is there when the thermo model has written it. There is no separate detonation-indicator field.

The checkout includes the initial directory `0/`. A result time folder exists only after **RDEsim: Run annulus2D** (`Allrun`) finishes. The open task stops with an error if that folder is not there yet, and it does not treat `0/` as a result. Build and run look for the OpenFOAM 14 bashrc at `/opt/openfoam14/etc/bashrc`, then `~/OpenFOAM/OpenFOAM-14/etc/bashrc`, and stop with an error if neither file exists.

`paraFoam` inside WSL is optional and secondary. From `cases/annulus2D`, after `. /opt/openfoam14/etc/bashrc`, `paraFoam` uses the reader shipped with OpenFOAM. The Windows ParaView task is the one to use with ParaView installed on Windows.

RDEsim writes an OpenFOAM 14 case, checks that the selected model stays a hyperbolic balance law, and can repeat that case while it changes wall steps and other geometry. The flow solver is `rdeFluid`, a small module on top of the Foundation `shockFluid` solver (Kurganov–Tadmor).

This is not a production RDE code. The gas is one perfect gas with one progress variable. The fuel numbers are illustrative. Nothing here has been compared with a measured detonation speed, thrust, or wave count. See [What is not physically validated](#what-is-not-physically-validated).

The checkout is an ordinary git tree. Clone it on the Ubuntu distro and open that Linux path in Cursor's WSL window. `.vscode/tasks.json` is only the four tasks above. The machine that prepared the tree does not share your WSL files.

## Environment

Locked target:

- Ubuntu 26.04.1 LTS (Resolute Raccoon) in WSL2, where every OpenFOAM command runs. Not Ubuntu 26.10.
- OpenFOAM Foundation 14 (released 14 July 2026), package `openfoam14`, sourced with `. /opt/openfoam14/etc/bashrc` in that WSL terminal. Not the OpenFOAM.com / Keysight line and not v2606.
- ParaView installed on Windows. It opens `\\wsl$\<distro>\...\cases\annulus2D\annulus2D.foam`. A `paraFoam` session inside WSL is optional.

The install commands below are the ones published at <https://openfoam.org/download/14-ubuntu/>. They were not executed while this tree was prepared. `wmake`, `blockMesh`, `foamRun`, and ParaView were not run here either: this environment has no OpenFOAM and no ParaView. The solver sources follow the OpenFOAM 14 API on the Foundation git line (`shockFluid`, `fieldMapper`, `foamRun` loading `librdeFluidSolver.so`). Foundation v14 stays compatible inside the v14 line; a later pack such as 20261008 was not compiled against.

### Install OpenFOAM 14 on Ubuntu 26.04

In the Ubuntu terminal:

```sh
sudo sh -c "wget -O - https://dl.openfoam.org/gpg.key > /etc/apt/trusted.gpg.d/openfoam.asc"
sudo rm -f /etc/apt/sources.list.d/*dl_openfoam_org*list
sudo add-apt-repository "http://dl.openfoam.org/ubuntu main dev"
```

The repository URL is `http`, not `https`, as on the Foundation page. On Ubuntu 26.04 only, if `apt` warns `Skipping acquire of configured file 'main/binary-amd64v3/Packages'`, pin the architecture:

```sh
sudo sed -i "s/^\(deb\).*\(http\)/\1 [arch=amd64] \2/" /etc/apt/sources.list.d/*dl_openfoam_org*list
```

Then:

```sh
sudo apt update
sudo apt -y install openfoam14
```

Confirm the pack date with `apt-cache policy openfoam14`. The July 2026 packs and the later `20260724` / `20261008` packs are the same v14 line. An optional newer source tarball, not required for the package install, is `http://dl.openfoam.org/source/14.20261008`. The Foundation page says this package also installs a `paraview` package inside Ubuntu when one is missing. That Ubuntu binary is not the viewer for this setup. Open results in the ParaView installed on Windows, using the `\\wsl$\` path at the top of this file.

Every new terminal that builds or runs a case needs the Foundation bashrc:

```sh
. /opt/openfoam14/etc/bashrc
foamRun -help
```

`foamRun -help` only checks that the package is on `PATH`. It does not run RDEsim.

## Build

From the repository root, after sourcing the bashrc:

```sh
./Allwmake
```

`Allwmake` runs `wmake -s libso applications/modules/rdeFluid` and installs `librdeFluidSolver.so` under `$FOAM_USER_LIBBIN`. OpenFOAM 14 exports `FOAM_MODULES` to `$WM_PROJECT_DIR/applications/modules`. `./Allclean` removes that build. `bin/rdeFoam` is `foamRun -solver rdeFluid`. The case `controlDict` already says `solver rdeFluid;`, so `foamRun` from the case directory is enough.

This module was not compiled here.

## Layout

| Path | Role |
| --- | --- |
| `cases/annulus2D/rdeConfig` | Every case setting, in SI |
| `fuels/<name>/fuelDict` | Fuel, oxidizer, heat release, one-step rate |
| `rdesim/` | Config reader, mesh, checks, optimization driver |
| `applications/modules/rdeFluid/` | OpenFOAM 14 solver and acoustic wall condition |
| `paraview/view_rde.py` | `pvbatch` script for the OpenFOAM reader |
| `cases/annulus2D/` | Example case. `python3 -m rdesim prepare` overwrites the generated dictionaries |

Python is the standard library only. There is no install step.

```sh
python3 -m unittest tests.test_rdesim -v
python3 -m rdesim check --config cases/annulus2D/rdeConfig
```

## Change a setting

Edit `cases/annulus2D/rdeConfig` or the fuel file, then regenerate. Values are SI. The reader does not convert units. A missing required entry is an error. Nothing in the list below is filled in by a hidden constant in the solver: chemistry coefficients, geometry, and scheme names are read back from the dictionaries `prepare` writes.

```sh
python3 -m rdesim prepare --config cases/annulus2D/rdeConfig --case cases/annulus2D
```

`prepare` exits 1 when the hyperbolicity check fails. It still writes `constant/hyperbolicityReport` so you can read the reason. Hand edits under `system/`, `constant/`, and `0/` are replaced by the next `prepare`.

### Domain and mesh — `domain` and `mesh` in `rdeConfig`

`domain.dimension` is 2 or 3. `domain.topology` is `unwrapped` or `meridional` (3D ignores the 2D topology name and builds an annular sector). `RInner`, `ROuter`, `axialLength`, `sectorAngle`, `theta0`, `emptyThickness`, `axialReference`, and `unwrapRadius` (`mean`, `inner`, or `outer`) set the annulus. Unwrap length is `sectorAngle` in radians times that reference radius at the inlet, before steps.

`mesh` sets `nCircumferential`, `nAxial`, `nRadial`, `nCircumferentialBlocks`, the three grading entries, `minCellsPerBlock`, and `stationTolerance`. Grading is the same `simpleGrading` on every block, not one distribution stretched over the whole domain. Cell counts use the largest remainder after each block receives `minCellsPerBlock`.

`prepare` writes `system/blockMeshDict`. Missed faces go to a patch named `unassigned`.

### Geometry — `geometry` in `rdeConfig`

The wall is built in this order:

1. Wall angles about `domain.axialReference`: `R(y) = R0 + (y - yRef) tan(angle)`.
2. `widthProfile`, piecewise linear, constant outside the tabulated stations. `widthAnchor` is `inner`, `outer`, or `centre`.
3. `nozzle`, which overwrites the gap on `[start, start+length]`. `conical` uses `throatFraction` in `(0, 1)`, `throatGap`, and `exitGap`. `piecewise` needs at least two profile points.
4. `steps`, applied last. A step shifts every downstream station, including the nozzle. The example has both a step and a nozzle, so the throat after the upstream step is smaller than `throatGap`.

`steps[].height` greater than 0 moves that wall into the channel (outer radius decreases, inner radius increases). `angle` is in `(0, 90]`. 90 degrees is a vertical face named `outerStep` or `innerStep`. A smaller angle is a ramp of length `|height|/tan(angle)` and does not add a step face. A negative gap or `RInner >= ROuter` is an error.

`injector.nSlots` of 0 means the whole inlet is the injector. Otherwise the slots are centred at `phase + (i+0.5)*spacing` and must lie fully inside the sector. `slotWidth` is smaller than the spacing. `recessDepth` is upstream of `y = 0`. `injectionAngle` is degrees from axial toward circumferential.

What each mesh actually cuts:

- **2D unwrapped (the example).** `x` circumferential, `y` axial, `z` the empty thickness. `nRadial` is 1. Injector slots and the recess are meshed. Radial steps and the nozzle contour are written to `constant/wallProfile` and are not cut into the boundary.
- **2D meridional.** `x` axial, `y` radial, `z` empty. `nCircumferential` is 1. Radial steps, the nozzle, and an annular recess are meshed. Slots are not split.
- **3D.** `z` axial, `x = r cos θ`, `y = r sin θ`. Each circumferential block spans less than 180 degrees. A 360 degree sector closes without cyclics; a smaller sector uses `periodicLeft` and `periodicRight`. Slots are not split. The recess is an annular ring. `injectionAngle` and any circumferential inlet velocity must be 0, because one uniform vector cannot carry swirl.

`constant/geometryReport` repeats those limits for the case you generated. The ignition box is in mesh coordinates. The example kernel is `(0.002 0.01 -0.001)` to `(0.008 0.03 0.002)`.

### Operating conditions — `operating`

`pInlet`, `TInlet`, `pOutlet`, `TOutlet`, `equivalenceRatio`, `massFlowRate` (kg/s through the meshed sector), `plenumPressure`, `plenumTemperature`, and `inletModel`.

- `fixedVelocity` uses the three `inletVelocity` components.
- `massFlow` sets the speed to `massFlowRate / (rho * A)`. `A` is the physical annular sector area, including slot coverage: `sectorFraction * 1/2 * (ROuter^2 - RInner^2)` times the open fraction. On a 2D mesh that speed is the physical speed. The flux through `emptyThickness` is not the physical mass flow.
- `isentropicPlenum` uses the perfect-gas critical pressure ratio from the mixture gamma. `pInlet` above `plenumPressure` is an error.

`injectionAngle` sets the direction for `massFlow` and `isentropicPlenum`. `fixedVelocity` keeps the components you wrote.

### Fuel and chemistry — `fuel`, `fuelOverrides`, and `fuels/<name>/fuelDict`

`fuel.name` selects `fuels/<name>/fuelDict`. `fuel.library` is a path relative to the config file, or a directory name found by walking parents. `fuelOverrides` is merged on top of that file. The library has `hydrogen-air` and `ethylene-air`. Copy a directory, edit `fuelDict`, and set `fuel.name`.

The fuel file holds the fuel name and molar mass, the oxidizer species (mole fractions must sum to 1), `nuFuel` and `nuOxygen`, `heatRelease.LHV` in J per kg of fuel, `thermo` (`hConst`, `perfectGas`, `frozenReactants` or `frozenProducts`, gamma, hf, product molar mass and gamma), `transport.mu` and `Pr`, and `chemistry` (`A`, `beta`, `Ta`, `TRef`, `TFloor`, `ignitionTemperature`, `lambdaMin`, `lambdaMax`). Those chemistry numbers are not in the C++ or Python sources.

`constants.RUniversal` in `rdeConfig` is the gas constant in J/(mol·K). The example uses 8.314462618. Specific gas constant is `RUniversal / W`. `Cp = gamma/(gamma-1) * R_specific`. Heat release in the energy equation is `q = Y_f * LHV`. `hf` is passed through; the example is 0, so the heat release is `q`, not `hf`.

`Y_f,st` comes from the stoichiometry. `Y_f = phi * Y_f,st / (1 + (phi-1) Y_f,st)`. `frozenReactants` mixes fuel and oxidizer molar masses and uses `thermo.gamma`. `frozenProducts` uses `productMolWeight` and `productGamma`. The thermo model is still one gamma. A `janaf` thermo or any equation of state other than `perfectGas` is rejected.

`numerics.flowModel euler` writes `mu` as 0 in `constant/physicalProperties` even though the fuel file still holds the viscosity you would use for `navierStokes`.

### Numerics and hyperbolicity — `numerics` and `hyperbolicity`

`flowModel` is `euler` or `navierStokes`. `fluxScheme` is `Kurganov` or `Tadmor`. `ddtScheme`, `gradScheme`, `laplacianScheme`, `reconstructScalar`, `reconstructVector`, `fluxInterpolation`, `CFL`, `deltaT`, `maxDeltaT`, `startTime`, `endTime`, `adjustTimeStep`, `deltaTFactor`, `writeFormat`, `writePrecision`, `timePrecision`, `runTimeModifiable`, and `sourceIntegration` (`pointImplicit` or `explicit`) are all written through.

`fvSolution` is copied from `numerics.fvSolution`. The example names `rho`, `lambda`, `U`, and `e` exactly, because this dictionary reader treats a parenthesis as the end of a word, so a solver regex such as `"(U|e)"` is not used. OpenFOAM falls back from `lambdaFinal` to `lambda` when it looks for a final solver; that fallback was not runtime-checked. `PIMPLE.nOuterCorrectors` is in the same dictionary.

`momentumTransport` and `thermophysicalTransport` are copied when `flowModel` is `navierStokes` (the transport dictionary is required then). The example laminar Fourier entry was not runtime-checked. Euler mode does not read it.

`hyperbolicity` holds `requireStrictlyHyperbolic`, the floors `gammaMin`, `rhoMin`, `pMin`, `TMin`, `c2Min`, the stability warning level `maxCoStability`, `sourceDependsOnGradients`, `boundaryCoupling`, and the allowed flux, equation-of-state, and source lists.

`prepare` writes `constant/rdeProperties` with the chemistry, the floors, and the objective list. `rdeFluid` reads that file. It does not contain a second copy of the fuel numbers.

### Boundaries, acoustics, ignition, output

`boundaries.inlet` chooses `p`, `T`, and `U` types. The inflow patch is named `injector`. Lands are `inletWall`. Slot side walls are `injectorSide`. The outlet is `outlet`. `boundaries.outlet.type` is `waveTransmissive`, `fixedValue`, or `zeroGradient`. `waveTransmissive` gets `gamma` from the mixture, `fieldInf` from `pOutlet`, and `lInf` from the config.

`boundaries.walls.velocityModel` is `slip` or `acoustic`. Acoustic patches are `boundaries.acoustics.targetPatches`. A name that the mesh does not have is an error. `mode` is `impedance` (`ZReal`, `ZImag`), `reflection` (`RReal`, `RImag`), or `microstructure`. A reflection of `+1` cannot be inverted to an impedance; use a large `ZReal` instead. `timeDomainModel` is `slip`, `resistive`, or `helmholtz`. `frequency`, `incidenceAngleDeg` in `(-90, 90)`, `pRef`, and `resistiveFloor` are in the same dictionary.

Microstructure entries live under `boundaries.acoustics.microstructure`: `quarterWave` and `groove` (`cavityDepth`, `viscousResistance`, and for a groove `grooveWidth` and `grooveLand`) and `helmholtz` (`neckLength`, `neckRadius`, `cavityVolume`, `endCorrection`, `cellArea`, `normalizedResistance`, `initialDisplacement`, `initialVelocity`). The example keeps the unused numbers filled in so changing `type` does not require new keys. `endCorrection` is your number; the example is 1.7. `constant/acousticReflection` stores the complex impedance and reflection coefficient that `prepare` computed.

`ignition` sets `type` (`uniform`, `hotKernel`, `progressKernel`, `hotProgressKernel`), the initial `p`, `T`, `lambda`, `inletLambda`, `initialVelocity`, the kernel box, `TKernel`, and `lambdaKernel`.

`output.fields`, `writeControl`, `writeInterval`, `purgeWrite`, and `writeCompression` control what `foamRun` writes. `output.objectives` lists scalars. A patch objective is an `areaAverage`. A volume objective is `max`, `min`, or `average`. The example outlet average is patch `outlet`. The inlet average is patch `injector`.

## Run a case

In the WSL Ubuntu terminal, from the repository root:

```sh
. /opt/openfoam14/etc/bashrc
cd cases/annulus2D
./Allrun
```

`Allrun` is `blockMesh` then `foamRun`. `./Allclean` removes `constant/polyMesh`, `postProcessing`, numbered VTK directories, `log.*`, and `processor*`. It leaves `VTK/geometryBlocks.vtk` and the initial `0/` directory. Result time folders are created by `foamRun` and exist only after `Allrun`.

The example `endTime` is `5e-6` s and `deltaT` is `1e-8` s. That is a short first window, not a burnout.

After a run, `rdeFluid` writes `postProcessing/rdeObjectives/objectives.dat` (latest) and a copy under the time name. It also writes `postProcessing/hyperbolicity/report` at start and stops if that check fails. The Python report in `constant/hyperbolicityReport` is the same kind of check before the solver starts.

## Optimization

`optimization` in `rdeConfig` sets `method` (`grid` or `coordinate`), `maxEvaluations`, `timeoutSeconds`, `coordinateLevels`, `workDirectory`, and `solverCommand`. The example command is `foamRun`. The command is split on whitespace and is not passed to a shell. Each evaluation is a case directory under `workDirectory`, and the command runs there.

`objective.expression` may use `+ - * /`, parentheses, and names from `output.objectives`. The example is `pMeanOutlet/pMeanInlet`, sense `maximize`. `designVariables` name an existing config path such as `geometry.steps[0].height`. `samples` includes both ends; one sample is the minimum only. `grid` is the Cartesian product. `coordinate` starts at the first sample of each variable, sweeps one variable at a time, keeps the best feasible point, and repeats for `coordinateLevels`.

A constraint whose quantity contains `.` or `[` is a config value checked before the run. Any other quantity must appear in `objectives.dat` after the run. Infeasible and failed evaluations are recorded. `python3 -m rdesim optimize --config cases/annulus2D/rdeConfig` exits 1 when no feasible point remains. The summary is `workDirectory/summary`.

The unit test drives the same loop with `tests/stub_solver.py` so it does not need OpenFOAM. The example command stays `foamRun`.

## ParaView on Windows

Use **RDEsim: Open latest result in Windows ParaView** after the run task. That task opens the case stub below in the ParaView installed on Windows. To open it yourself, File → Open and choose the same path. Replace `<distro>` with a name printed by `wsl.exe -l -q`. Replace `<you>` with the Ubuntu home directory:

```text
\\wsl$\<distro>\home\<you>\RDEsim\cases\annulus2D\annulus2D.foam
```

That file is the OpenFOAM reader stub. It is in the checkout before any run. The mesh and fields come from the case next to it. Colour these cell fields, in this order when you want to see the wave:

1. **p** — pressure
2. **T** — temperature
3. **rho** — density
4. **U** — velocity
5. **lambda** — progress variable (`0` unburned, `1` burned)

The initial directory `0/` is already in the checkout. A later time folder, the one you step to after the start, exists only after **RDEsim: Run annulus2D** (`Allrun`) in WSL. The open task will not start ParaView until that folder is there. **e** (sensible internal energy) is there when the solver has written it. There is no extra detonation-indicator field. The wave, if the short window produces one, is the jump in **p**, **T**, **rho**, and **lambda**.

`VTK/geometryBlocks.vtk` is one hex per block (VTK cell type 12). It is the block outline, not the mesh `blockMesh` builds. Use `annulus2D.foam` for the fields.

`paraFoam` inside WSL is a secondary option, only if you want the reader that ships with OpenFOAM. From `cases/annulus2D`, after `. /opt/openfoam14/etc/bashrc`:

```sh
paraFoam
```

`pvbatch paraview/view_rde.py cases/annulus2D` is the same kind of optional WSL path. It was not run here. There is no `.pvsm` state file, because those files embed absolute paths. The Windows `\\wsl$\` path above is the setup this repository expects.

## Equations

Default model, `flowModel euler`: the reactive Euler equations, a hyperbolic balance law. `u` is velocity, `rho` density, `p` pressure, `e` sensible internal energy, `lambda` the progress variable, `omega` the progress rate.

```text
d(rho)/dt + div(rho u) = 0
d(rho u)/dt + div(rho u u + p I) = 0
d(rho e)/dt + div((rho e + p) u) = rho omega q
d(rho lambda)/dt + div(rho lambda u) = rho omega
```

The momentum and energy fluxes are the Kurganov–Tadmor fluxes from `shockFluid`, including the kinetic coupling that solver already applies. `p = rho R_specific T` with `R_specific = RUniversal / W`. The eigenvalues of the inviscid flux are `u-c`, `u`, and `u+c` when `c^2 = gamma p / rho > 0`.

`omega = A (T/TRef)^beta exp(-Ta/T) (lambdaMax - lambda)` when `T` is at least `ignitionTemperature`, and 0 otherwise. `TFloor` only keeps the exponential argument away from zero. `pointImplicit` treats the `lambda` sink implicitly. `explicit` uses `max(lambdaMax - lambda, 0)`. `lambda` is clamped to `[lambdaMin, lambdaMax]` after the solve. The energy source uses the same `omega` and adds `rho omega q`.

`navierStokes` turns on viscosity and the thermal flux. Those terms make the system incompletely parabolic. With `requireStrictlyHyperbolic yes` (the example), that choice fails the check. The Euler path forces `mu = 0`.

The check fails when gamma is not above `gammaMin`, when density, pressure, temperature, or `c^2` are not above their floors, when the equation of state, flux, or source integration is outside your list or outside what this version implements (`perfectGas`, `Kurganov`/`Tadmor`, `pointImplicit`/`explicit`), when `sourceDependsOnGradients` is yes, when `boundaryCoupling` is not `local`, when a strictly hyperbolic run asks for Navier–Stokes or `steadyState`, or when the acoustic time-domain model is not `slip`, `resistive`, or `helmholtz`. `CFL` above `maxCoStability` is a stability warning, not a loss of hyperbolicity. The Python check also rejects gamma at or below 1 while it builds the mixture, before the report.

## Acoustics

The wall is locally reacting, with the `exp(+i ω t)` convention. Normal incidence is the default. Oblique reflection is `(ζ cos θ - 1) / (ζ cos θ + 1)` with `ζ = Z / (rho c)`.

- Impedance uses `ZReal` and `ZImag`.
- Reflection inverts `R` to `Z`. `R = +1` is rejected.
- A quarter-wave tube is `Z = R_v - i rho c cot(k L)` with `k = 2 π f / c` and open-area ratio 1.
- A groove uses `sigma = width / (width + land)` and the same closed pipe.
- A Helmholtz resonator uses neck area `π neckRadius^2`, `l_eff = neckLength + endCorrection*neckRadius`, `sigma = neckArea/cellArea` in `(0, 1]`, and `Z_H = rho c [R_v + i (k l_eff - neckArea/(k V))]`, then `Z_w = Z_H/sigma`. The resonance frequency is `c/(2π) sqrt(neckArea/(l_eff V))`.

The time-domain condition `acousticWallVelocity` is local:

- `slip` copies the tangential interior velocity and sets the normal component to 0.
- `resistive` sets `u_n = (p - pRef) / ZReal` and stops if `|ZReal|` is below `resistiveFloor`. A near-zero real impedance is a pressure-release wall.
- `helmholtz` integrates `m ξ'' + r ξ' + s ξ = p'` with explicit Euler and `u_face = sigma ξ'`. The coefficients `m`, `r`, `s`, and `sigma` are the `inertia`, `resistance`, `stiffness`, and `porosity` written into `0/U`.

Quarter-wave and groove reactive impedances are reported in `constant/acousticReflection`. The time-domain condition does not enforce `Im(Z)` for those shapes. A Helmholtz restart does not restore the per-face oscillator state; `ξ` and `ξ'` start again from `initialDisplacement` and `initialVelocity`.

A hard wall is `|R| → 1` for a large real `Z`. A soft wall is `R = -1` at `Z = 0`. `Z = rho c` at normal incidence gives `R = 0`. The same `Z` at 60 degrees gives `R = -1/3`.

## What is not physically validated

- One-step, single-gamma, perfect-gas kinetics. Not a detailed mechanism, and not a variable-gamma product gas.
- Hydrogen and ethylene entries are starting numbers, not a calibrated detonation chemistry.
- No comparison with Chapman–Jouguet speed, cell size, wave count, or thrust.
- The 2D unwrapped mesh does not resolve radial steps or the nozzle wall.
- The 3D mesh does not resolve discrete injector slots.
- Acoustics are lumped. Cavities are not meshed. The Helmholtz ODE is explicit Euler and is not a proven stable integrator for every `deltaT`.
- The example end time is a smoke window.
- `blockMesh` face order, 3D arc edges, the zonal kernel syntax, the laminar Fourier dictionary, and this solver were not executed against OpenFOAM 14 here.
- A CFL warning is about the explicit flux, not about whether a detonation will form.

## License

GPL-3.0-or-later, because `rdeFluid` links OpenFOAM. See `LICENSE`.
