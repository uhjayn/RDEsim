# RDEsim

OpenFOAM runs in WSL. ParaView runs on Windows and opens the case through the WSL path.

RDEsim writes an OpenFOAM 14 case for a rotating-detonation-engine annulus, checks that the selected model stays a hyperbolic balance law, and can repeat the case while it changes wall steps and other geometry. The flow solver is `rdeFluid`, a small module on top of the Foundation `shockFluid` solver (Kurganov–Tadmor).

This is not a production RDE code. The model is the reactive Euler equations with one perfect gas. The default reaction model has one progress variable, and a two-step model with an induction variable can be selected (see [Reaction models](#reaction-models)). The fuel numbers are illustrative. Nothing here has been validated against measured detonation data (speed, cell size, wave count, or thrust); the only wave-speed check compares the model with its own Chapman–Jouguet speed. See [What is not physically validated](#what-is-not-physically-validated).

## Quick start

You need Windows with WSL2 and an Ubuntu distro (its name shows as `Ubuntu-26.04` in `wsl.exe -l -q`), OpenFOAM Foundation 14 installed in that distro ([install steps](#install-openfoam-14-on-ubuntu-2604)), ParaView 6.2.0 for Windows, and Cursor with its WSL extension.

1. **Get the repo into WSL.** In the Ubuntu terminal:

   ```sh
   cd ~
   git clone https://github.com/uhjayn/RDEsim.git
   cd RDEsim
   git switch main
   git pull
   ```

   For an existing clone, run only `cd ~/RDEsim`, `git switch main`, and `git pull`. Then open `~/RDEsim` in Cursor through the WSL extension, so the folder is the Linux path and not a Windows copy.
2. **Build and run.** Command Palette → **Tasks: Run Task** → **RDEsim: Run and visualize**. That one task runs `./bin/rdesim-wsl simulate`: it prepares `cases/annulus2D`, builds the solver on the first run (or when the solver sources change), runs `blockMesh` and `foamRun`, then starts ParaView. The short case ends at `5e-6` s.
3. **Open the result.** The task prints `Windows reports paraview.exe is running.` only after Windows confirms it with `tasklist.exe`. If it does not, open ParaView yourself and use File → Open with the path the script prints, `\\wsl$\<distro>\home\<you>\RDEsim\cases\annulus2D\annulus2D.foam`. [ParaView on Windows](#paraview-on-windows) has the details.
4. **Colour a field.** In ParaView select `annulus2D.foam` in the Pipeline Browser and click Apply (choose the OpenFOAM reader if ParaView asks). Use Last Frame to go to the result time. In the Coloring drop-down pick **p**, then **T**, **rho**, **U**, and **lambda**.

For a longer window run **RDEsim: Run and visualize (long)** instead. Its limit is described under [Run a case](#run-a-case).

## Cursor tasks

Edit `cases/annulus2D/rdeConfig` and `fuels/<name>/fuelDict` in Cursor. Those two files are the case. There is no separate GUI. Then start **RDEsim: Run and visualize** from Command Palette → **Tasks: Run Task**. That one task runs `./bin/rdesim-wsl simulate`. The script sources the OpenFOAM 14 bashrc, prepares `cases/annulus2D` from `rdeConfig`, builds with `Allwmake` only when `librdeFluidSolver.so` is missing or the solver sources are newer, runs `Allrun`, and asks Windows to start ParaView. It stops at the first failure. It says ParaView is running only when `tasklist.exe` lists `paraview.exe`. If that launch fails, the case stays on disk and the script prints a path to paste into ParaView's File → Open.

The same command from a WSL terminal is `./bin/rdesim-wsl simulate`.

**RDEsim: Run and visualize (long)** runs `./bin/rdesim-wsl simulate --case cases/annulus2D-long`. That case is the same geometry and the same model as `cases/annulus2D`. Only the run window is longer, so a wave can move along the annulus. The short case stays the smoke test.

**RDEsim: Run validation** runs `./bin/rdesim-wsl validate`, the H2-air detonation-tube comparison, see [Validation](#validation-h2-air-detonation-tube). **RDEsim: Run validation (full, resolution study)** runs `validate-full`.

Four separate tasks remain if you want one step only: **RDEsim: Prepare case**, **RDEsim: Build solver**, **RDEsim: Run annulus2D**, and **RDEsim: Open latest result in Windows ParaView**. From a terminal those are `./bin/rdesim-wsl prepare`, `build`, `run`, and `open`. `prepare`, `run`, and `open` take the same `--case` option. `build` and `test` do not. **RDEsim: Cell structure view** runs `./bin/rdesim-wsl cells`.

Further tasks, each a thin wrapper over `./bin/rdesim-wsl`:

- **RDEsim: Run and visualize (run folder, resumes)** is `simulate --runs`. See [Run folders, restart, and checkpoints](#run-folders-restart-and-checkpoints).
- **RDEsim: Run and visualize (choose reaction model)** is `simulate --runs --model NAME` and asks for `oneStepProgress` or `twoStep`. See [Reaction models](#reaction-models).
- **RDEsim: 1D CJ check (choose reaction model)** is `cjcheck --model NAME`.
- **RDEsim: Cell structure view (fine)** is `cells --case cases/annulus2D-cells`. See [Cell structure and mesh resolution](#cell-structure-and-mesh-resolution).
- **RDEsim: Optimize (resume, with budget)** is `optimize --budget N` and asks for the budget. See [Optimizer budget and resume](#optimizer-budget-and-resume).

Options after the command: `--case PATH` (not with `build` or `test`), `--model NAME`, `--runs`, `--fresh`, and `--budget N` (`optimize` only). `./bin/rdesim-wsl` with no command prints the usage.

**RDEsim: Run tests** runs `./bin/rdesim-wsl test`. That sources the OpenFOAM 14 bashrc, runs `Allwmake`, runs `python3 -m unittest tests.test_rdesim tests.test_runs tests.test_reaction_models`, then `blockMesh` and `foamRun` on `cases/annulus2D`. It removes an earlier result time folder first and keeps `0/`. It stops at the first failure and prints the time folder `foamRun` wrote. `sh tests/test_paraview_launch.sh` checks the ParaView launch messages with stub `powershell.exe` and `tasklist.exe`. It does not start a real ParaView.

The scripts look for the OpenFOAM 14 bashrc at `/opt/openfoam14/etc/bashrc`, then `~/OpenFOAM/OpenFOAM-14/etc/bashrc`, and stop if neither file exists.

The checkout is an ordinary git tree. Clone it on the Ubuntu distro and open that Linux path in Cursor's WSL window. `.vscode/tasks.json` holds all the **RDEsim:** tasks listed above. The machine that prepared the tree does not share your WSL files.

## Environment

Target setup:

- Ubuntu 26.04.1 LTS (Resolute Raccoon) in WSL2, where every OpenFOAM command runs. Not Ubuntu 26.10. The distro name shown by `wsl.exe -l -q` is `Ubuntu-26.04`.
- OpenFOAM Foundation 14 (released 14 July 2026), package `openfoam14`, sourced with `. /opt/openfoam14/etc/bashrc` in that WSL terminal. Not the OpenFOAM.com / Keysight line and not v2606.
- ParaView 6.2.0 installed on Windows. The script starts it in the Windows session and opens `\\wsl.localhost\<distro>\...\cases\annulus2D\annulus2D.foam`. A `paraFoam` session inside WSL is optional.
- Cursor connected to the distro through its WSL extension, so the tasks run in Ubuntu.

The install commands below are the ones published at <https://openfoam.org/download/14-ubuntu/>. With OpenFOAM 14 installed, `./Allwmake` builds the solver and `blockMesh` plus `foamRun` run `cases/annulus2D`. `./bin/rdesim-wsl test` repeats that. The Windows ParaView launch is covered by `tests/test_paraview_launch.sh` with stub Windows tools only. The solver is the OpenFOAM 14 `shockFluid` module loaded by `foamRun` as `librdeFluidSolver.so`.

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

Confirm the pack date with `apt-cache policy openfoam14`. The July 2026 packs and the later `20260724` / `20261008` packs are the same v14 line. An optional newer source tarball, not required for the package install, is `http://dl.openfoam.org/source/14.20261008`. The Foundation page says this package also installs a `paraview` package inside Ubuntu when one is missing. That Ubuntu binary is not the viewer for this setup. Open results in the ParaView installed on Windows. The fallback is the `\\wsl$\` (or `\\wsl.localhost\`) path under [ParaView on Windows](#paraview-on-windows).

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

## Layout

| Path | Role |
| --- | --- |
| `cases/annulus2D/rdeConfig` | Every case setting, in SI |
| `fuels/<name>/fuelDict` | Fuel, oxidizer, heat release, one-step rate |
| `bin/rdesim-wsl` | Script behind the Cursor tasks: simulate, prepare, build, run, open, test |
| `.vscode/tasks.json` | The **RDEsim:** Cursor tasks |
| `rdesim/` | Config reader, mesh, checks, run folders (`runs.py`), mesh-resolution report (`resolution.py`), 1D CJ check (`cjcheck.py`), optimization driver |
| `applications/modules/rdeFluid/` | OpenFOAM 14 solver and acoustic wall condition |
| `paraview/view_rde.py` | `pvbatch` script for the OpenFOAM reader |
| `cases/detonationTubeH2/` | 1D H2-air validation: `validationConfig`, `referenceData` |
| `fuels/hydrogen-air-validation/` | Hydrogen fuel calibrated to the reference, written by `rdesim.validation calibrate` |
| `rdesim/validation.py`, `rdesim/detonation.py` | Validation driver and steady ZND/CJ theory |
| `tools/validation_reference.py` | Regenerates `referenceData` with Cantera and CEA |
| `cases/annulus2D/` | Example case. `python3 -m rdesim prepare` overwrites the generated dictionaries |
| `cases/annulus2D-long/` | Same case with a longer run window |
| `cases/annulus2D-cells/` | Fine-mesh sector, about 10 cells across the one-step chemical length |
| `runs/` | Run folders made by `--runs` and the 1D CJ tube. Ignored by git |
| `docs/reaction-models-design.md` | Design note for the reaction models |
| `tests/` | Unit tests and the stub-based ParaView launch test |

Python is the standard library only. There is no install step.

```sh
python3 -m unittest tests.test_rdesim tests.test_runs tests.test_reaction_models -v
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

Fuel molar masses stay in kg/mol. OpenFOAM 14's `specie` entry is kg/kmol and `R() = 8314.462618/molWeight`. `prepare` writes that kg/kmol value so OpenFOAM's gas constant equals `RUniversal / W`. With the example `RUniversal`, that is 1000 times the SI molar mass. `Cp` stays in J/(kg·K).

`Y_f,st` comes from the stoichiometry. `Y_f = phi * Y_f,st / (1 + (phi-1) Y_f,st)`. `frozenReactants` mixes fuel and oxidizer molar masses and uses `thermo.gamma`. `frozenProducts` uses `productMolWeight` and `productGamma`. The thermo model is still one gamma. A `janaf` thermo or any equation of state other than `perfectGas` is rejected.

`numerics.flowModel euler` writes `mu` as 0 in `constant/physicalProperties` even though the fuel file still holds the viscosity you would use for `navierStokes`.

### Numerics and hyperbolicity — `numerics` and `hyperbolicity`

`flowModel` is `euler` or `navierStokes`. `fluxScheme` is `Kurganov` or `Tadmor`. `ddtScheme`, `gradScheme`, `laplacianScheme`, `reconstructScalar`, `reconstructVector`, `fluxInterpolation`, `CFL`, `deltaT`, `maxDeltaT`, `startTime`, `endTime`, `adjustTimeStep`, `deltaTFactor`, `writeFormat`, `writePrecision`, `timePrecision`, `runTimeModifiable`, and `sourceIntegration` (`pointImplicit` or `explicit`) are all written through.

`fvSolution` is copied from `numerics.fvSolution`. The example names `rho`, `lambda`, `U`, and `e` exactly, and repeats them as `rhoFinal`, `lambdaFinal`, `UFinal`, and `eFinal` with the same solver settings. OpenFOAM 14 looks up the `Final` entry on the last outer corrector and does not fall back to the name without `Final`. A solver regex such as `"(U|e)"` is not used, because this dictionary reader treats a parenthesis as the end of a word. `PIMPLE.nOuterCorrectors` is in the same dictionary. The example uses 2.

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

The smoke-test `endTime` is `5e-6` s and `deltaT` is `1e-8` s. That is a short first window, not a burnout. On that window the disturbance stays next to the hot kernel.

`cases/annulus2D-long` is the same geometry and the same model. `endTime` is `2e-5` s and `writeInterval` is `2e-6` s, so the run writes 10 time folders. The annulus is `0.12` m long and the inlet sound speed is about `409` m/s. In the smoke case, pressure above `1.1` times the inlet value only reaches about `y = 0.035` m. At `2e-5` s that front is at about `y = 0.048` m. Colour **p** and step the time slider from `2e-06` to `2e-05`. **U** moves with the same front; **p** is the field that was checked. The shipped long case stops at `2e-5` s and has `positivity` off. Without the guard, a run past about `2.4e-5` s makes a cell temperature non-positive, and `shockFluid` then hits a floating-point exception in the sound-speed square root. Tightening the Courant number to `0.05` and switching the reconstruction to `vanLeer` did not move that failure. The `positivity` switch does; see [Temperature positivity](#temperature-positivity). On the authoring machine `blockMesh` and `foamRun` finished this window in 2 seconds of wall time. `foamRun` reported an execution time of 1.1 seconds. The maximum Courant number stayed at about 0.20. Temperatures and densities stayed positive, and there were no NaNs.

`cases/annulus2D-cells` is a straight 5 degree sector, 5.6 mm long, with 180 by 260 cells (46800). The nozzle and the wall step are off, and there is no injector recess, so the spacing stays uniform. At the 2500 K kernel the one-step chemical length is the inlet sound speed divided by the rate, not a ZND induction length. That length is 0.22 mm. The edges are about 0.022 mm, so the count is about 10.4 axial cells and 10.2 circumferential cells across it. That count is not a validated cell size. `numerics.startFrom` is `latestTime`, `writeInterval` is `5e-7` s, and `endTime` is `3e-6` s. A later `foamRun` continues from the last time folder. Run `blockMesh` again only when the mesh changes. **RDEsim: Cell structure view (fine)** runs `./bin/rdesim-wsl cells --case cases/annulus2D-cells`. `positivity` is on, with `TMin` 50.

On the authoring machine `blockMesh` and `foamRun` reached `3e-6` s in 81 s of wall time. `foamRun` reported an execution time of 80 s and exited 0. The maximum Courant number was 0.202. No field went NaN, and density stayed positive. The lowest temperature was 225 K, so the positivity guard recorded 0 cell-steps. At `3e-6` s the pressure above 1.1 times the inlet value reaches about `y = 0.0029` m. `pMax` is one bright patch at the injector end, the hot kernel spread across part of the sector. Circumferential slices of that patch are one hump, not a repeated cell pattern. That picture is a pressure history. It is not a validated cell size. `prepare` writes the `peakPressure` function object from `output.cellStructure`, on the same `5e-7` s interval as the flow fields.

After a run, `rdeFluid` writes `postProcessing/rdeObjectives/objectives.dat` (latest) and a copy under the time name. It also writes `postProcessing/hyperbolicity/report` at start and stops if that check fails. The Python report in `constant/hyperbolicityReport` is the same kind of check before the solver starts.

## Optimization

`optimization` in `rdeConfig` sets `method` (`grid` or `coordinate`), `maxEvaluations`, `timeoutSeconds`, `coordinateLevels`, `workDirectory`, and `solverCommand`. The example command is `./Allrun`, because each evaluation is a new case directory with no mesh yet. `Allrun` is `blockMesh` then `foamRun`. The command is split on whitespace and is not passed to a shell, so `blockMesh && foamRun` is not a valid `solverCommand`. The command runs in that case directory.

`objective.expression` may use `+ - * /`, parentheses, and names from `output.objectives`. The example is `pMeanOutlet/pMeanInlet`, sense `maximize`. `designVariables` name an existing config path such as `geometry.steps[0].height`. `samples` includes both ends; one sample is the minimum only. `grid` is the Cartesian product. `coordinate` starts at the first sample of each variable, sweeps one variable at a time, keeps the best feasible point, and repeats for `coordinateLevels`.

A constraint whose quantity contains `.` or `[` is a config value checked before the run. Any other quantity must appear in `objectives.dat` after the run. Infeasible and failed evaluations are recorded. `python3 -m rdesim optimize --config cases/annulus2D/rdeConfig` exits 1 when no feasible point remains. The summary is `workDirectory/summary`.

The unit test drives the same loop with `tests/stub_solver.py` so it does not need OpenFOAM. The example command stays `foamRun`.

### Optimizer budget and resume

Every evaluation lives in `workDirectory/eval-<fingerprint>` where the fingerprint hashes the physics of that evaluation's config (not `endTime`, `output`, or `optimization`). A finished evaluation writes `evaluation.json`. Running `optimize` again reuses finished evaluations and only runs the missing ones, so an interrupted or budget-limited search continues where it stopped. `--budget N` (or `optimization.maxNewEvaluations`) caps how many new solver runs one invocation may start; reused evaluations are free. When the budget stops the search the summary lists the deferred points and the exit status is 0 if a feasible best exists. `--fresh` ignores the stored evaluations. From the terminal: `./bin/rdesim-wsl optimize [--case PATH] [--budget N] [--fresh]`, or the Cursor task **RDEsim: Optimize (resume, with budget)**, which asks for the budget. Setting `optimization.resume` to `false` does the same.

## Run folders, restart, and checkpoints

`python3 -m rdesim prepare --config cases/annulus2D/rdeConfig --run-folder` (or `./bin/rdesim-wsl simulate --runs`, `prepare --runs`, `run --runs`, `open --runs`, `cells --runs`; all take `--case` too, for example `simulate --case cases/annulus2D-long --runs`; the Cursor task **RDEsim: Run and visualize (run folder, resumes)** does this for the default case; add `--fresh` to delete the earlier result times first) writes a self-contained case to `runs/<case>-<model>-<fingerprint>/`. The fingerprint hashes the physics (geometry, mesh, operating point, fuel, ignition, boundaries, numerics except `endTime`). Re-running the same config reuses the folder and sets `startFrom latestTime`, so an interrupted run or a longer `endTime` continues from the last written time. A changed physics setting makes a new folder, so earlier results are never overwritten. `python3 -m rdesim run-dir --config <rdeConfig>` prints the folder path. `runs/` is ignored by git.

`output.writeInterval` is the normal write cadence. Pick a value that divides `endTime`: OpenFOAM writes only on that cadence, so the final state is lost if `endTime` is not a multiple of it. `output.checkpoint` (`{"enabled": true, "writeControl": "clockTime", "writeInterval": 600}`; `timeStep`, `runTime`, `cpuTime` also work) adds a `writeObjects` function object that writes restart fields at that cadence in addition. A checkpoint-only time folder has the solution fields and restarts correctly; it need not have every output field. `output.startFrom` accepts `startTime` (default), `latestTime`, or `firstTime`. `numerics.startFrom` (used by the fine-channel case) takes the same values and wins when both are set. Neither changes the run folder fingerprint.

## Reaction models

`reaction.model` in `rdeConfig` (or `--model NAME` on `prepare`, `check`, `optimize`, `./bin/rdesim-wsl`, or the Cursor task **RDEsim: Run and visualize (choose reaction model)**) selects the source terms. Without a `reaction` dictionary the model is `oneStepProgress`, the original single progress variable, and the case files and results are unchanged.

```
reaction
{
    model   twoStep;    // oneStepProgress | twoStep | detailed
}
```

`twoStep` adds an induction variable `xi` (0 to 1, no heat) followed by heat release `lambda` that starts where `xi = 1`. Its five numbers are in `fuels/<name>/fuelDict` under `twoStep` (`inductionA`, `inductionTa`, `releaseA`, `releaseTa`, `releaseOrder`); override them with `fuelOverrides`. The `hydrogen-air` values come from a fit to Burke 2012 ignition delays (20 bar, 1400 to 1900 K, within 7 %) and are not a validated detonation-structure model; the comment in that file lists exactly what was fitted. Cost: one extra transported scalar, so roughly the one-step cost per time step. It resolves a much shorter induction zone than the one-step example, so a cell size that is fine for one-step may be too coarse for it.

`detailed` carries the nine H2-air species of the Burke et al. (2012) H2/O2 mechanism (H, H2, O, OH, H2O, O2, HO2, H2O2, N2; 23 reactions, NASA-7 thermo, all numbers copied from the published data, see `docs/reaction-models-design.md` for source and the one deliberate change to N2's lower temperature limit). The species are advected with the same shock-capturing flux and a stiff ODE solver (`seulex`) advances the chemistry of every cell at constant volume after each transport step (first-order Lie splitting). Select it with `reaction.model detailed` or `--model detailed`; the Cursor task **RDEsim: Run and visualize (detailed chemistry)** runs `simulate --runs --model detailed`. Settings are in the `detailed{}` block of the fuelDict: `mechanism` (a directory with `thermo` and `reactions`), `TMin` (cells colder than this skip chemistry, 500 K), `initialChemicalTimeStep`, `ode{solver absTol relTol}`, and `log yes` to print the largest difference between the ODE temperature and the thermo temperature each step. Fields `H`, `H2`, ... `N2` are written next to `lambda`, which stays as `1 - Y_fuel/Y_fuel,unburned` so objectives and the cell-structure view keep working. The gas is no longer one `gamma`: its equilibrium CJ speed is the real one (1975 m/s, Cantera), not the 2627 m/s of the reduced models.

Cost and limits of `detailed`: Euler only (the mechanism files carry no transport data, so `navierStokes` is refused), no `progressKernel` ignition, no NOx, perfect-gas equation of state. On the fine `annulus2D-cells` channel (46800 cells) the first 3 microseconds took 174 s against 82 s for the one-step model on the same machine, one core; the cost grows with the number of hot cells because the ODE runs in every cell above `TMin`, so a long run with a large burnt region is several times more expensive than that. `foamRun` runs one core; detailed runs are the case for `--runs` (resume) and `output.checkpoint`. The mesh-resolution report of `cells` is based on the one-step rate and says so for other models. What was checked, all on OpenFOAM 14:
- `./bin/rdesim-wsl chemcheck`: one closed cell at 1200 K, 1 atm, phi 1 against Cantera's constant-volume reactor on the same mechanism. Ignition delay 49.5 microseconds for both, final T 2896 K against 2894 K (0.08 %). At 1000 K and 2 atm the delays were 497 and 499 microseconds and at 1500 K 14.25 microseconds for both (the write interval limits the resolution to 1 to 2 microseconds).
- `./bin/rdesim-wsl cjcheck --model detailed` needs `pip install cantera` for the reference. In a 0.12 m tube ignited by a 2500 K, 4 MPa kernel at 1 atm, 300 K, phi 1, the lead shock ran at 2005 m/s at 20 micrometre cells (6000 cells, 24.5 minutes on one core) and 2004 m/s at 60 micrometre cells (2000 cells, 3 minutes), against 1975 m/s for the equilibrium CJ speed of the same species, so 1.5 % above it. The wave is a driven one that has not necessarily relaxed to CJ, so 1.5 % above is not an error estimate. The distance from the shock to half fuel consumption was 77 micrometres (3.8 cells) at 20 micrometres and 155 micrometres at 60 micrometres, so the reaction zone is not converged at these resolutions: only the wave speed is a result, not a structure.
- Not validated: cell size or regularity, any 2D structure, the hyperbolicity of the burnt composition (checked for the unburned gas only), and anything above 3500 K (the janaf range limit of most species).

Limits that apply to `oneStepProgress` and `twoStep`: one perfect gas with one `gamma` and a constant `q`. For the example mixture the model's own Chapman-Jouguet speed is 2627 m/s; the real H2-air value is about 1975 m/s (see `docs/reaction-models-design.md`). They describe a hotter, faster wave than the real mixture.

Sanity check of the heat release and flux scheme in 1D: `./bin/rdesim-wsl cjcheck --model twoStep` (or `python3 -m rdesim.cjcheck --config cases/annulus2D/rdeConfig --model twoStep --case runs/cjcheck`). It writes a straight tube, ignites one end, runs `foamRun` and compares the lead-shock speed with the perfect-gas Chapman-Jouguet speed. With the defaults (6000 cells of 20 micrometres, 0.12 m) it takes about 2.5 minutes. It is a self-consistency check at one resolution, not a validation against measurements and not a cell-size check.

## ParaView on Windows

**RDEsim: Run and visualize** asks Windows to start ParaView after the run. The separate task **RDEsim: Open latest result in Windows ParaView** does only that launch. Both print the `paraview.exe` path, the exact `Start-Process` or `cmd.exe /c start` command, and its exit code. After a short wait they run `tasklist.exe /FI "IMAGENAME eq paraview.exe"`. The line `Windows reports paraview.exe is running.` is printed only when that list contains `paraview.exe`. A WSL pid is not a running ParaView.

The script tries the path from `wslpath -w` first, then `\\wsl.localhost\<distro>\...`, then `\\wsl$\<distro>\...`, and stops at the first one that produces a Windows `paraview.exe` process. A folder named `ParaView 6.2.0` matches. Set `RDESIM_PARAVIEW_EXE` if the executable is not under `Program Files`, `Program Files (x86)`, or `AppData\Local`.

If Task Manager still has no `paraview.exe`, the script says the launch failed, prints the command it ran, and prints this fallback (first the `\\wsl.localhost\` form, then the `\\wsl$\` form). If `paraview.exe` starts and no window appears, reinstall ParaView on Windows. Paste the path into ParaView's File → Open. Replace `<distro>` with the name from `wsl.exe -l -q` (for example `Ubuntu-26.04`). Replace `<you>` with the Ubuntu home directory:

```text
\\wsl.localhost\<distro>\home\<you>\RDEsim\cases\annulus2D\annulus2D.foam
\\wsl$\<distro>\home\<you>\RDEsim\cases\annulus2D\annulus2D.foam
```

That file is the OpenFOAM reader stub. It is in the checkout before any run. The mesh and fields come from the case next to it. Colour these cell fields, in this order when you want to see the wave:

1. **p** — pressure
2. **T** — temperature
3. **rho** — density
4. **U** — velocity
5. **lambda** — progress variable (`0` unburned, `1` burned)

The initial directory `0/` is already in the checkout. A later time folder, the one you step to after the start, exists only after **RDEsim: Run and visualize** or **RDEsim: Run annulus2D** (`Allrun`) in WSL. The open task will not start ParaView until that folder is there. **e** (sensible internal energy) is there when the solver has written it. There is no extra detonation-indicator field. The wave, if the short window produces one, is the jump in **p**, **T**, **rho**, and **lambda**.

`VTK/geometryBlocks.vtk` is one hex per block (VTK cell type 12). It is the block outline, not the mesh `blockMesh` builds. Use `annulus2D.foam` for the fields.

`paraFoam` inside WSL is a secondary option, only if you want the reader that ships with OpenFOAM. From `cases/annulus2D`, after `. /opt/openfoam14/etc/bashrc`:

```sh
paraFoam
```

`pvbatch paraview/view_rde.py cases/annulus2D` is the same kind of optional WSL path. It was not run here. There is no `.pvsm` state file, because those files embed absolute paths. The Windows `\\wsl.localhost\` path above is the manual fallback. This environment cannot start a real Windows ParaView; the launch text was checked with stub `powershell.exe` and `tasklist.exe` only.

## Cell structure and mesh resolution

`output.cellStructure` in `rdeConfig` sets the write control for a `peakPressure` function object. `cases/annulus2D` and `cases/annulus2D-cells` enable it; `cases/annulus2D-long` does not. Each cell stores `pMax`, the running maximum of pressure, and `tAtPMax`, the time that maximum occurred. Change `writeControl` and `writeInterval` there if you want those fields on a different schedule from the flow fields.

**RDEsim: Cell structure view** prints the mesh-resolution report and runs `pvbatch` on `paraview/cell_structure.py`. The script writes a grey image of `pMax` and a wireframe of the mesh to `postProcessing/cellStructure/`. If Windows `paraview.exe` is installed, the task uses `pvbatch.exe` beside it. Otherwise it uses `pvbatch` on `PATH`. `./bin/rdesim-wsl cells --case PATH` does the same for another case, for example `cases/annulus2D-cells` (the **fine** task), and takes `--model` and `--runs`. The mesh-resolution report alone is `python3 -m rdesim resolution --config PATH/rdeConfig`.

The length in `constant/meshResolution` is the inlet sound speed divided by the one-step rate. It is not a ZND induction length and not a cell width from a detonation database. At the 2500 K kernel that length is 0.22 mm. The median axial edge on the example mesh is 2.6 mm, which is 0.086 cells across that length. The report warns that the grid is too coarse for detonation cells. A reactive Euler model with one progress variable tends to make weak or irregular cells even when the mesh is fine.

The 5e-6 s smoke run wrote `pMax` and `tAtPMax` in each time folder. `pMax` runs from about 1.01e5 Pa to 6.52e5 Pa and stays high only inside the hot kernel, about `y` ≤ 0.033 m. Around the circumference it is two broad humps, not a repeated cell pattern. The grey image shows that same patch. It is a pressure history. It is not a measured cell size. Those images were rendered here with Linux `/usr/bin/pvbatch`. Windows `pvbatch.exe` was not run.

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

## Temperature positivity

`shockFluid::fluxPredictor` evaluates `sqrt(gamma R T)`. A non-positive cell temperature is a floating-point exception, not a failed convergence message. On `annulus2D-long` that happened just after `2.4e-5` s. The failing cell sat on the injector face: pressure had fallen to about 1 kPa, speed was about 1100 m/s, and the one-step source in that cell was zero. The momentum update had put more kinetic energy in the cell than `e + 1/2|U|^2` could support at positive temperature. Tightening the Courant number or switching the limiter from vanAlbada to vanLeer did not move that time.

`positivity` in `rdeConfig` is the switch. With `enabled yes`, a cell that would fall below `TMin` has `|U|` reduced and sensible energy raised by the same amount, so `e + 1/2|U|^2` is unchanged. Energy is added only when the cell cannot reach `TMin` even at rest. The log line `Positivity guard` and `postProcessing/positivity/activations.dat` count those cell-steps. Omit the dictionary and the Kurganov update is unchanged. `cases/annulus2D` and `cases/annulus2D-cells` turn it on. `cases/annulus2D-long` does not, so it still ends at `2e-5` s. `chemistry.TFloor` is only the argument of the Arrhenius exponential. It is not this floor.

The example uses `TMin 50`. That is below the inlet temperature and below the 244 K dip seen on the coarse mesh before the crash. It is a numerical floor, not a measured gas temperature.

With that switch, `annulus2D-long` was continued from the healthy `2e-5` s fields to `1e-4` s. `foamRun` took 3.75 s and exited 0. The maximum Courant number stayed at about 0.20. No field went NaN, and density stayed positive. The guard scaled `|U|` in 78 cell-steps, two cells at a time, from `2.86e-5` s to `3.01e-5` s, and added energy in none of them. At `4e-5` s the hot products reached about `y = 0.09` m with fresh gas still ahead. By `1e-4` s those products filled the 0.12 m channel. Final temperature was 247 K to 4590 K. This is not a Chapman–Jouguet speed.

## Validation: H2-air detonation tube

`cases/detonationTubeH2` checks the solver and the one-step model against numbers computed from detailed chemistry for stoichiometric H2-air at 101325 Pa and 300 K, and calibrates the one-step fuel to them. It is a 1D tube, closed at one end, with a burnt-gas kernel that starts the wave.

```sh
./bin/rdesim-wsl validate         # about 10 minutes: build if needed, run calibrated-N10, print the table
./bin/rdesim-wsl validate-full    # about 65 minutes: also original fuel, 5 and 20 cells per half-release length
python3 -m rdesim.validation analytic   # theory only, a few seconds, no OpenFOAM
python3 -m rdesim.validation compare    # print the table from saved results
```

The Cursor tasks are **RDEsim: Run validation** and **RDEsim: Run validation (full, resolution study)**. `validate` exits with status 1 when a gated row fails. Everything adjustable is in `cases/detonationTubeH2/validationConfig` (tube, kernel, numerics, runs, tolerance bands, calibration target). The run directories under `cases/detonationTubeH2/run/` are ignored by git. Runtimes are single-core on a 4-core Linux VM. `calibrated-N5`, `calibrated-N20` and `original-N10` were run at the same time as each other, so their times include some contention. Two more separate runs of `validate`, the second after merging the later `main` (new solver build), gave the same numbers to all printed digits in 579 s and 473 s.

### Reference data

Written by `tools/validation_reference.py` to `cases/detonationTubeH2/referenceData`. Cantera 3.2, NASA CEA 3.3.4 and scipy are only needed to regenerate that file, not to run the validation.

| Quantity | Value | Source |
| --- | --- | --- |
| CJ speed | 1969.0 m/s | Cantera, `gri30.yaml`, equilibrium CJ (minimum-speed tangency). NASA CEA 3.3.4 independently gives 1968.4 m/s. `h2o2.yaml` gives 1976.3 m/s |
| CJ pressure, temperature | 15.49 atm, 2945 K | same Cantera calculation. CEA: 15.48 atm, 2942 K |
| von Neumann state | 27.46 atm, 1532 K | frozen shock at D_CJ, Cantera |
| ZND induction length | 0.2440 mm (max dT/dx), 0.2450 mm (max thermicity), `h2o2.yaml` 0.2318 mm | ZND ODE with the same mechanism, Cantera plus scipy. The induction time agrees with a constant-pressure reactor behind the von Neumann state: 0.6514 us against 0.6485 us |
| Length to 50% and 90% heat release | 0.2734 mm, 0.8130 mm | same ZND profile, heat release from the formation-enthalpy sum |
| Reduced activation energy Ea/(R T_vN) | 10.05 from induction time, 9.35 from the half-release time | shock-speed sensitivity of the ZND time, speeds D_CJ times 0.98 and 1.02 |
| Detonation cell width | 8.0 to 15.1 mm | Shepherd detonation database (Kaneshige and Shepherd, GALCIT FM97-8), five H2-air entries at 100 to 101.3 kPa, 293 to 300 K, phi 1.0 to 1.02, listed with their papers in `referenceData` |

Not verified, and not invented: a measured ZND induction length and a measured CJ pressure for this mixture. The induction length is a model result and depends on the mechanism (`gri30.yaml` and `h2o2.yaml` differ by 5 percent). The cell-width entries are measurements at slightly different conditions, so they are given as a range. The cell width is not computed here, because a 1D tube has no cells.

### What was calibrated

`fuels/hydrogen-air` is not edited. `python3 -m rdesim.validation calibrate --write` writes `fuels/hydrogen-air-validation/fuelDict` from it. Only three entries change, and the file records the original and the new value of each in a `calibration` block at its end.

- `heatRelease.LHV` is set so that the CJ speed of the one-gamma gas equals 1969.0 m/s. This is closed form.
- `chemistry.Ta` is set so that the reduced activation energy of the model's own ZND profile equals the target in `validationConfig` (5, see below).
- `chemistry.A` is set so that the model's length to half heat release equals the reference 0.2734 mm.

The `twoStep` block of the original file is copied unchanged into the calibrated file. The calibration is for the `oneStepProgress` model only, and `twoStep` is not calibrated or validated here. Nothing is hardcoded in C++ or Python. The heat release, activation temperature and pre-exponential factor are read from the fuel file by both the solver and the theory module `rdesim/detonation.py`.

Why the half-release length and not the induction length: a first-order one-step model with a reduced activation energy below about 6.5 has no maximum of dT/dx behind the shock, so its induction length is not defined. The half-release length is defined for every case.

Why the reduced activation energy is 5 and not the 9 to 10 of the real mixture: with this model and this solver, the 1D wave is not stable at the real value. At 10 cells per half-release length and 600 half-lengths of tube, one initiation kernel, a hand-run scan gave: Ea/(R T_vN) of 4 and 5 give a wave that travels at the CJ speed on average, 6 and 7 gallop with large pressure swings, 8 and 9.35 quench or fail. This is the known behaviour of one-step Arrhenius detonations. So the activation energy is a stability compromise and is not matched. To repeat the scan, set `calibration.reducedActivationEnergy` in `validationConfig`, run `calibrate --write` and `run`, then restore the fuel with `git checkout fuels/hydrogen-air-validation`.

### Results

Steady ZND/CJ theory of the one-step model, printed by `analytic` (no CFD). Error is against the Cantera `gri30` reference.

| Quantity | Reference | Before (`hydrogen-air`) | After (`hydrogen-air-validation`) |
| --- | --- | --- | --- |
| CJ speed (m/s) | 1969.0 | 2626.6 (+33.4%) | 1969.0 (+0.0%) |
| CJ pressure (atm) | 15.49 | 24.52 (+58.3%) | 13.96 (-9.9%) |
| CJ temperature (K) | 2945 | 4365 (+48.2%) | 2518 (-14.5%) |
| von Neumann pressure (atm) | 27.46 | 48.04 (+74.9%) | 26.92 (-2.0%) |
| Length to 50% heat release (mm) | 0.2734 | 0.0966 (-64.7%) | 0.2734 (-0.0%) |
| Induction length, max dT/dx (mm) | 0.2440 | 0.0307 (-87.4%) | 0.2123 (-13.0%, not a target) |
| Length to 90% heat release (mm) | 0.8130 | 0.2740 (-66.3%) | 0.5805 (-28.6%) |
| Ea/(R T_vN), half-release | 9.35 | 2.72 (-70.9%) | 5.00 (-46.5%) |

CFD, `calibrated-N10`, the gated run: 15000 cells of 27.3 um, 410 mm tube, 10 cells per model half-release length, 154 us, 599 s. The wave speed is a least-squares fit of the shock position over the 60 to 95 percent part of the tube. The pressure and temperature are read at the sonic plane behind the shock (relative Mach number 1), the CJ point.

| Quantity | Reference | CFD | Error | Band |
| --- | --- | --- | --- | --- |
| Wave speed (m/s) | 1969.0 | 1971.1 | +0.1% | 2% |
| Pressure at sonic plane (atm) | 15.49 | 14.05 | -9.3% | 15% |
| Temperature at sonic plane (K) | 2945 | 2560 | -13.1% | 20% |
| Length to 50% heat release (mm) | 0.2734 | 0.3299 | +20.7% | 25% |

The bands were chosen before the runs. They are a regression check, not a statement about physics. The pressure and temperature errors are the one-gamma limit shown in the theory table above (-9.9% and -14.5%). They cannot be removed by calibrating `A`, `Ta` and `LHV`, because the CJ speed and the CJ pressure cannot both be matched by one gamma.

Resolution and fuel, `validate-full`, same tube and kernel:

| Run | Cells | Wave speed (m/s) | p at sonic plane (atm) | T (K) | Half-release length (mm) | Runtime |
| --- | --- | --- | --- | --- | --- | --- |
| `calibrated-N5` | 7500 | 1981.5 (+0.6%) | 14.15 (-8.6%) | 2543 (-13.7%) | 0.367 (+34%, FAIL) | 168 s |
| `calibrated-N10` | 15000 | 1971.1 (+0.1%) | 14.05 (-9.3%) | 2560 (-13.1%) | 0.330 (+21%) | 599 s |
| `calibrated-N20` | 30000 | 1987.1 (+0.9%) | 14.46 (-6.6%) | 2526 (-14.2%) | 0.267 (-2.4%) | 2420 s |

| `original-N10`, original fuel, 9.7 um cells | 15000 | 2640.3 (+34.1%) | 24.76 (+59.9%) | 4384 (+48.9%) | 0.112 (-59%) | 650 s |

The original fuel is not a failure of the solver. The same table also prints a check against the model's own CJ theory: the original-fuel CFD wave is within 0.5% in speed, 1.0% in pressure and 0.4% in temperature of that theory, and the calibrated runs are within 0.6%, 1.4% and 1.6%. The solver reproduces the CJ state of the model it is given. The original numbers are simply far from H2-air.

### What still does not match

- Before calibration the original hydrogen entry gave a CJ speed 33% too high, a CJ pressure 58% too high, and a reaction zone 65% too short.
- After calibration the wave speed matches by construction. The pressure (-9% to -10%) and temperature (-13% to -15%) at the CJ plane do not, and cannot with one gamma of 1.4. The real burnt gas has a lower effective gamma.
- The reduced activation energy is 5 against 9.35. The wave in the CFD is not steady at this value. The shock speed between written frames ranges from about 1750 to 2400 m/s and the half-release length from 0.12 to 0.71 mm (median 0.33 mm, target 0.27 mm). A real H2-air detonation is also unstable, but this is not the same instability and nothing here says it is.
- The reaction-zone tail is too short. 90 percent of the heat is released after 0.58 mm in the model, against 0.81 mm in the reference. A one-step first-order model has one rate for the whole zone, while the real one has a short induction zone and a longer recombination zone.
- The model's heat release is 54 percent of the real lower heating value of the fuel (`LHV` 6.49e7 against 1.2e8 J/kg). It is an effective number that compensates for the single gamma, not a measured heat of reaction.
- The result is not grid-converged. The median half-release length is +34%, +21% and -2.4% off at 5, 10 and 20 cells per half-release length, and the spread of the half-release length between frames grows with resolution (0.21 to 0.52 mm at 5, 0.12 to 0.71 mm at 10, 0.07 to 0.89 mm at 20). The finer mesh resolves a larger oscillation of the wave, so the 20-cell agreement of the median is not evidence of convergence. The average wave speed stays within 0.1% to 0.9% of the reference at all three resolutions, and the speed does not trend with resolution.
- Cell width, the transverse structure, and the induction length itself (a target only through the half-release length) are not computed. The calibration is valid only for phi = 1, 300 K and 1 atm. The calibrated fuel has not been tested at another equivalence ratio, pressure or temperature, and nothing about the annulus cases is validated by it.

A match to the reference at the calibration point is a fit. It does not show that the model predicts H2-air detonations at other conditions.

## What is not physically validated

- The default and `twoStep` models are one-step-style, single-gamma, perfect-gas kinetics, not a detailed mechanism and not a variable-gamma product gas. `detailed` uses the Burke 2012 H2/O2 mechanism with a perfect-gas equation of state per species, Euler only; it has been checked against Cantera in 0D and against the equilibrium CJ speed in 1D, not against measured detonation data or cell widths.
- Hydrogen and ethylene entries are starting numbers, not a calibrated detonation chemistry.
- `cjcheck` compares a 1D wave with the model's own Chapman–Jouguet speed only. The 1D tube in `cases/detonationTubeH2` compares a calibrated `hydrogen-air-validation` fuel with Chapman–Jouguet speed, pressure, temperature and half-release length computed with Cantera for stoichiometric H2-air at 1 atm and 300 K, see [Validation](#validation-h2-air-detonation-tube). That is a calibration at one point, not a prediction. The annulus cases, cell size, wave count, thrust, the `twoStep` model, other fuels and other conditions are not compared with measurements. The original `hydrogen-air` one-step numbers give a CJ speed 33% too high.
- `pMax` is a pressure history. On the example mesh it is the kernel footprint, not a measured cell size. The chemical length in `constant/meshResolution` is inlet sound speed over the one-step rate.
- The `twoStep` constants are a fit to ignition delay and release duration, not a calibrated structure model. The one-step and two-step models share one gamma and one `q`, so both have a CJ speed (2627 m/s for the example) far above the real H2-air value.
- The 2D unwrapped mesh does not resolve radial steps or the nozzle wall.
- The 3D mesh does not resolve discrete injector slots.
- Acoustics are lumped. Cavities are not meshed. The Helmholtz ODE is explicit Euler and is not a proven stable integrator for every `deltaT`.
- The example end time is a smoke window.
- The 2D unwrapped `annulus2D` case was meshed with `blockMesh` and advanced with `foamRun` to `5e-6` s on OpenFOAM 14. 3D arc edges and the laminar Fourier dictionary were not executed.
- A CFL warning is about the explicit flux, not about whether a detonation will form.
- `positivity.TMin` stops a cell temperature from crossing that floor by scaling velocity. It is not a physical wall temperature, and a run that reports activations has been altered in those cells.

## License

GPL-3.0-or-later, because `rdeFluid` links OpenFOAM. See `LICENSE`.
