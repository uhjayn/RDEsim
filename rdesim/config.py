"""Load ``rdeConfig``, merge the selected fuel, and check required entries."""

from __future__ import annotations

from pathlib import Path

from .foamdict import FoamDict, FoamList, parse_file


class ConfigError(Exception):
    def __init__(self, messages: str | list[str]):
        if isinstance(messages, str):
            messages = [messages]
        self.messages = messages
        super().__init__("\n".join(messages))


def load_config(path: str | Path) -> dict:
    path = Path(path).resolve()
    raw = parse_file(path)
    if "FoamFile" in raw:
        raw = FoamDict((key, value) for key, value in raw.items() if key != "FoamFile")
    errors: list[str] = []
    _require_mapping(raw, "rdeConfig", errors)
    for section in (
        "constants",
        "domain",
        "mesh",
        "geometry",
        "operating",
        "fuel",
        "numerics",
        "hyperbolicity",
        "boundaries",
        "ignition",
        "optimization",
        "output",
    ):
        if section not in raw:
            errors.append(f"Missing required dictionary '{section}'")
        elif not isinstance(raw[section], dict):
            errors.append(f"'{section}' must be a dictionary")
    if errors:
        raise ConfigError(errors)

    fuel_raw, fuel_path = _load_fuel(path, raw["fuel"], errors)
    overrides = raw.get("fuelOverrides", FoamDict())
    if "fuelOverrides" in raw and not isinstance(overrides, dict):
        errors.append("'fuelOverrides' must be a dictionary")
        overrides = FoamDict()
    fuel = _deep_merge(fuel_raw, overrides)
    _validate_sections(raw, fuel, errors)
    if errors:
        raise ConfigError(errors)
    return {
        "path": path,
        "raw": raw,
        "fuel": fuel,
        "fuelPath": fuel_path,
    }


REACTION_MODELS = ("oneStepProgress",)


def reaction_model(config: dict) -> str:
    """Selected reaction model. Absent means the one-step default."""
    reaction = config["raw"].get("reaction")
    if isinstance(reaction, dict) and "model" in reaction:
        return str(reaction["model"])
    return "oneStepProgress"


def set_path(root: dict, path: str, value) -> None:
    parent, last = _walk(root, path, create=False)
    parent[last] = value


def get_path(root: dict, path: str):
    parent, last = _walk(root, path, create=False)
    return parent[last]


def _walk(root: dict, path: str, create: bool):
    tokens = _tokenize_path(path)
    if not tokens:
        raise ConfigError(f"Empty configuration path '{path}'")
    current = root
    for token in tokens[:-1]:
        current = _step(current, token, path)
    last = tokens[-1]
    if isinstance(last, tuple):
        current = _step(current, (last[0], None), path)
        index = last[1]
        if not isinstance(current, list) or index >= len(current):
            raise ConfigError(f"Index out of range in '{path}'")
        # The caller assigns into the list element, which must be a dict
        # only when the path ends at the index. Return the list and index.
        return current, index
    if not isinstance(current, dict) or last not in current:
        raise ConfigError(f"Missing entry '{path}'")
    return current, last


def _step(current, token, path):
    if isinstance(token, tuple):
        key, index = token
        if not isinstance(current, dict) or key not in current:
            raise ConfigError(f"Missing entry '{key}' while resolving '{path}'")
        sequence = current[key]
        if not isinstance(sequence, list):
            raise ConfigError(f"'{key}' is not a list in '{path}'")
        if index is None:
            return sequence
        if index < 0 or index >= len(sequence):
            raise ConfigError(f"Index {index} out of range in '{path}'")
        return sequence[index]
    if not isinstance(current, dict) or token not in current:
        raise ConfigError(f"Missing entry '{token}' while resolving '{path}'")
    return current[token]


def _tokenize_path(path: str) -> list:
    tokens = []
    for part in path.split("."):
        if not part:
            raise ConfigError(f"Invalid configuration path '{path}'")
        if "[" in part:
            name, rest = part.split("[", 1)
            if not rest.endswith("]"):
                raise ConfigError(f"Invalid configuration path '{path}'")
            tokens.append((name, int(rest[:-1])))
        else:
            tokens.append(part)
    return tokens


def _load_fuel(config_path: Path, fuel_dict: dict, errors: list[str]):
    if "name" not in fuel_dict:
        errors.append("Missing required entry 'fuel.name'")
    if "library" not in fuel_dict:
        errors.append("Missing required entry 'fuel.library'")
    if errors:
        return FoamDict(), Path()
    name = str(fuel_dict["name"])
    library = _resolve_library(config_path, str(fuel_dict["library"]))
    fuel_path = library / name / "fuelDict"
    if not fuel_path.is_file():
        errors.append(
            f"Fuel '{name}' was not found at {fuel_path}. "
            "Add a directory under the fuel library containing fuelDict."
        )
        return FoamDict(), fuel_path
    data = parse_file(fuel_path)
    if "FoamFile" in data:
        data = FoamDict((key, value) for key, value in data.items() if key != "FoamFile")
    return data, fuel_path


def _resolve_library(config_path: Path, library: str) -> Path:
    candidate = Path(library)
    if candidate.is_absolute() and candidate.is_dir():
        return candidate
    relative = config_path.parent / library
    if relative.is_dir():
        return relative.resolve()
    for parent in config_path.parents:
        found = parent / library
        if found.is_dir():
            return found.resolve()
    raise ConfigError(
        f"Fuel library '{library}' was not found relative to {config_path} "
        "or any parent directory."
    )


def _deep_merge(base, override):
    if not isinstance(override, dict):
        return override
    merged = FoamDict(base) if isinstance(base, dict) else FoamDict()
    for key, value in override.items():
        if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def _require_mapping(value, label: str, errors: list[str]) -> bool:
    if not isinstance(value, dict):
        errors.append(f"'{label}' must be a dictionary")
        return False
    return True


def _need(mapping, key: str, label: str, errors: list[str]):
    if not isinstance(mapping, dict) or key not in mapping:
        errors.append(f"Missing required entry '{label}'")
        return None
    return mapping[key]


def _need_number(mapping, key: str, label: str, errors: list[str]):
    value = _need(mapping, key, label, errors)
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        errors.append(f"'{label}' must be a number")
        return None
    return float(value)


def _need_word(mapping, key: str, label: str, errors: list[str], choices=None):
    value = _need(mapping, key, label, errors)
    if value is None:
        return None
    rendered = _as_words(value)
    if rendered is None:
        errors.append(f"'{label}' must be a word or a sequence of words")
        return None
    if choices is not None and rendered not in choices:
        errors.append(
            f"'{label}' is '{rendered}'. Implemented choices: {', '.join(choices)}"
        )
    return rendered


def _as_words(value):
    if isinstance(value, str):
        return value
    if isinstance(value, list) and value and all(isinstance(item, str) for item in value):
        return " ".join(value)
    return None


def _need_bool(mapping, key: str, label: str, errors: list[str]):
    value = _need(mapping, key, label, errors)
    if value is None:
        return None
    if isinstance(value, str) and value.lower() in ("yes", "true", "on"):
        return True
    if isinstance(value, str) and value.lower() in ("no", "false", "off"):
        return False
    errors.append(f"'{label}' must be yes or no")
    return None


def _validate_sections(raw: dict, fuel: dict, errors: list[str]) -> None:
    constants = raw["constants"]
    _need_number(constants, "RUniversal", "constants.RUniversal", errors)

    domain = raw["domain"]
    dimension = _need_number(domain, "dimension", "domain.dimension", errors)
    if dimension is not None and int(dimension) not in (2, 3):
        errors.append("'domain.dimension' must be 2 or 3")
    topology = _need_word(
        domain,
        "topology",
        "domain.topology",
        errors,
        ("unwrapped", "meridional"),
    )
    for key in (
        "RInner",
        "ROuter",
        "axialLength",
        "sectorAngle",
        "emptyThickness",
        "axialReference",
        "theta0",
    ):
        _need_number(domain, key, f"domain.{key}", errors)
    _need_word(
        domain,
        "unwrapRadius",
        "domain.unwrapRadius",
        errors,
        ("mean", "inner", "outer"),
    )
    if dimension == 2 and topology == "unwrapped":
        pass
    if dimension == 3 and topology == "meridional":
        errors.append("'domain.topology' meridional is a 2D mesh. Set dimension 2.")

    mesh = raw["mesh"]
    for key in (
        "nCircumferential",
        "nAxial",
        "nRadial",
        "nCircumferentialBlocks",
        "minCellsPerBlock",
    ):
        value = _need_number(mesh, key, f"mesh.{key}", errors)
        if value is not None and int(value) != value:
            errors.append(f"'mesh.{key}' must be an integer")
        if value is not None and value < 1 and key != "nCircumferentialBlocks":
            errors.append(f"'mesh.{key}' must be at least 1")
    for key in (
        "axialGrading",
        "circumferentialGrading",
        "radialGrading",
        "stationTolerance",
    ):
        _need_number(mesh, key, f"mesh.{key}", errors)

    geometry = raw["geometry"]
    _need_number(geometry, "innerWallAngle", "geometry.innerWallAngle", errors)
    _need_number(geometry, "outerWallAngle", "geometry.outerWallAngle", errors)
    _need_word(
        geometry,
        "widthAnchor",
        "geometry.widthAnchor",
        errors,
        ("inner", "outer", "centre"),
    )
    if "widthProfile" not in geometry or not isinstance(geometry["widthProfile"], list):
        errors.append("'geometry.widthProfile' must be a list")
    else:
        for index, entry in enumerate(geometry["widthProfile"]):
            if not isinstance(entry, dict):
                errors.append(f"'geometry.widthProfile[{index}]' must be a dictionary")
                continue
            _need_number(entry, "location", f"geometry.widthProfile[{index}].location", errors)
            _need_number(entry, "gap", f"geometry.widthProfile[{index}].gap", errors)
    if "steps" not in geometry or not isinstance(geometry["steps"], list):
        errors.append("'geometry.steps' must be a list")
    else:
        for index, entry in enumerate(geometry["steps"]):
            if not isinstance(entry, dict):
                errors.append(f"'geometry.steps[{index}]' must be a dictionary")
                continue
            _need_word(
                entry,
                "side",
                f"geometry.steps[{index}].side",
                errors,
                ("inner", "outer"),
            )
            _need_number(entry, "location", f"geometry.steps[{index}].location", errors)
            _need_number(entry, "height", f"geometry.steps[{index}].height", errors)
            _need_number(entry, "angle", f"geometry.steps[{index}].angle", errors)
    nozzle = _need(geometry, "nozzle", "geometry.nozzle", errors)
    if isinstance(nozzle, dict):
        _need_bool(nozzle, "enabled", "geometry.nozzle.enabled", errors)
        for key in ("start", "length", "throatGap", "exitGap", "throatFraction"):
            _need_number(nozzle, key, f"geometry.nozzle.{key}", errors)
        _need_word(
            nozzle,
            "shape",
            "geometry.nozzle.shape",
            errors,
            ("conical", "piecewise"),
        )
        _need_word(
            nozzle,
            "anchor",
            "geometry.nozzle.anchor",
            errors,
            ("inner", "outer", "centre"),
        )
        if "profile" not in nozzle or not isinstance(nozzle["profile"], list):
            errors.append("'geometry.nozzle.profile' must be a list")
    injector = _need(geometry, "injector", "geometry.injector", errors)
    if isinstance(injector, dict):
        _need_number(injector, "nSlots", "geometry.injector.nSlots", errors)
        for key in ("slotWidth", "recessDepth", "injectionAngle", "phase"):
            _need_number(injector, key, f"geometry.injector.{key}", errors)

    operating = raw["operating"]
    for key in (
        "pInlet",
        "TInlet",
        "pOutlet",
        "TOutlet",
        "equivalenceRatio",
        "massFlowRate",
        "plenumPressure",
        "plenumTemperature",
    ):
        _need_number(operating, key, f"operating.{key}", errors)
    _need_word(
        operating,
        "inletModel",
        "operating.inletModel",
        errors,
        ("fixedVelocity", "massFlow", "isentropicPlenum"),
    )
    velocity = _need(operating, "inletVelocity", "operating.inletVelocity", errors)
    if isinstance(velocity, dict):
        for key in ("circumferential", "axial", "radial"):
            _need_number(velocity, key, f"operating.inletVelocity.{key}", errors)

    _validate_fuel(fuel, errors)

    numerics = raw["numerics"]
    _need_word(
        numerics,
        "flowModel",
        "numerics.flowModel",
        errors,
        ("euler", "navierStokes"),
    )
    _need_word(numerics, "fluxScheme", "numerics.fluxScheme", errors)
    for key in (
        "ddtScheme",
        "gradScheme",
        "laplacianScheme",
        "reconstructScalar",
        "reconstructVector",
        "fluxInterpolation",
        "writeFormat",
    ):
        _need_word(numerics, key, f"numerics.{key}", errors)
    for key in (
        "CFL",
        "deltaT",
        "maxDeltaT",
        "endTime",
        "deltaTFactor",
        "startTime",
        "writePrecision",
        "timePrecision",
    ):
        _need_number(numerics, key, f"numerics.{key}", errors)
    _need_bool(numerics, "adjustTimeStep", "numerics.adjustTimeStep", errors)
    _need_bool(numerics, "runTimeModifiable", "numerics.runTimeModifiable", errors)
    if "startFrom" in numerics:
        _need_word(
            numerics,
            "startFrom",
            "numerics.startFrom",
            errors,
            ("startTime", "latestTime"),
        )
    _need_word(
        numerics,
        "sourceIntegration",
        "numerics.sourceIntegration",
        errors,
        ("pointImplicit", "explicit"),
    )
    solution = _need(numerics, "fvSolution", "numerics.fvSolution", errors)
    if isinstance(solution, dict):
        _need(solution, "solvers", "numerics.fvSolution.solvers", errors)
        _need(solution, "PIMPLE", "numerics.fvSolution.PIMPLE", errors)

    hyper = raw["hyperbolicity"]
    _need_bool(
        hyper,
        "requireStrictlyHyperbolic",
        "hyperbolicity.requireStrictlyHyperbolic",
        errors,
    )
    for key in ("gammaMin", "rhoMin", "pMin", "TMin", "c2Min", "maxCoStability"):
        _need_number(hyper, key, f"hyperbolicity.{key}", errors)
    for key in (
        "hyperbolicFluxSchemes",
        "hyperbolicEquationsOfState",
        "localSourceIntegrations",
    ):
        value = _need(hyper, key, f"hyperbolicity.{key}", errors)
        if value is not None and not isinstance(value, list):
            errors.append(f"'hyperbolicity.{key}' must be a list")
    _need_bool(
        hyper,
        "sourceDependsOnGradients",
        "hyperbolicity.sourceDependsOnGradients",
        errors,
    )
    _need_word(hyper, "boundaryCoupling", "hyperbolicity.boundaryCoupling", errors)
    _validate_positivity(raw, errors)

    _validate_reaction(raw, errors)
    _validate_boundaries(raw["boundaries"], errors)
    _validate_ignition(raw["ignition"], errors)
    _validate_output(raw["output"], errors)
    _validate_optimization(raw["optimization"], errors)


def _validate_positivity(raw: dict, errors: list[str]) -> None:
    block = raw.get("positivity")
    if block is None:
        return
    if not isinstance(block, dict):
        errors.append("'positivity' must be a dictionary")
        return
    enabled = _need_bool(block, "enabled", "positivity.enabled", errors)
    if enabled:
        minimum = _need_number(block, "TMin", "positivity.TMin", errors)
        if minimum is not None and float(minimum) <= 0:
            errors.append("'positivity.TMin' must be positive")


def _validate_reaction(raw: dict, errors: list[str]) -> None:
    reaction = raw.get("reaction")
    if reaction is None:
        return
    if not isinstance(reaction, dict):
        errors.append("'reaction' must be a dictionary")
        return
    if "model" in reaction:
        _need_word(reaction, "model", "reaction.model", errors, REACTION_MODELS)


def _validate_fuel(fuel: dict, errors: list[str]) -> None:
    fuel_specie = _need(fuel, "fuel", "fuel.fuel", errors)
    if isinstance(fuel_specie, dict):
        _need_number(fuel_specie, "molWeight", "fuel.fuel.molWeight", errors)
    oxidizer = _need(fuel, "oxidizer", "fuel.oxidizer", errors)
    if isinstance(oxidizer, dict):
        _need_number(
            oxidizer,
            "oxygenMoleFraction",
            "fuel.oxidizer.oxygenMoleFraction",
            errors,
        )
        species = _need(oxidizer, "species", "fuel.oxidizer.species", errors)
        if isinstance(species, list):
            if not species:
                errors.append("'fuel.oxidizer.species' must list at least one species")
            for index, entry in enumerate(species):
                if not isinstance(entry, dict):
                    errors.append(f"'fuel.oxidizer.species[{index}]' must be a dictionary")
                    continue
                _need_word(entry, "name", f"fuel.oxidizer.species[{index}].name", errors)
                _need_number(
                    entry,
                    "molWeight",
                    f"fuel.oxidizer.species[{index}].molWeight",
                    errors,
                )
                _need_number(
                    entry,
                    "moleFraction",
                    f"fuel.oxidizer.species[{index}].moleFraction",
                    errors,
                )
    stoich = _need(fuel, "stoichiometry", "fuel.stoichiometry", errors)
    if isinstance(stoich, dict):
        _need_number(stoich, "nuFuel", "fuel.stoichiometry.nuFuel", errors)
        _need_number(stoich, "nuOxygen", "fuel.stoichiometry.nuOxygen", errors)
    heat = _need(fuel, "heatRelease", "fuel.heatRelease", errors)
    if isinstance(heat, dict):
        _need_number(heat, "LHV", "fuel.heatRelease.LHV", errors)
    thermo = _need(fuel, "thermo", "fuel.thermo", errors)
    if isinstance(thermo, dict):
        _need_word(thermo, "model", "fuel.thermo.model", errors, ("hConst",))
        _need_word(
            thermo,
            "equationOfState",
            "fuel.thermo.equationOfState",
            errors,
            ("perfectGas",),
        )
        _need_word(
            thermo,
            "composition",
            "fuel.thermo.composition",
            errors,
            ("frozenReactants", "frozenProducts"),
        )
        _need_number(thermo, "gamma", "fuel.thermo.gamma", errors)
        _need_number(thermo, "hf", "fuel.thermo.hf", errors)
        _need_number(thermo, "productMolWeight", "fuel.thermo.productMolWeight", errors)
        _need_number(thermo, "productGamma", "fuel.thermo.productGamma", errors)
    transport = _need(fuel, "transport", "fuel.transport", errors)
    if isinstance(transport, dict):
        _need_number(transport, "mu", "fuel.transport.mu", errors)
        _need_number(transport, "Pr", "fuel.transport.Pr", errors)
    chemistry = _need(fuel, "chemistry", "fuel.chemistry", errors)
    if isinstance(chemistry, dict):
        _need_word(
            chemistry,
            "model",
            "fuel.chemistry.model",
            errors,
            ("oneStepProgress",),
        )
        for key in (
            "A",
            "beta",
            "Ta",
            "TRef",
            "TFloor",
            "ignitionTemperature",
            "lambdaMin",
            "lambdaMax",
        ):
            _need_number(chemistry, key, f"fuel.chemistry.{key}", errors)


def _validate_boundaries(boundaries: dict, errors: list[str]) -> None:
    inlet = _need(boundaries, "inlet", "boundaries.inlet", errors)
    if isinstance(inlet, dict):
        for key in ("p", "T", "U"):
            _need_word(inlet, key, f"boundaries.inlet.{key}", errors)
    outlet = _need(boundaries, "outlet", "boundaries.outlet", errors)
    if isinstance(outlet, dict):
        _need_word(
            outlet,
            "type",
            "boundaries.outlet.type",
            errors,
            ("waveTransmissive", "zeroGradient", "fixedValue"),
        )
        _need_number(outlet, "lInf", "boundaries.outlet.lInf", errors)
    walls = _need(boundaries, "walls", "boundaries.walls", errors)
    if isinstance(walls, dict):
        _need_word(
            walls,
            "velocityModel",
            "boundaries.walls.velocityModel",
            errors,
            ("slip", "acoustic"),
        )
    acoustics = _need(boundaries, "acoustics", "boundaries.acoustics", errors)
    if isinstance(acoustics, dict):
        patches = _need(acoustics, "targetPatches", "boundaries.acoustics.targetPatches", errors)
        if patches is not None and not isinstance(patches, list):
            errors.append("'boundaries.acoustics.targetPatches' must be a list")
        _need_word(
            acoustics,
            "mode",
            "boundaries.acoustics.mode",
            errors,
            ("impedance", "reflection", "microstructure"),
        )
        _need_word(
            acoustics,
            "timeDomainModel",
            "boundaries.acoustics.timeDomainModel",
            errors,
            ("slip", "resistive", "helmholtz"),
        )
        for key in (
            "frequency",
            "incidenceAngleDeg",
            "ZReal",
            "ZImag",
            "RReal",
            "RImag",
            "pRef",
            "resistiveFloor",
        ):
            _need_number(acoustics, key, f"boundaries.acoustics.{key}", errors)
        micro = _need(acoustics, "microstructure", "boundaries.acoustics.microstructure", errors)
        if isinstance(micro, dict):
            _need_word(
                micro,
                "type",
                "boundaries.acoustics.microstructure.type",
                errors,
                ("helmholtz", "quarterWave", "groove"),
            )
            for key in (
                "neckLength",
                "neckRadius",
                "cavityVolume",
                "endCorrection",
                "cellArea",
                "normalizedResistance",
                "cavityDepth",
                "grooveWidth",
                "grooveLand",
                "viscousResistance",
                "initialDisplacement",
                "initialVelocity",
            ):
                _need_number(
                    micro,
                    key,
                    f"boundaries.acoustics.microstructure.{key}",
                    errors,
                )


def _validate_ignition(ignition: dict, errors: list[str]) -> None:
    _need_word(
        ignition,
        "type",
        "ignition.type",
        errors,
        ("none", "hotKernel", "progressKernel", "hotProgressKernel"),
    )
    for key in ("boxMin", "boxMax"):
        value = _need(ignition, key, f"ignition.{key}", errors)
        if isinstance(value, list) and len(value) != 3:
            errors.append(f"'ignition.{key}' must contain three coordinates")
    for key in (
        "TKernel",
        "lambdaKernel",
        "initialPressure",
        "initialTemperature",
        "initialLambda",
        "inletLambda",
    ):
        _need_number(ignition, key, f"ignition.{key}", errors)
    velocity = _need(ignition, "initialVelocity", "ignition.initialVelocity", errors)
    if isinstance(velocity, dict):
        for key in ("circumferential", "axial", "radial"):
            _need_number(velocity, key, f"ignition.initialVelocity.{key}", errors)


def _validate_output(output: dict, errors: list[str]) -> None:
    fields = _need(output, "fields", "output.fields", errors)
    if fields is not None and not isinstance(fields, list):
        errors.append("'output.fields' must be a list")
    _need_word(output, "writeControl", "output.writeControl", errors)
    _need_number(output, "writeInterval", "output.writeInterval", errors)
    _need_number(output, "purgeWrite", "output.purgeWrite", errors)
    _need_bool(output, "writeCompression", "output.writeCompression", errors)
    if "startFrom" in output:
        _need_word(
            output,
            "startFrom",
            "output.startFrom",
            errors,
            ("startTime", "latestTime", "firstTime"),
        )
    checkpoint = output.get("checkpoint")
    if checkpoint is not None:
        if not isinstance(checkpoint, dict):
            errors.append("'output.checkpoint' must be a dictionary")
        else:
            _need_bool(checkpoint, "enabled", "output.checkpoint.enabled", errors)
            _need_word(
                checkpoint,
                "writeControl",
                "output.checkpoint.writeControl",
                errors,
                ("clockTime", "cpuTime", "runTime", "adjustableRunTime", "timeStep"),
            )
            _need_number(
                checkpoint, "writeInterval", "output.checkpoint.writeInterval", errors
            )
    objectives = _need(output, "objectives", "output.objectives", errors)
    if isinstance(objectives, list):
        for index, entry in enumerate(objectives):
            if not isinstance(entry, dict):
                errors.append(f"'output.objectives[{index}]' must be a dictionary")
                continue
            _need_word(entry, "name", f"output.objectives[{index}].name", errors)
            _need_word(entry, "field", f"output.objectives[{index}].field", errors)
            region = _need_word(
                entry,
                "region",
                f"output.objectives[{index}].region",
                errors,
                ("patch", "volume"),
            )
            _need_word(entry, "operation", f"output.objectives[{index}].operation", errors)
            if region == "patch":
                _need_word(entry, "patch", f"output.objectives[{index}].patch", errors)
    _validate_cell_structure(output, errors)


def _validate_cell_structure(output: dict, errors: list[str]) -> None:
    block = output.get("cellStructure")
    if block is None:
        return
    if not isinstance(block, dict):
        errors.append("'output.cellStructure' must be a dictionary")
        return
    if "enabled" in block:
        _need_bool(block, "enabled", "output.cellStructure.enabled", errors)
    if "writeControl" in block:
        _need_word(
            block,
            "writeControl",
            "output.cellStructure.writeControl",
            errors,
            ("timeStep", "writeTime", "adjustableRunTime", "runTime"),
        )
    if "writeInterval" in block:
        _need_number(block, "writeInterval", "output.cellStructure.writeInterval", errors)


def _validate_optimization(opt: dict, errors: list[str]) -> None:
    _need_word(
        opt,
        "method",
        "optimization.method",
        errors,
        ("grid", "coordinate"),
    )
    _need_number(opt, "maxEvaluations", "optimization.maxEvaluations", errors)
    _need_number(opt, "timeoutSeconds", "optimization.timeoutSeconds", errors)
    _need_number(opt, "coordinateLevels", "optimization.coordinateLevels", errors)
    for key in ("maxNewEvaluations",):
        if key in opt:
            value = _need_number(opt, key, f"optimization.{key}", errors)
            if value is not None and value < 0:
                errors.append(f"'optimization.{key}' must not be negative")
    if "resume" in opt:
        _need_bool(opt, "resume", "optimization.resume", errors)
    _need_word(opt, "workDirectory", "optimization.workDirectory", errors)
    _need_word(opt, "solverCommand", "optimization.solverCommand", errors)
    objective = _need(opt, "objective", "optimization.objective", errors)
    if isinstance(objective, dict):
        _need_word(objective, "name", "optimization.objective.name", errors)
        _need_word(
            objective,
            "sense",
            "optimization.objective.sense",
            errors,
            ("maximize", "minimize"),
        )
        _need_word(objective, "expression", "optimization.objective.expression", errors)
    if "constraints" not in opt or not isinstance(opt["constraints"], list):
        errors.append("'optimization.constraints' must be a list")
    else:
        for index, entry in enumerate(opt["constraints"]):
            if not isinstance(entry, dict):
                errors.append(f"'optimization.constraints[{index}]' must be a dictionary")
                continue
            _need_word(entry, "quantity", f"optimization.constraints[{index}].quantity", errors)
            if "min" not in entry and "max" not in entry:
                errors.append(
                    f"'optimization.constraints[{index}]' needs min, max, or both"
                )
    variables = _need(opt, "designVariables", "optimization.designVariables", errors)
    if isinstance(variables, list):
        if not variables:
            errors.append("'optimization.designVariables' must contain at least one variable")
        for index, entry in enumerate(variables):
            if not isinstance(entry, dict):
                errors.append(
                    f"'optimization.designVariables[{index}]' must be a dictionary"
                )
                continue
            _need_word(entry, "name", f"optimization.designVariables[{index}].name", errors)
            _need_number(entry, "min", f"optimization.designVariables[{index}].min", errors)
            _need_number(entry, "max", f"optimization.designVariables[{index}].max", errors)
            samples = _need_number(
                entry,
                "samples",
                f"optimization.designVariables[{index}].samples",
                errors,
            )
            if samples is not None and samples < 1:
                errors.append(
                    f"'optimization.designVariables[{index}].samples' must be at least 1"
                )


def as_bool(value) -> bool:
    if isinstance(value, str):
        return value.lower() in ("yes", "true", "on")
    return bool(value)


def as_list(value) -> list:
    if isinstance(value, FoamList):
        return list(value)
    if isinstance(value, list):
        return value
    return [value]
