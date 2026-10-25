"""Write an OpenFOAM 14 case from a resolved rdeConfig."""

from __future__ import annotations

from pathlib import Path

from .acoustics import evaluate_acoustics
from .config import ConfigError, as_bool
from .foamdict import format_data, format_number, render_entries, write_file
from .geometry import build_profile
from .hyperbolicity import check_hyperbolicity
from .mesh import build_mesh
from .mixture import inlet_speed, mesh_velocity, mixture_from_config, velocity_components
from .resolution import resolution_report
from .vtkwrite import write_blocks


def prepare_case(config: dict, case_dir: str | Path) -> dict:
    case_dir = Path(case_dir)
    state = mixture_from_config(config)
    profile = build_profile(config)
    mesh = build_mesh(config, profile)
    acoustics = evaluate_acoustics(config, state.density, state.sound_speed)
    report = check_hyperbolicity(config, state)
    prepared = dict(config)
    prepared["_injectorArea"] = profile.injector_area
    components = velocity_components(
        prepared,
        inlet_speed(prepared, state, profile.injector_area),
    )
    prepared["_inletVelocity"] = mesh_velocity(
        components,
        mesh.topology,
        mesh.dimension,
        False,
    )
    prepared["_boundaryNames"] = _collect_boundary_names(mesh)
    _write_tree(case_dir, prepared, state, profile, mesh, acoustics, report)
    return {
        "case": case_dir,
        "mixture": state,
        "profile": profile,
        "mesh": mesh,
        "acoustics": acoustics,
        "hyperbolicity": report,
    }


def _write_tree(case_dir, config, state, profile, mesh, acoustics, report) -> None:
    system = case_dir / "system"
    constant = case_dir / "constant"
    zero = case_dir / "0"
    for directory in (system, constant, zero, case_dir / "VTK"):
        directory.mkdir(parents=True, exist_ok=True)
    (case_dir / f"{case_dir.name}.foam").write_text("", encoding="utf-8")
    write_file(system / "blockMeshDict", "dictionary", "blockMeshDict", mesh.block_mesh_dict())
    write_file(system / "controlDict", "dictionary", "controlDict", _control(config))
    write_file(system / "fvSchemes", "dictionary", "fvSchemes", _schemes(config))
    write_file(system / "fvSolution", "dictionary", "fvSolution", _solution(config))
    write_file(
        constant / "physicalProperties",
        "dictionary",
        "physicalProperties",
        _thermo(config, state),
    )
    write_file(
        constant / "momentumTransport",
        "dictionary",
        "momentumTransport",
        _momentum(config),
    )
    if "thermophysicalTransport" in config["raw"]["numerics"]:
        write_file(
            constant / "thermophysicalTransport",
            "dictionary",
            "thermophysicalTransport",
            render_entries(config["raw"]["numerics"]["thermophysicalTransport"]) + "\n",
        )
    elif config["raw"]["numerics"]["flowModel"] == "navierStokes":
        raise ConfigError(
            "numerics.flowModel navierStokes requires numerics.thermophysicalTransport"
        )
    write_file(constant / "fvModels", "dictionary", "fvModels", "// No fvModels in this case.\n")
    write_file(
        constant / "fvConstraints",
        "dictionary",
        "fvConstraints",
        "// No fvConstraints in this case.\n",
    )
    write_file(constant / "rdeProperties", "dictionary", "rdeProperties", _properties(config, state))
    write_file(constant / "wallProfile", "dictionary", "wallProfile", _wall(profile))
    write_file(constant / "geometryReport", "dictionary", "geometryReport", _geometry_report(profile, mesh))
    write_file(
        constant / "acousticReflection",
        "dictionary",
        "acousticReflection",
        _acoustic(acoustics),
    )
    _fields(zero, config, state, mesh, acoustics)
    write_blocks(case_dir / "VTK" / "geometryBlocks.vtk", mesh.hexes)
    _scripts(case_dir)
    (case_dir / "constant" / "hyperbolicityReport").write_text(report.text() + "\n", encoding="utf-8")
    resolution = resolution_report(config)
    (case_dir / "constant" / "meshResolution").write_text(resolution.text + "\n", encoding="utf-8")


def _control(config) -> str:
    numerics = config["raw"]["numerics"]
    output = config["raw"]["output"]
    adjust = "yes" if as_bool(numerics["adjustTimeStep"]) else "no"
    compression = "on" if as_bool(output["writeCompression"]) else "off"
    modifiable = "yes" if as_bool(numerics["runTimeModifiable"]) else "no"
    lines = [
        "solver          rdeFluid;",
        "",
        "startFrom       startTime;",
        f"startTime       {format_number(numerics['startTime'])};",
        "stopAt          endTime;",
        f"endTime         {format_number(numerics['endTime'])};",
        f"deltaT          {format_number(numerics['deltaT'])};",
        f"adjustTimeStep  {adjust};",
        f"maxCo           {format_number(numerics['CFL'])};",
        f"maxDeltaT       {format_number(numerics['maxDeltaT'])};",
        f"deltaTFactor    {format_number(numerics['deltaTFactor'])};",
        "",
        f"writeControl    {output['writeControl']};",
        f"writeInterval   {format_number(output['writeInterval'])};",
        f"purgeWrite      {format_number(output['purgeWrite'])};",
        f"writeFormat     {numerics['writeFormat']};",
        f"writePrecision  {format_number(numerics['writePrecision'])};",
        f"writeCompression {compression};",
        "timeFormat      general;",
        f"timePrecision   {format_number(numerics['timePrecision'])};",
        f"runTimeModifiable {modifiable};",
        "",
        "functions",
        "{",
    ]
    for objective in config["raw"]["output"]["objectives"]:
        lines.extend(_function_object(objective, output))
    peak = _peak_pressure_function(output)
    if peak:
        lines.append(peak)
    lines.append("}")
    return "\n".join(lines) + "\n"


def _peak_pressure_function(output) -> str:
    block = output.get("cellStructure")
    if not isinstance(block, dict):
        block = {}
    enabled = as_bool(block["enabled"]) if "enabled" in block else True
    if not enabled:
        return ""
    control = block.get("writeControl", output["writeControl"])
    interval = format_number(block.get("writeInterval", output["writeInterval"]))
    return f"""    peakPressure
    {{
        type            peakPressure;
        libs            ("librdeFluidSolver.so");
        enabled         yes;
        executeAtStart  no;
        executeControl  timeStep;
        executeInterval 1;
        writeControl    {control};
        writeInterval   {interval};
    }}"""


def _function_object(objective, output) -> list[str]:
    name = objective["name"]
    field = objective["field"]
    operation = objective["operation"]
    write_control = output["writeControl"]
    interval = format_number(output["writeInterval"])
    if objective["region"] == "patch":
        body = f"""    {name}
    {{
        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");
        writeControl    {write_control};
        writeInterval   {interval};
        writeFields     false;
        log             false;
        operation       {operation};
        patch           {objective["patch"]};
        fields          ({field});
    }}"""
    else:
        body = f"""    {name}
    {{
        type            volFieldValue;
        libs            ("libfieldFunctionObjects.so");
        writeControl    {write_control};
        writeInterval   {interval};
        writeFields     false;
        log             false;
        operation       {operation};
        cellZone        all;
        fields          ({field});
    }}"""
    return [body]


def _scheme(value) -> str:
    return format_data(value)


def _schemes(config) -> str:
    numerics = config["raw"]["numerics"]
    scalar = _scheme(numerics["reconstructScalar"])
    vector = _scheme(numerics["reconstructVector"])
    return f"""fluxScheme      {_scheme(numerics["fluxScheme"])};

ddtSchemes
{{
    default         {_scheme(numerics["ddtScheme"])};
}}

gradSchemes
{{
    default         {_scheme(numerics["gradScheme"])};
}}

divSchemes
{{
    default         none;
}}

laplacianSchemes
{{
    default         {_scheme(numerics["laplacianScheme"])};
}}

interpolationSchemes
{{
    default         none;
    reconstruct(rho) {scalar};
    reconstruct(U)  {vector};
    reconstruct(T)  {scalar};
    reconstruct(e)  {scalar};
    reconstruct(lambda) {scalar};
    flux(rhoU)      {_scheme(numerics["fluxInterpolation"])};
}}

snGradSchemes
{{
    default         orthogonal;
}}
"""


def _solution(config) -> str:
    return render_entries(config["raw"]["numerics"]["fvSolution"]) + "\n"


def _openfoam_mol_weight(state) -> float:
    """kg/kmol entry that makes OpenFOAM's specie::R() match this mixture.

    Foundation OpenFOAM 14 stores molWeight in kg/kmol and uses
    R = 8314.462618/molWeight. Fuel files stay in kg/mol. Multiplying by
    1000 matches that gas constant when constants.RUniversal is
    8.314462618 J/(mol·K). Dividing the OpenFOAM constant by the mixture
    specific gas constant keeps the same gamma if RUniversal differs.
    """

    openfoam_rr = 8314.462618
    return openfoam_rr / state.specific_gas


def _thermo(config, state) -> str:
    mu = state.mu
    return f"""thermoType
{{
    type            hePsiThermo;
    mixture         pureMixture;
    transport       const;
    thermo          hConst;
    equationOfState {state.equation_of_state};
    specie          specie;
    energy          sensibleInternalEnergy;
}}

mixture
{{
    specie
    {{
        molWeight   {format_number(_openfoam_mol_weight(state))};
    }}
    thermodynamics
    {{
        Cp          {format_number(state.cp)};
        hf          {format_number(state.hf)};
    }}
    transport
    {{
        mu          {format_number(mu)};
        Pr          {format_number(state.prandtl)};
    }}
}}
"""


def _momentum(config) -> str:
    numerics = config["raw"]["numerics"]
    if "momentumTransport" in numerics:
        return render_entries(numerics["momentumTransport"]) + "\n"
    return "simulationType  laminar;\n"


def _properties(config, state) -> str:
    chemistry = state.chemistry
    hyper = config["raw"]["hyperbolicity"]
    numerics = config["raw"]["numerics"]
    lines = [
        "flowModel          " + numerics["flowModel"] + ";",
        "equationOfState    " + state.equation_of_state + ";",
        "ddtScheme          " + numerics["ddtScheme"] + ";",
        "",
        "chemistry",
        "{",
        "    model                  " + chemistry["model"] + ";",
        "    A                      " + format_number(chemistry["A"]) + ";",
        "    beta                   " + format_number(chemistry["beta"]) + ";",
        "    Ta                     " + format_number(chemistry["Ta"]) + ";",
        "    TRef                   " + format_number(chemistry["TRef"]) + ";",
        "    TFloor                 " + format_number(chemistry["TFloor"]) + ";",
        "    ignitionTemperature    " + format_number(chemistry["ignitionTemperature"]) + ";",
        "    lambdaMin              " + format_number(chemistry["lambdaMin"]) + ";",
        "    lambdaMax              " + format_number(chemistry["lambdaMax"]) + ";",
        "    q                      " + format_number(state.heat_release) + ";",
        "    sourceIntegration      " + numerics["sourceIntegration"] + ";",
        "    sourceDependsOnGradients " + hyper["sourceDependsOnGradients"] + ";",
        "}",
        "",
        "hyperbolicity",
        "{",
        render_entries(hyper, indent=1),
        "}",
        "",
        _positivity_block(config),
        "objectives",
        "(",
    ]
    for objective in config["raw"]["output"]["objectives"]:
        lines.append("    {")
        lines.append(render_entries(objective, indent=2))
        lines.append("    }")
    lines.append(");")
    return "\n".join(lines) + "\n"


def _positivity_block(config) -> str:
    block = config["raw"].get("positivity")
    if not isinstance(block, dict):
        return ""
    enabled = "yes" if as_bool(block["enabled"]) else "no"
    lines = [
        "positivity",
        "{",
        "    enabled    " + enabled + ";",
    ]
    if "TMin" in block:
        lines.append("    TMin       " + format_number(block["TMin"]) + ";")
    lines.extend(["}", ""])
    return "\n".join(lines)


def _wall(profile) -> str:
    lines = [
        "unwrapRadius            " + format_number(profile.unwrap_radius) + ";",
        "circumferentialLength   " + format_number(profile.circumferential_length) + ";",
        "injectorArea            " + format_number(profile.injector_area) + ";",
        "recessDepth             " + format_number(profile.recess) + ";",
        "stations",
        "(",
    ]
    for segment in profile.segments:
        lines.append(
            "    ("
            + " ".join(
                format_number(value)
                for value in (
                    segment.y0,
                    segment.r_in0,
                    segment.r_out0,
                    segment.y1,
                    segment.r_in1,
                    segment.r_out1,
                )
            )
            + ")"
        )
    lines.append(");")
    lines.append("stepFaces")
    lines.append("(")
    for face in profile.step_faces:
        lines.append(
            "    { side "
            + face.side
            + "; location "
            + format_number(face.y)
            + "; r0 "
            + format_number(face.r0)
            + "; r1 "
            + format_number(face.r1)
            + "; }"
        )
    lines.append(");")
    return "\n".join(lines) + "\n"


def _geometry_report(profile, mesh) -> str:
    resolved_steps = "yes" if profile.radial_steps_resolved else "no"
    resolved_slots = "yes" if profile.injector_slots_resolved else "no"
    lines = [
        f"topology                 {mesh.topology};",
        f"dimension                {mesh.dimension};",
        f"radialStepsResolved      {resolved_steps};",
        f"injectorSlotsResolved    {resolved_slots};",
        "notes",
        "(",
    ]
    for note in profile.notes:
        lines.append(f'    "{note}"')
    lines.append(");")
    return "\n".join(lines) + "\n"


def _acoustic(acoustics) -> str:
    warning = acoustics.warning.replace('"', "'")
    return f"""ZReal              {format_number(acoustics.z_real)};
ZImag              {format_number(acoustics.z_imag)};
RReal              {format_number(acoustics.r_real)};
RImag              {format_number(acoustics.r_imag)};
magnitude          {format_number(acoustics.magnitude)};
phaseDeg           {format_number(acoustics.phase_deg)};
absorption         {format_number(acoustics.absorption)};
timeDomainModel    {acoustics.time_domain_model};
inertia            {format_number(acoustics.inertia)};
resistance         {format_number(acoustics.resistance)};
stiffness          {format_number(acoustics.stiffness)};
porosity           {format_number(acoustics.porosity)};
warning            "{warning}";
"""


def _fields(zero, config, state, mesh, acoustics) -> None:
    raw = config["raw"]
    ignition = raw["ignition"]
    operating = raw["operating"]
    p_internal = _scalar_internal(
        float(ignition["initialPressure"]),
        ignition,
        "TKernel",
        False,
    )
    # Pressure is not zoned. Temperature and progress are.
    write_file(zero / "p", "volScalarField", "p", _scalar_field(
        "[pressure]",
        f"uniform {format_number(ignition['initialPressure'])}",
        _pressure_boundaries(config, state),
    ))
    t_internal = _zoned_scalar(
        float(ignition["initialTemperature"]),
        ignition,
        "TKernel",
        ignition["type"] in ("hotKernel", "hotProgressKernel"),
    )
    write_file(zero / "T", "volScalarField", "T", _scalar_field(
        "[temperature]",
        t_internal,
        _temperature_boundaries(config),
    ))
    initial = {
        "circumferential": float(ignition["initialVelocity"]["circumferential"]),
        "axial": float(ignition["initialVelocity"]["axial"]),
        "radial": float(ignition["initialVelocity"]["radial"]),
    }
    internal_u = mesh_velocity(initial, mesh.topology, mesh.dimension, False)
    inlet_u = config["_inletVelocity"]
    write_file(zero / "U", "volVectorField", "U", _vector_field(
        internal_u,
        _velocity_boundaries(config, mesh, inlet_u, acoustics),
    ))
    lambda_internal = _zoned_scalar(
        float(ignition["initialLambda"]),
        ignition,
        "lambdaKernel",
        ignition["type"] in ("progressKernel", "hotProgressKernel"),
    )
    write_file(zero / "lambda", "volScalarField", "lambda", _scalar_field(
        "[0 0 0 0 0 0 0]",
        lambda_internal,
        _lambda_boundaries(config, mesh),
    ))
    _ = p_internal
    _ = operating


def _collect_boundary_names(mesh) -> list[str]:
    names = []
    for hex_block in mesh.hexes:
        for name in hex_block.faces:
            if name not in names:
                names.append(name)
    return names


def _scalar_internal(value, ignition, key, enabled):
    return f"uniform {format_number(value)}"


def _zoned_scalar(default, ignition, value_key, enabled) -> str:
    if not enabled:
        return f"uniform {format_number(default)}"
    box0 = " ".join(format_number(v) for v in ignition["boxMin"])
    box1 = " ".join(format_number(v) for v in ignition["boxMax"])
    return f"""{{
    type            zonal;
    defaultValue    {format_number(default)};
    zones
    {{
        kernel
        {{
            type            box;
            box             ({box0}) ({box1});
            value           {format_number(ignition[value_key])};
        }}
    }}
}}"""


def _scalar_field(dimensions, internal, boundaries: dict) -> str:
    lines = [
        f"dimensions      {dimensions};",
        "",
        "internalField   " + _inline(internal) + ";",
        "",
        "boundaryField",
        "{",
    ]
    for name, body in boundaries.items():
        lines.append(f"    {name}")
        lines.append("    {")
        lines.append(body)
        lines.append("    }")
    lines.append("}")
    return "\n".join(lines) + "\n"


def _inline(internal: str) -> str:
    if internal.startswith("{"):
        return internal
    return internal


def _vector_field(internal, boundaries) -> str:
    vector = "(" + " ".join(format_number(c) for c in internal) + ")"
    return _scalar_field("[velocity]", f"uniform {vector}", boundaries)


def _pressure_boundaries(config, state) -> dict:
    raw = config["raw"]
    names = _boundary_names(config)
    result = {}
    for name in names:
        if name in ("front", "back"):
            result[name] = "        type            empty;"
        elif name in ("periodicLeft", "periodicRight"):
            result[name] = "        type            cyclic;"
        elif name == "injector":
            result[name] = f"        type            {raw['boundaries']['inlet']['p']};"
        elif name == "outlet":
            result[name] = _outlet_pressure(config, state)
        else:
            result[name] = "        type            zeroGradient;"
    return result


def _outlet_pressure(config, state) -> str:
    outlet = config["raw"]["boundaries"]["outlet"]
    operating = config["raw"]["operating"]
    if outlet["type"] == "waveTransmissive":
        return "\n".join([
            "        type            waveTransmissive;",
            f"        gamma           {format_number(state.gamma)};",
            f"        fieldInf        {format_number(operating['pOutlet'])};",
            f"        lInf            {format_number(outlet['lInf'])};",
            f"        value           uniform {format_number(operating['pOutlet'])};",
        ])
    if outlet["type"] == "fixedValue":
        return "\n".join([
            "        type            fixedValue;",
            f"        value           uniform {format_number(operating['pOutlet'])};",
        ])
    return "        type            zeroGradient;"


def _temperature_boundaries(config) -> dict:
    raw = config["raw"]
    result = {}
    for name in _boundary_names(config):
        if name in ("front", "back"):
            result[name] = "        type            empty;"
        elif name in ("periodicLeft", "periodicRight"):
            result[name] = "        type            cyclic;"
        elif name == "injector":
            kind = raw["boundaries"]["inlet"]["T"]
            if kind == "fixedValue":
                result[name] = "\n".join([
                    "        type            fixedValue;",
                    f"        value           uniform {format_number(raw['operating']['TInlet'])};",
                ])
            else:
                result[name] = f"        type            {kind};"
        elif name == "outlet" and raw["boundaries"]["outlet"]["type"] == "fixedValue":
            result[name] = "\n".join([
                "        type            fixedValue;",
                f"        value           uniform {format_number(raw['operating']['TOutlet'])};",
            ])
        else:
            result[name] = "        type            zeroGradient;"
    return result


def _velocity_boundaries(config, mesh, inlet_u, acoustics) -> dict:
    raw = config["raw"]
    targets = {str(name) for name in raw["boundaries"]["acoustics"]["targetPatches"]}
    acoustic = raw["boundaries"]["walls"]["velocityModel"] == "acoustic"
    result = {}
    vector = "(" + " ".join(format_number(c) for c in inlet_u) + ")"
    for name in _boundary_names(config):
        if name in ("front", "back"):
            result[name] = "        type            empty;"
        elif name in ("periodicLeft", "periodicRight"):
            result[name] = "        type            cyclic;"
        elif name == "injector":
            kind = raw["boundaries"]["inlet"]["U"]
            if kind == "fixedValue":
                result[name] = "\n".join([
                    "        type            fixedValue;",
                    f"        value           uniform {vector};",
                ])
            else:
                result[name] = f"        type            {kind};"
        elif name == "outlet":
            result[name] = "        type            zeroGradient;"
        elif acoustic and name in targets:
            result[name] = _acoustic_bc(config, acoustics)
        else:
            result[name] = "        type            slip;"
    missing = targets - set(result)
    if missing and acoustic:
        raise ConfigError(
            "boundaries.acoustics.targetPatches names patches that this mesh "
            "does not have: " + ", ".join(sorted(missing))
        )
    return result


def _acoustic_bc(config, acoustics) -> str:
    model = config["raw"]["boundaries"]["acoustics"]["timeDomainModel"]
    micro = config["raw"]["boundaries"]["acoustics"]["microstructure"]
    p_ref = format_number(config["raw"]["boundaries"]["acoustics"]["pRef"])
    lines = [
        "        type            acousticWallVelocity;",
        f"        model           {model};",
        "        p               p;",
        f"        pRef            {p_ref};",
        f"        value           uniform (0 0 0);",
    ]
    if model == "resistive":
        lines.append(f"        ZReal           {format_number(acoustics.z_real)};")
        lines.append(
            "        resistiveFloor  "
            + format_number(config["raw"]["boundaries"]["acoustics"]["resistiveFloor"])
            + ";"
        )
    elif model == "helmholtz":
        lines.append(f"        inertia         {format_number(acoustics.inertia)};")
        lines.append(f"        resistance      {format_number(acoustics.resistance)};")
        lines.append(f"        stiffness       {format_number(acoustics.stiffness)};")
        lines.append(f"        porosity        {format_number(acoustics.porosity)};")
        lines.append(
            "        initialDisplacement "
            + format_number(micro["initialDisplacement"])
            + ";"
        )
        lines.append(
            "        initialVelocity " + format_number(micro["initialVelocity"]) + ";"
        )
    return "\n".join(lines)


def _lambda_boundaries(config, mesh) -> dict:
    inlet = format_number(config["raw"]["ignition"]["inletLambda"])
    result = {}
    for name in _boundary_names(config):
        if name in ("front", "back"):
            result[name] = "        type            empty;"
        elif name in ("periodicLeft", "periodicRight"):
            result[name] = "        type            cyclic;"
        elif name == "injector":
            result[name] = "\n".join([
                "        type            fixedValue;",
                f"        value           uniform {inlet};",
            ])
        else:
            result[name] = "        type            zeroGradient;"
    return result


def _boundary_names(config) -> list[str]:
    return list(config["_boundaryNames"])


def _scripts(case_dir: Path) -> None:
    allrun = case_dir / "Allrun"
    allclean = case_dir / "Allclean"
    allrun.write_text(
        "#!/bin/sh\ncd \"$(dirname \"$0\")\" || exit 1\nblockMesh\nfoamRun\n",
        encoding="utf-8",
    )
    allclean.write_text(
        "#!/bin/sh\ncd \"$(dirname \"$0\")\" || exit 1\nrm -rf constant/polyMesh postProcessing VTK/[0-9]* log.* processor*\n",
        encoding="utf-8",
    )
    allrun.chmod(0o755)
    allclean.chmod(0o755)
