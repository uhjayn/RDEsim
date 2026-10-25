#!/usr/bin/env python3
"""Reference detonation numbers for stoichiometric H2-air from Cantera and NASA CEA.

This is the only place the reference values are computed. It needs packages the
rest of RDEsim does not (cantera, scipy, numpy; cea is optional). Run it in a
scratch virtualenv:

    python3 -m venv /tmp/refenv
    /tmp/refenv/bin/pip install cantera scipy numpy cea
    /tmp/refenv/bin/python tools/validation_reference.py --write

Without --write the numbers are printed and nothing is saved. With --write the
whole file cases/detonationTubeH2/referenceData is rewritten. The cell-width
block below (LITERATURE_CELL_WIDTH) is transcribed by hand from the Shepherd
detonation database data files named in each entry; the script does not
download or recompute it.

Definitions (also used for the one-step model, see rdesim/detonation.py):
  CJ        equilibrium Chapman-Jouguet state, minimum-speed tangency of the
            Rayleigh line and the equilibrium Hugoniot.
  vN        frozen-composition shock at the CJ speed.
  induction distance from the vN state to the maximum of dT/dx in the steady
            ZND profile. The distance to the maximum of thermicity is also given.
  theta     reduced effective activation energy, Ea/(R T_vN) =
            (1/T_vN) d ln(tau_ind) / d(1/T_vN), by central difference over
            shock speeds (1 -/+ delta) D_CJ. tau_ind is the induction time.
  half      distance and time to half the chemical heat release, from the
            formation-enthalpy sum of the species, normalised by the CJ
            equilibrium products (thetaHalf is theta for that time).
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import cantera as ct
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "cases" / "detonationTubeH2" / "referenceData"

T0 = 300.0
P0 = 101325.0
X0 = "H2:2,O2:1,N2:3.76"
DELTA = 0.02
RU = ct.gas_constant


def cj(mech: str):
    g1 = ct.Solution(mech)
    g1.TPX = T0, P0, X0
    r1 = g1.density
    v1 = 1.0 / r1
    h1 = g1.enthalpy_mass
    g2 = ct.Solution(mech)

    def resid(k, D):
        v2 = v1 / k
        p = P0 + r1**2 * D**2 * (v1 - v2)
        h = h1 + 0.5 * D**2 * (1.0 - (v2 / v1) ** 2)
        g2.TPX = 2500.0, p, g1.X
        g2.HP = h, p
        g2.equilibrate("HP")
        return (1.0 / g2.density - v2) * r1

    def gmin(D):
        r = minimize_scalar(
            lambda k: resid(k, D), bounds=(1.5, 2.3), method="bounded",
            options={"xatol": 1e-10},
        )
        return r.fun, r.x

    D = brentq(lambda d: gmin(d)[0], 1900.0, 2100.0, xtol=1e-7)
    resid(gmin(D)[1], D)
    return g1, g2, D


def frozen_shock(g1, D, mech):
    """Frozen-composition shock at speed D. Returns a Solution at the post-shock state."""
    r1 = g1.density
    h1 = g1.enthalpy_mass
    p1 = g1.P
    g = ct.Solution(mech)
    g.TPX = g1.T, g1.P, g1.X
    wbar = g1.mean_molecular_weight

    def resid(r2):
        p2 = p1 + r1 * D**2 * (1.0 - r1 / r2)
        t2 = p2 * wbar / (r2 * RU)
        g.TPX = t2, p2, g1.X
        return g.enthalpy_mass - (h1 + 0.5 * D**2 * (1.0 - (r1 / r2) ** 2))

    r2 = brentq(resid, 1.2 * r1, 12.0 * r1, xtol=1e-12)
    resid(r2)
    return g, r1 * D


def znd(g1, D, mech, g_end, x_end=None, return_profile=False):
    gv, m = frozen_shock(g1, D, mech)
    p_tot = gv.P + m * m / gv.density
    u0 = m / gv.density
    gas = ct.Solution(mech)
    gas.TDY = gv.T, gv.density, gv.Y
    n = gas.n_species
    wk = gas.molecular_weights

    def rhs(x, y):
        u = y[0]
        T = y[1]
        Y = np.maximum(y[2:], 0.0)
        rho = m / u
        gas.TDY = T, rho, Y
        wdot = gas.net_production_rates
        dYdx = wdot * wk / (rho * u)
        hk = gas.partial_molar_enthalpies / wk
        cp = gas.cp_mass
        wbar = gas.mean_molecular_weight
        dYdt = dYdx * u
        sig = np.sum((wbar / wk - hk / (cp * T)) * dYdt)
        a2 = gas.cp_mass / gas.cv_mass * ct.gas_constant * T / wbar
        M2 = u * u / a2
        dudx = sig / (1.0 - M2)
        dTdx = -(u * dudx + np.sum(hk * dYdx)) / cp
        out = np.empty(n + 2)
        out[0] = dudx
        out[1] = dTdx
        out[2:] = dYdx
        return out

    y0 = np.concatenate(([u0, gv.T], gv.Y))
    x_max = x_end if x_end is not None else 5e-3

    def stop(x, y):
        u = y[0]
        gas.TDY = y[1], m / u, np.maximum(y[2:], 0.0)
        a2 = gas.cp_mass / gas.cv_mass * ct.gas_constant * y[1] / gas.mean_molecular_weight
        return 1.0 - u * u / a2 - 2e-3

    stop.terminal = True
    stop.direction = -1
    sol = solve_ivp(
        rhs, (0.0, x_max), y0, method="BDF", rtol=1e-8, atol=1e-12,
        events=stop, dense_output=True, first_step=1e-10,
    )
    xs = np.concatenate((np.linspace(0, 1, 2001)[1:] ** 4 * sol.t[-1],))
    xs = np.concatenate(([0.0], xs))
    ys = sol.sol(xs)
    u = ys[0]
    T = ys[1]
    Y = np.maximum(ys[2:].T, 0.0)
    rho = m / u
    p = p_tot - m * u
    dT = np.gradient(T, xs)
    sig = np.zeros_like(xs)
    for i in range(len(xs)):
        gas.TDY = T[i], rho[i], Y[i]
        wdot = gas.net_production_rates
        dYdt = wdot * wk / rho[i]
        hk = gas.partial_molar_enthalpies / wk
        sig[i] = np.sum((gas.mean_molecular_weight / wk - hk / (gas.cp_mass * T[i])) * dYdt)
    t = np.concatenate(([0.0], np.cumsum(0.5 * (1 / u[1:] + 1 / u[:-1]) * np.diff(xs))))
    i_T = int(np.argmax(dT))
    i_s = int(np.argmax(sig))
    lam = heat_release_progress(dict(Y=Y), mech, g_end)
    res = dict(
        D=D, Tvn=gv.T, pvn=gv.P, rhovn=gv.density, uvn=u0,
        x_ind_T=float(xs[i_T]), t_ind_T=float(t[i_T]),
        x_ind_sigma=float(xs[i_s]), t_ind_sigma=float(t[i_s]),
        x50=float(np.interp(0.5, lam, xs)), t50=float(np.interp(0.5, lam, t)),
        x90=float(np.interp(0.9, lam, xs)),
        terminated=bool(sol.t_events[0].size),
    )
    if return_profile:
        res["profile"] = dict(x=xs, T=T, p=p, u=u, Y=Y, rho=rho, sig=sig, lam=lam)
    return res


def heat_release_progress(profile, mech, g_end):
    """Progress of chemical heat release from the formation-enthalpy sum."""
    g0 = ct.Solution(mech)
    g0.TPX = 298.15, ct.one_atm, X0
    hf = g0.partial_molar_enthalpies / g0.molecular_weights
    chem = profile["Y"] @ hf
    chem_end = float(g_end.Y @ hf)
    return (chem[0] - chem) / (chem[0] - chem_end)


def theta_eff(g1, D_cj, mech, g_end):
    """Ea/R in kelvin from the time to max dT/dx and from the time to half heat release."""
    lo = znd(g1, D_cj * (1 - DELTA), mech, g_end)
    hi = znd(g1, D_cj * (1 + DELTA), mech, g_end)
    dinv = 1.0 / hi["Tvn"] - 1.0 / lo["Tvn"]
    ea_ind = (math.log(hi["t_ind_T"]) - math.log(lo["t_ind_T"])) / dinv
    ea_50 = (math.log(hi["t50"]) - math.log(lo["t50"])) / dinv
    return ea_ind, ea_50


def cea_cj():
    try:
        import cea
    except ImportError:
        return None
    reac_names = ["H2", "O2", "N2"]
    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(reac_names, products_from_reactants=True)
    solver = cea.DetonationSolver(prod, reactants=reac)
    sol = cea.DetonationSolution(solver)
    w_h2, w_o2, w_n2 = 2.01588, 31.998, 28.0134
    mass = np.array([2.0 * w_h2, 1.0 * w_o2, 3.76 * w_n2])
    solver.solve(sol, mass / mass.sum(), T0, P0 / 1e5)
    return dict(
        D=float(sol.velocity), p=float(sol.P) * 1e5, T=float(sol.T),
        rho=float(sol.density), version=cea.__version__,
    )


# Transcribed from https://shepherd.caltech.edu/detn_db/data/plotdata/<file>
# (listed on https://shepherd.caltech.edu/detn_db/html/db_121.html, "Cell Width - H2 Fuel").
# Only H2-air with no diluent, 100-101.3 kPa, 293-300 K, phi at or next to 1.
# (file, phi, P kPa, T K, cell width mm, original reference)
LITERATURE_CELL_WIDTH = [
    ("at14.txt", 1.0233, 100.0, 300.0, 8.1875,
     "Ciccarelli, Ginsberg, Boccio, Economos, Sato, Kinoshita, Combust. Flame 99(2):212-220, 1994"),
    ("ja20a.txt", 1.0, 100.0, 300.0, 8.0,
     "Ciccarelli et al., NUREG/CR-6391 (BNL-NUREG-52482), Brookhaven National Laboratory, 1997"),
    ("H2-Air1.txt", 1.0, 101.3, 293.0, 15.1,
     "Guirao, Knystautas, Lee, Benedick, Berman, 19th Symp. (Int.) Combust., pp. 583-590, 1982"),
    ("at5a.txt", 1.0, 100.834, 293.0, 11.0429,
     "Stamps and Tieszen, Combust. Flame 83(3):353-364, 1991"),
    ("at5a.txt", 1.0, 101.452, 293.0, 9.2068,
     "Stamps and Tieszen, Combust. Flame 83(3):353-364, 1991"),
]


def const_p_delay(mech, T, p):
    g = ct.Solution(mech)
    g.TPX = T, p, X0
    reactor = ct.IdealGasConstPressureReactor(g)
    net = ct.ReactorNet([reactor])
    t, temp = [], []
    while net.time < 1e-5:
        net.step()
        t.append(net.time)
        temp.append(reactor.T)
    t = np.array(t)
    return float(t[np.argmax(np.gradient(np.array(temp), t))])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()

    lines = []
    out = lines.append
    out("// Reference data for the stoichiometric H2-air detonation tube.")
    out("// Written by tools/validation_reference.py. Edit that script, not this file.")
    out(f"// cantera {ct.__version__}; numpy {np.__version__}")
    out("//")
    out("// computed.<mech>   Cantera: equilibrium CJ state, frozen vN state, steady ZND profile.")
    out("// computed.cea      NASA CEA CJ state: an independent check of the CJ numbers only.")
    out("// literature        Cell width, transcribed from the Shepherd detonation database.")
    out("// Not available and not invented: a measured induction length, a measured CJ pressure.")
    out("")
    out("computed")
    out("{")
    out(f"    initialTemperature    {T0:g};")
    out(f"    initialPressure       {P0:g};")
    out('    composition           "H2:2, O2:1, N2:3.76 (mole ratio; phi = 1)";')
    out(f"    thetaShockSpeedDelta  {DELTA:g};")

    primary = "gri30.yaml"
    for mech, tag in (("gri30.yaml", "gri30"), ("h2o2.yaml", "h2o2")):
        g1, g2, D = cj(mech)
        d = znd(g1, D, mech, g2, return_profile=True)
        ea_r, ea_50 = theta_eff(g1, D, mech, g2)
        th = ea_r / d["Tvn"]
        th50 = ea_50 / d["Tvn"]
        out(f"    {tag}")
        out("    {")
        out(f"        mechanism         {mech};")
        out(f"        rho0              {g1.density:.6g};")
        out(f"        c0                {g1.sound_speed:.6g};")
        out(f"        W0                {g1.mean_molecular_weight:.6g};")
        out(f"        gamma0            {g1.cp_mass / g1.cv_mass:.6g};")
        out(f"        DCJ               {D:.6g};")
        out(f"        pCJ               {g2.P:.6g};")
        out(f"        TCJ               {g2.T:.6g};")
        out(f"        rhoCJ             {g2.density:.6g};")
        out(f"        WCJ               {g2.mean_molecular_weight:.6g};")
        out(f"        pvN               {d['pvn']:.6g};")
        out(f"        TvN               {d['Tvn']:.6g};")
        out(f"        rhovN             {d['rhovn']:.6g};")
        out(f"        xIndTmax          {d['x_ind_T']:.6g};")
        out(f"        tIndTmax          {d['t_ind_T']:.6g};")
        out(f"        xIndThermicity    {d['x_ind_sigma']:.6g};")
        out(f"        tIndThermicity    {d['t_ind_sigma']:.6g};")
        out(f"        EaOverR           {ea_r:.6g};")
        out(f"        theta             {th:.6g};")
        out(f"        tIndConstPReactor {const_p_delay(mech, d['Tvn'], d['pvn']):.6g};")
        out(f"        EaOverRHalf       {ea_50:.6g};")
        out(f"        thetaHalf         {th50:.6g};")
        out(f"        xHeatRelease50    {d['x50']:.6g};")
        out(f"        tHeatRelease50    {d['t50']:.6g};")
        out(f"        xHeatRelease90    {d['x90']:.6g};")
        out("    }")
        print(
            f"{tag}: D_CJ={D:.2f} m/s p_CJ={g2.P/101325:.3f} atm T_CJ={g2.T:.1f} K "
            f"| vN p={d['pvn']/101325:.2f} atm T={d['Tvn']:.1f} K "
            f"| x_ind(max dT/dx)={d['x_ind_T']*1e3:.4f} mm t={d['t_ind_T']*1e6:.3f} us "
            f"x_ind(max thermicity)={d['x_ind_sigma']*1e3:.4f} mm "
            f"| Ea/R={ea_r:.0f} K theta={th:.3f} theta50={th50:.3f} "
            f"| x50={d['x50']*1e3:.4f} mm x90={d['x90']*1e3:.4f} mm (terminated={d['terminated']})"
        )
    c = cea_cj()
    if c is not None:
        out("    cea")
        out("    {")
        out(f"        version           {c['version']};")
        out(f"        DCJ               {c['D']:.6g};")
        out(f"        pCJ               {c['p']:.6g};")
        out(f"        TCJ               {c['T']:.6g};")
        out(f"        rhoCJ             {c['rho']:.6g};")
        out("    }")
        print(f"cea {c['version']}: D_CJ={c['D']:.2f} p_CJ={c['p']/101325:.3f} atm T_CJ={c['T']:.1f} K")
    out("}")
    out("")
    out("literature")
    out("{")
    out("    database    \"M. Kaneshige and J.E. Shepherd, Detonation database, GALCIT Report FM97-8, 1997;")
    out("                 https://shepherd.caltech.edu/detn_db/html/db_121.html (Cell Width - H2 Fuel)\";")
    out("    cellWidthUnit    mm;")
    out("    cellWidth")
    out("    (")
    for fname, phi, pk, tk, width, ref in LITERATURE_CELL_WIDTH:
        out(f"        {{ file {fname}; phi {phi:g}; pressureKPa {pk:g}; temperatureK {tk:g}; width {width:g}; reference \"{ref}\"; }}")
    out("    );")
    out("}")
    text = "\n".join(lines) + "\n"
    if args.write:
        OUT.parent.mkdir(parents=True, exist_ok=True)
        OUT.write_text(text, encoding="utf-8")
        print(f"Wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
