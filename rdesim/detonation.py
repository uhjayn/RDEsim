"""Analytic detonation theory for the rdeFluid one-step model.

Standard library only. The model is one perfect gas (gamma, specific gas
constant R) with a progress variable lambda. Heat release q is added to the
energy equation as q * dlambda/dt and the rate is

    dlambda/dt = A (T/TRef)^beta exp(-Ta/T) (lambdaMax - lambda),   T >= Tign

with T floored at TFloor inside the exponential, as in rdeFluid. Everything
here is computed from those numbers. Nothing is fitted to a hidden constant.

Steady one-dimensional detonation in the shock frame, with gas speed u,
density rho and mass flux m = rho0 D:

    p + m u = P            h + u^2/2 - q lambda = H           h = gamma/(gamma-1) p/rho

Eliminating p and rho gives a quadratic in u for each lambda. The subsonic
root is the reacting flow behind the von Neumann (vN) shock.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass
class Gas:
    gamma: float
    R: float
    T0: float
    p0: float
    q: float

    @property
    def rho0(self) -> float:
        return self.p0 / (self.R * self.T0)

    @property
    def c0(self) -> float:
        return math.sqrt(self.gamma * self.R * self.T0)

    @property
    def cp(self) -> float:
        return self.gamma / (self.gamma - 1.0) * self.R

    @property
    def cv(self) -> float:
        return self.cp / self.gamma


@dataclass
class Chemistry:
    A: float
    beta: float
    Ta: float
    TRef: float
    TFloor: float
    Tign: float
    lambdaMax: float = 1.0

    def rate(self, T: float) -> float:
        """The coefficient of (lambdaMax - lambda). Zero below Tign."""
        if T < self.Tign:
            return 0.0
        safe = max(T, self.TFloor)
        return self.A * (safe / self.TRef) ** self.beta * math.exp(-self.Ta / safe)


def cj_state(gas: Gas, q: float | None = None) -> dict:
    """Chapman-Jouguet speed and burnt state for the one-gamma gas.

    p2/p0 = (1 + gamma M^2)/(gamma + 1) and rho2/rho0 = (gamma+1) M^2/(1 + gamma M^2)
    with M = D/c0. The heat release that goes with a given M follows from the
    energy equation, so the speed is found by inverting it.
    """
    q = gas.q if q is None else q

    def heat_for_mach(mach: float) -> float:
        g = gas.gamma
        pr = (1.0 + g * mach**2) / (g + 1.0)
        rr = (g + 1.0) * mach**2 / (1.0 + g * mach**2)
        d = mach * gas.c0
        u2 = d / rr
        t2 = gas.T0 * pr / rr
        return gas.cp * t2 + 0.5 * u2**2 - gas.cp * gas.T0 - 0.5 * d**2

    lo, hi = 1.0 + 1e-9, 50.0
    if heat_for_mach(hi) < q:
        raise ValueError("heat release too large for the bracket")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if heat_for_mach(mid) < q:
            lo = mid
        else:
            hi = mid
    mach = 0.5 * (lo + hi)
    g = gas.gamma
    pr = (1.0 + g * mach**2) / (g + 1.0)
    rr = (g + 1.0) * mach**2 / (1.0 + g * mach**2)
    return {
        "D": mach * gas.c0,
        "M": mach,
        "p": pr * gas.p0,
        "T": gas.T0 * pr / rr,
        "rho": rr * gas.rho0,
        "u": mach * gas.c0 * (1.0 - 1.0 / rr),
    }


def heat_for_cj_speed(gas: Gas, speed: float) -> float:
    """Heat release q (J/kg) that puts the CJ speed at ``speed`` for this gas."""
    g = gas.gamma
    mach = speed / gas.c0
    if mach <= 1.0:
        raise ValueError("target speed is not supersonic for the reactant gas")
    pr = (1.0 + g * mach**2) / (g + 1.0)
    rr = (g + 1.0) * mach**2 / (1.0 + g * mach**2)
    u2 = speed / rr
    t2 = gas.T0 * pr / rr
    return gas.cp * t2 + 0.5 * u2**2 - gas.cp * gas.T0 - 0.5 * speed**2


def von_neumann(gas: Gas, speed: float) -> dict:
    g = gas.gamma
    mach = speed / gas.c0
    pr = 2.0 * g * mach**2 / (g + 1.0) - (g - 1.0) / (g + 1.0)
    rr = (g + 1.0) * mach**2 / ((g - 1.0) * mach**2 + 2.0)
    return {
        "p": pr * gas.p0,
        "rho": rr * gas.rho0,
        "T": gas.T0 * pr / rr,
        "u": speed / rr,
    }


def _flow_state(gas: Gas, speed: float, lam: float, total_p: float, mass: float, total_h: float):
    """Subsonic-root gas speed, temperature and pressure at progress lambda."""
    g = gas.gamma
    disc = (g * total_p / mass) ** 2 - 2.0 * (g + 1.0) * (g - 1.0) * (total_h + gas.q * lam)
    if disc < 0.0:
        raise ValueError("no steady ZND solution: speed is below the CJ speed")
    u = (g * total_p / mass - math.sqrt(disc)) / (g + 1.0)
    p = total_p - mass * u
    t = p * u / (mass * gas.R)
    return u, t, p


def znd_structure(gas: Gas, chem: Chemistry, speed: float, points: int = 20001, s_max: float = 14.0) -> dict:
    """Induction length and time from the steady ZND profile.

    The integration variable is s = -ln(1 - lambda/lambdaMax), so the first-order
    tail does not diverge. Induction is where dT/dx peaks, the same definition
    used for the Cantera reference.
    """
    vn = von_neumann(gas, speed)
    mass = gas.rho0 * speed
    total_p = gas.p0 + mass * speed
    total_h = gas.cp * gas.T0 + 0.5 * speed**2
    lam_max = chem.lambdaMax

    if vn["T"] < chem.Tign:
        return {
            "vn": vn, "x_ind": math.inf, "t_ind": math.inf, "x50": math.inf,
            "t50": math.inf, "x90": math.inf, "peak": None,
        }

    ds = s_max / (points - 1)
    xs = [0.0]
    ts = [0.0]
    temps = []
    us = []
    fx = []
    ft = []
    for i in range(points):
        s = i * ds
        lam = lam_max * (1.0 - math.exp(-s))
        try:
            u, t, _ = _flow_state(gas, speed, lam, total_p, mass, total_h)
        except ValueError:
            u, t = us[-1], temps[-1]
        temps.append(t)
        us.append(u)
        r = chem.rate(t)
        if r <= 0.0:
            fx.append(math.inf)
            ft.append(math.inf)
        else:
            fx.append(u / r)
            ft.append(1.0 / r)
    for i in range(1, points):
        xs.append(xs[-1] + 0.5 * ds * (fx[i] + fx[i - 1]))
        ts.append(ts[-1] + 0.5 * ds * (ft[i] + ft[i - 1]))

    dtdx = []
    for i in range(points):
        s = i * ds
        lam = lam_max * (1.0 - math.exp(-s))
        j0, j1 = max(i - 1, 0), min(i + 1, points - 1)
        dlam = lam_max * (math.exp(-j0 * ds) - math.exp(-j1 * ds))
        dt_dlam = (temps[j1] - temps[j0]) / dlam if dlam != 0.0 else 0.0
        dlam_dx = chem.rate(temps[i]) * (lam_max - lam) / us[i]
        dtdx.append(dt_dlam * dlam_dx)
    peak = max(range(points), key=lambda i: dtdx[i])
    if 0 < peak < points - 1:
        y0, y1, y2 = dtdx[peak - 1], dtdx[peak], dtdx[peak + 1]
        denom = y0 - 2.0 * y1 + y2
        shift = 0.5 * (y0 - y2) / denom if denom != 0.0 else 0.0
        shift = max(-1.0, min(1.0, shift))
    else:
        shift = 0.0
    frac = peak + shift
    i0 = int(math.floor(frac))
    i1 = min(i0 + 1, points - 1)
    w = frac - i0
    x_ind = xs[i0] * (1 - w) + xs[i1] * w
    t_ind = ts[i0] * (1 - w) + ts[i1] * w

    def at(fraction: float, table: list[float]) -> float:
        s_target = -math.log(1.0 - fraction)
        k = s_target / ds
        k0 = int(math.floor(k))
        k1 = min(k0 + 1, points - 1)
        return table[k0] + (k - k0) * (table[k1] - table[k0])

    if peak == 0:
        # dT/dx is largest at the vN state: the rate has no thermal-runaway onset.
        x_ind = math.nan
        t_ind = math.nan

    return {
        "vn": vn, "x_ind": x_ind, "t_ind": t_ind,
        "x50": at(0.5, xs), "t50": at(0.5, ts), "x90": at(0.9, xs),
        "peak": {"dTdx": dtdx[peak], "T": temps[peak], "lambda": lam_max * (1 - math.exp(-peak * ds))},
    }


def reduced_activation_energy(
    gas: Gas, chem: Chemistry, speed: float, delta: float, key: str = "t50", points: int = 8001
) -> float:
    """Ea/(R T_vN) from ln(tau) against 1/T_vN over speeds speed*(1 -/+ delta).

    ``key`` picks the time: "t50" is the time to half the heat release, "t_ind"
    the time to the maximum of dT/dx.
    """
    lo = znd_structure(gas, chem, speed * (1.0 - delta), points=points)
    hi = znd_structure(gas, chem, speed * (1.0 + delta), points=points)
    dln = math.log(hi[key]) - math.log(lo[key])
    dinv = 1.0 / hi["vn"]["T"] - 1.0 / lo["vn"]["T"]
    ea_over_r = dln / dinv
    mid = von_neumann(gas, speed)
    return ea_over_r / mid["T"]


def _bisect(func, lo: float, hi: float, tol: float = 1e-12, limit: int = 200) -> float:
    flo = func(lo)
    fhi = func(hi)
    if flo * fhi > 0:
        raise ValueError(f"root not bracketed on [{lo}, {hi}]: f = {flo}, {fhi}")
    for _ in range(limit):
        mid = 0.5 * (lo + hi)
        fmid = func(mid)
        if flo * fmid <= 0:
            hi, fhi = mid, fmid
        else:
            lo, flo = mid, fmid
        if abs(hi - lo) < tol * max(abs(lo), abs(hi), 1e-300):
            break
    return 0.5 * (lo + hi)


def calibrate(gas: Gas, chem: Chemistry, target_speed: float, target_x50: float, target_theta: float, delta: float) -> dict:
    """Solve for q, Ta and A that meet the three targets.

    q          from the CJ speed (closed form, one gamma).
    Ta         bisection on the reduced activation energy at the CJ speed, from
               the time to half heat release.
    A          bisection on the distance to half heat release at the CJ speed.
    A is solved inside every Ta trial, because Ta moves that distance.
    """
    q = heat_for_cj_speed(gas, target_speed)
    cal_gas = Gas(gas.gamma, gas.R, gas.T0, gas.p0, q)
    cj = cj_state(cal_gas)
    speed = cj["D"]

    def solve_A(Ta: float) -> float:
        trial = Chemistry(1.0, chem.beta, Ta, chem.TRef, chem.TFloor, chem.Tign, chem.lambdaMax)

        def f(logA):
            c = Chemistry(math.exp(logA), trial.beta, trial.Ta, trial.TRef, trial.TFloor, trial.Tign, trial.lambdaMax)
            x = znd_structure(cal_gas, c, speed, points=4001)["x50"]
            if not math.isfinite(x):
                return 700.0
            return (math.log(x) if x > 0.0 else -700.0) - math.log(target_x50)

        return math.exp(_bisect(f, math.log(1e0), math.log(1e40), tol=1e-11))

    t_vn = von_neumann(cal_gas, speed)["T"]

    def theta_error(Ta: float) -> float:
        c = Chemistry(solve_A(Ta), chem.beta, Ta, chem.TRef, chem.TFloor, chem.Tign, chem.lambdaMax)
        return reduced_activation_energy(cal_gas, c, speed, delta, "t50", 4001) - target_theta

    ta = _bisect(theta_error, 0.5 * target_theta * t_vn, 2.0 * target_theta * t_vn, tol=1e-9)
    a = solve_A(ta)
    return {"q": q, "A": a, "Ta": ta, "speed": speed}
