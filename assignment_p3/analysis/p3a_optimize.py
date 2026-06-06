"""ET4260 Part 3a.3 -- dimension optimization of the phase shifter.

Maximize the achievable peak low-frequency modulation sensitivity
    S_max = (2 pi L / lambda) * dn_eff * exp(-b g0 - 1/2) / V_DC_MAX
subject to the assignment constraints
    f_r > 500 kHz,  V_PI > 200 V,  V_DC_MAX < V_PI,  a0_max(w_r) > 2e4 m/s^2
and the engineering constraints that Part 2's FEM taught us are binding:
    - side (lateral) pull-in of fingers + shuttle rocking (calibrated x0.90)
    - bar/spine tip rigidity (internal modes >= 2.5 f_r)
    - comb must fit under the spine of length L
    - features in [200 nm, 1 mm], end fingers 3Wc >= 600 nm, 2 < N < 200.

Strategy: vectorized grid search using a precomputed shape function for the
finger side pull-in, then exact re-evaluation (p3a_model.evaluate) of the
best candidates, then manual rounding to clean values.
"""
import numpy as np
from itertools import product
import p3a_model as M

EPS0, E_SI, RHO, T = M.EPS0, M.E_SI, M.RHO, M.T_SI
B, LAM, DN = M.B_OPT, M.LAMBDA, M.DNEFF
X_OPT = M.X_OPT

F_TARGET = 540e3           # analytic target -> FEM ~ 518 kHz (x0.96)
MARGIN_PI = 1.25           # min pull-in over operating voltage
V_PI_MIN = 230.0           # >200 V with FEM margin
A_RES_MIN = 2.2e4          # >2e4 with margin
MODE_SEP = 2.5             # secondary modes >= 2.5 f_r


# ---------- precompute side pull-in shape function Phi(a/L) ----------
# V_side^2 = E I h^3 / (2 eps0 t L^4) * Phi(a/L),  I = t Wc^3/12
def _phi(ahat, n=400):
    xi = np.linspace(0, 1, n)
    xa = np.linspace(ahat, 1, n)
    w = np.empty(n)
    for i, x in enumerate(xi):
        c = np.where(x <= xa, x**2 * (3 * xa - x), xa**2 * (3 * x - xa))
        w[i] = np.trapezoid(c, xa) / 6.0
    wpp = np.gradient(np.gradient(w, xi), xi)
    num = np.trapezoid(wpp**2, xi)
    den = np.trapezoid(np.where(xi >= ahat, w**2, 0.0), xi)
    return num / den


AHAT = np.linspace(0.0, 0.95, 96)
PHI = np.array([_phi(a) for a in AHAT])


def v_side_fast(Lc, Wc, h, d):
    ahat = np.clip(1.0 - d / Lc, 0.0, 0.95)
    phi = np.interp(ahat, AHAT, PHI)
    I = T * Wc**3 / 12.0
    return np.sqrt(E_SI * I * h**3 / (2 * EPS0 * T * Lc**4) * phi)


def u_tip(gt0, h, Wc):
    """Vectorized Newton solve of g_t0 = 1.5 u + u^3/(2 h Wc)."""
    c = 1.0 / (2 * h * Wc)
    u = np.minimum(gt0 / 1.5, (gt0 / c) ** (1 / 3))
    for _ in range(60):
        f = 1.5 * u + c * u**3 - gt0
        u = u - f / (1.5 + 3 * c * u**2)
    return u


def design(L, Wc, h, Lc, d0, Lm, Lb, Wm):
    """Build a full candidate from the scan variables; returns dict or None."""
    # bar/spine tip-rigidity floor: f_bar = (3.516/2pi) c W/(ell^2) >= 2.5 f_r
    cbar = np.sqrt(E_SI / (12 * RHO))
    ell = (L - Wm) / 2
    f_bar = (3.516 / (2 * np.pi)) * cbar * Wm / ell**2
    if f_bar < MODE_SEP * F_TARGET or Wm > 0.4 * L:
        return None

    # comb fill
    p = 2 * (Wc + h)
    n_mov = int(min(np.floor((0.98 * L - 6 * Wc) / p), 99))
    N = 2 * n_mov
    if N < 4:
        return None

    # solve Wb for k = m_eff w^2 (fixed point on beam mass)
    w2 = (2 * np.pi * F_TARGET) ** 2
    m_fix = RHO * T * (2 * L * Wm + Wm * Lm + n_mov * Lc * Wc)
    Wb = 2e-6
    for _ in range(8):
        m_eff = m_fix + M.BEAM_MASS_EFF * RHO * T * 4 * Lb * Wb
        k_req = m_eff * w2
        Wb = (k_req * Lb**3 / (4 * E_SI * T)) ** (1 / 3)
    if not (1e-6 <= Wb <= 20e-6):
        return None
    k = 4 * E_SI * T * Wb**3 / Lb**3

    # voltages
    V_max = np.sqrt(k * h / (B * N * EPS0 * T))
    gt0 = Lc - d0
    if gt0 <= 0:
        return None
    u = u_tip(gt0, h, Wc)
    x_pi = gt0 - u
    if x_pi < 1.2 * X_OPT:                        # travel margin
        return None
    V_pi = np.sqrt(k * u**3 / (N * EPS0 * Wc * T))
    d_op = d0 + X_OPT
    kneg = 2 * EPS0 * T * d_op / h**3
    V_sf = v_side_fast(Lc, Wc, h, d_op)
    k_ax = E_SI * T * Wb / Lb
    k_tr = E_SI * T * Wb**3 / Lb**3
    s_pair = max(Lm - Wb - 2e-6, 1e-6)
    I_col = T * Wm**3 / 12.0
    k_col = 12 * E_SI * I_col / s_pair**3
    k_lvl = 1.0 / (1.0 / (2 * k_ax) + 1.0 / k_col)
    K_th = 2 * k_lvl * (s_pair / 2)**2 + 4 * k_tr * (Wm / 2)**2
    sum_x2 = p**2 * n_mov * (n_mov**2 - 1) / 12.0
    V_rk = np.sqrt(K_th / (kneg * sum_x2)) if sum_x2 > 0 else np.inf
    I_spine = T * Wm**3 / 12.0
    ell_s = (L - Wm) / 2
    V_sp = np.sqrt((1.875 / ell_s)**4 * E_SI * I_spine * p / kneg)
    V_lt = np.sqrt(4 * k_ax / (n_mov * kneg))
    V_side = M.V_SIDE_CAL / np.sqrt(1 / V_sf**2 + 1 / V_rk**2
                                    + 1 / V_sp**2 + 1 / V_lt**2)
    V_PI = min(V_pi, V_side)
    if V_PI < max(V_PI_MIN, MARGIN_PI * V_max):
        return None

    # damping/acceleration constraint (V=0 worst case)
    A_m = 2 * L * Wm + Wm * Lm + n_mov * Lc * Wc + 2 * Lb * Wb
    b_air = M.C_AIR * A_m
    Q = np.sqrt(k * m_eff) / (b_air + np.sqrt(k * m_eff) / M.Q_CLAMP)
    a_res = 200e-9 * w2 / Q
    if a_res < A_RES_MIN:
        return None

    # secondary modes
    f_f = (3.516 / (2 * np.pi)) * cbar * Wc / Lc**2
    f_b = (22.37 / (2 * np.pi)) * cbar * Wb / Lb**2
    if f_f < MODE_SEP * F_TARGET or f_b < MODE_SEP * F_TARGET:
        return None

    S = (2 * np.pi * L / LAM) * DN * np.exp(-B * 200e-9 - 0.5) / V_max
    return dict(L=L, Wc=Wc, h=h, Lc=Lc, d0=d0, Lm=Lm, Lb=Lb, Wm=Wm, Wb=Wb,
                N=N, k=k, V_max=V_max, V_PI=V_PI, V_pi_ax=V_pi,
                V_side=V_side, V_sf=V_sf, V_rk=V_rk, V_sp=V_sp, V_lt=V_lt,
                Q=Q, a_res=a_res, S=S)


def main(grids=None):
    if grids is None:
        grids = dict(
            L=np.arange(100, 301, 25) * 1e-6,
            Wc=np.arange(0.4, 1.21, 0.2) * 1e-6,
            h=np.arange(0.4, 1.21, 0.1) * 1e-6,
            Lc=np.array([5, 6, 7, 8]) * 1e-6,
            d0=np.array([1.5, 2.5]) * 1e-6,
            Lm=np.array([20, 40]) * 1e-6,
            Lb=np.array([40, 60, 80]) * 1e-6,
            Wm=np.array([4, 6, 8, 10, 12, 15, 18, 22, 26, 30]) * 1e-6,
        )
    best = []
    for L, Wc, h, Lc, d0, Lm, Lb, Wm in product(*grids.values()):
        d = design(L, Wc, h, Lc, d0, Lm, Lb, Wm)
        if d:
            best.append(d)
    best.sort(key=lambda d: -d['S'])
    print(f"feasible designs: {len(best)}")
    for d in best[:15]:
        print(f"S={d['S']*1e3:7.2f} mrad/V | L={d['L']*1e6:5.0f} Wm={d['Wm']*1e6:5.1f} "
              f"Wc={d['Wc']*1e6:4.2f} h={d['h']*1e6:4.2f} Lc={d['Lc']*1e6:4.1f} "
              f"d0={d['d0']*1e6:3.1f} Lm={d['Lm']*1e6:3.0f} Lb={d['Lb']*1e6:3.0f} "
              f"Wb={d['Wb']*1e6:4.2f} N={d['N']:4d} | k={d['k']:6.2f} "
              f"Vmax={d['V_max']:5.0f} Vpi_ax={d['V_pi_ax']:5.0f} Vside={d['V_side']:5.0f} "
              f"(f{d['V_sf']:5.0f}/r{d['V_rk']:5.0f}/s{d['V_sp']:5.0f}/l{d['V_lt']:5.0f}) "
              f"Q={d['Q']:5.1f} a_res={d['a_res']:.2e}")
    return best


if __name__ == '__main__':
    main()
