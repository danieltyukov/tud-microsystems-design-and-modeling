"""ET4260 Part 3a -- analytical model of the comb-drive optical phase shifter.

Lumped-element model from Part 1, extended with the side-instability physics
learned from the Part-2 FEM (rocking/finger lateral pull-in), used both to
derive V_DC_MAX (a.1), the reachability constraints (a.2) and to drive the
dimension optimization (a.3).

Conventions (all SI):
  - actuation displacement x > 0 = shuttle pulled INTO the comb, i.e. AWAY
    from the waveguide (Part-2 FEM verified sign), so g(V) = g0 + x(V).
  - figure-faithful mass bookkeeping (Part-2 lesson): TWO plates (modulation
    bar + comb spine, both L x Wm), column Wm x Lm, N/2 moving fingers,
    4 beams (37% effective for the resonant mode, 50% for static body load).
"""
import numpy as np

# ---------- constants ----------
EPS0 = 8.854187817e-12
E_SI = 170e9
RHO = 2330.0
LAMBDA = 1550e-9
B_OPT = 0.3e6            # absorption / index decay constant [1/m]
NEFF0, NEFFINF = 2.6, 2.4
DNEFF = NEFF0 - NEFFINF
C_AIR = 16.15            # air damping coefficient [N s/m^3] (Part 1)
Q_CLAMP = 500.0
T_SI = 200e-9            # suspended thickness, fixed in part 3a

# FEM-calibration factors from Part 2 (analytical -> expected FEM)
K_FEM_FACTOR = 0.92      # FEM spring 8-9% softer (plate compliance)
BEAM_MASS_EFF = 0.37     # modal beam-mass fraction (P2: m_eff = shuttle+0.37 m_b)
BEAM_MASS_STATIC = 0.50  # body-load effective beam fraction (P2 Q5)

X_OPT = 1.0 / (2 * B_OPT)   # optimum DC displacement = 1/(2b) = 1.667 um
V_SIDE_CAL = 222.0 / 247.0  # Rayleigh estimate -> P2-FEM-calibrated side pull-in


def tip_gap_pullin_u(gt0, h, Wc):
    """Solve g_t0 = 3u/2 + u^3/(2 h Wc) for the tip gap u at pull-in (P1 Q7)."""
    coef = 1.0 / (2 * h * Wc)
    u = np.roots([coef, 0.0, 1.5, -gt0])
    u = u[np.isreal(u)].real
    u = u[(u > 0) & (u <= gt0)]
    return float(u.min())


def finger_side_pullin_V(Lc, Wc, h, d):
    """Lateral (side) pull-in voltage of one comb finger, Rayleigh quotient.

    Cantilever finger (length Lc, width Wc) loaded by the electrostatic
    negative stiffness 2*eps0*t/h^3 per unit engaged length over the engaged
    tip segment [Lc-d, Lc]. Trial shape = static deflection under uniform
    load on that segment. t cancels.
    """
    I = T_SI * Wc**3 / 12.0
    n = 2001
    s = np.linspace(0, Lc, n)
    a = max(Lc - d, 0.0)
    # static deflection of cantilever under UDL q on [a, Lc] (unit q)
    # superposition of point loads: delta(x) integral over load positions
    # delta_P(x; xa) for unit point load at xa:
    #   x <= xa: x^2 (3 xa - x) / 6EI ; x > xa: xa^2 (3x - xa) / 6EI
    xa = np.linspace(a, Lc, n)
    w = np.zeros(n)
    for i, x in enumerate(s):
        contrib = np.where(x <= xa, x**2 * (3 * xa - x), xa**2 * (3 * x - xa))
        w[i] = np.trapezoid(contrib, xa) / (6 * E_SI * I)
    # curvature via finite differences
    wpp = np.gradient(np.gradient(w, s), s)
    U_bend = 0.5 * E_SI * I * np.trapezoid(wpp**2, s)
    mask = s >= a
    # negative-stiffness work per V^2: 0.5 * (2 eps0 t / h^3) * integral w^2
    W_es = 0.5 * (2 * EPS0 * T_SI / h**3) * np.trapezoid(w[mask]**2, s[mask])
    return float(np.sqrt(U_bend / W_es))


def evaluate(dims, f_label=''):
    """Evaluate all model quantities for a dimension dict.

    dims keys: L, g0, Lb, Wb, Lm, Wm, h, d0, Lc, Wc, N  (SI units)
    """
    L, g0, Lb, Wb = dims['L'], dims['g0'], dims['Lb'], dims['Wb']
    Lm, Wm, h, d0 = dims['Lm'], dims['Wm'], dims['h'], dims['d0']
    Lc, Wc, N = dims['Lc'], dims['Wc'], dims['N']
    t = T_SI

    r = dict(dims)
    # ---- stiffness (4 fixed-guided beams) ----
    k = 4 * E_SI * t * Wb**3 / Lb**3
    r['k'] = k
    r['k_fem_est'] = K_FEM_FACTOR * k

    # ---- masses (figure-faithful) ----
    n_mov = N // 2
    m_plates = RHO * t * (2 * L * Wm)
    m_col = RHO * t * (Wm * Lm)
    m_fing = RHO * t * (n_mov * Lc * Wc)
    m_beams = RHO * t * (4 * Lb * Wb)
    m_shuttle = m_plates + m_col + m_fing
    m_eff = m_shuttle + BEAM_MASS_EFF * m_beams
    m_static = m_shuttle + BEAM_MASS_STATIC * m_beams
    r.update(m_shuttle=m_shuttle, m_eff=m_eff, m_static=m_static,
             m_beams=m_beams)

    # ---- resonance ----
    w0 = np.sqrt(k / m_eff)
    r['f_r'] = w0 / (2 * np.pi)
    r['f_r_fem_est'] = np.sqrt(K_FEM_FACTOR) * r['f_r']

    # ---- damping / Q (Part-1 model) ----
    A_m = 2 * L * Wm + Wm * Lm + n_mov * Lc * Wc + 0.5 * 4 * Lb * Wb
    b_air = C_AIR * A_m
    b_clamp = np.sqrt(k * m_eff) / Q_CLAMP
    Q = np.sqrt(k * m_eff) / (b_air + b_clamp)
    r.update(Q=Q, b_air=b_air, b_clamp=b_clamp)

    # ---- electrostatics: ideal comb ----
    dCdx = N * EPS0 * t / h
    r['dCdx'] = dCdx

    # ---- V_DC_MAX (a.1) ----
    V_max = np.sqrt(k * h / (B_OPT * N * EPS0 * t))
    r['V_DC_MAX'] = V_max
    r['x_at_Vmax'] = X_OPT                       # = 1/(2b), geometry-independent
    r['g_op'] = g0 + X_OPT

    # ---- peak sensitivity (operate at V_DC_MAX) ----
    S_phi = (2 * np.pi * L / LAMBDA) * DNEFF * B_OPT * np.exp(-B_OPT * r['g_op'])
    r['S_phi'] = S_phi
    r['S_max'] = S_phi * V_max * dCdx / k        # rad/V at V_DC_MAX
    # closed form check: (2 pi L/lambda) dn e^{-b g0 - 1/2} / V_max
    r['S_max_closed'] = (2 * np.pi * L / LAMBDA) * DNEFF \
        * np.exp(-B_OPT * g0 - 0.5) / V_max

    # ---- axial (tip-gap) pull-in (P1 Q7 model) ----
    gt0 = Lc - d0
    r['gt0'] = gt0
    u = tip_gap_pullin_u(gt0, h, Wc)
    x_pi = gt0 - u
    V_pi = np.sqrt(k * u**3 / (N * EPS0 * Wc * t))
    r.update(u_pi=u, x_PI=x_pi, V_PI_axial=V_pi)

    # ---- side instability (Part-2 lesson): four channels, Dunkerley ----
    d_op = d0 + X_OPT                            # engagement at operating point
    kneg_f_perV2 = 2 * EPS0 * t * d_op / h**3    # per finger, both gaps
    p = 2 * (Wc + h)                             # comb pitch per mover

    # (1) individual finger bending (Rayleigh quotient, engaged-tip load)
    V_side_f = finger_side_pullin_V(Lc, Wc, h, d_op)

    # (2) shuttle rigid rocking: beam axial stiffness on lever arms, with the
    #     column bending between the two beam levels in series
    k_ax = E_SI * t * Wb / Lb
    k_tr = E_SI * t * Wb**3 / Lb**3
    s_pair = max(Lm - Wb - 2e-6, 1e-6)           # beam pair separation
    I_col = t * Wm**3 / 12.0
    k_col = 12 * E_SI * I_col / s_pair**3        # guided-guided column segment
    k_lvl = 1.0 / (1.0 / (2 * k_ax) + 1.0 / k_col)
    K_theta = 2 * k_lvl * (s_pair / 2)**2 + 4 * k_tr * (Wm / 2)**2
    sum_x2 = p**2 * n_mov * (n_mov**2 - 1) / 12.0
    V_rock = np.sqrt(K_theta / (kneg_f_perV2 * sum_x2)) if sum_x2 > 0 else np.inf

    # (3) spine lateral bending on the negative elastic foundation of the
    #     fingers: cantilever half (length ell_s) buckles when
    #     kneg/p = (1.875/ell_s)^4 EI_s
    I_spine = t * Wm**3 / 12.0
    ell_s = (L - Wm) / 2
    V_spine = np.sqrt((1.875 / ell_s)**4 * E_SI * I_spine * p
                      / (kneg_f_perV2 / 1.0)) if n_mov > 0 else np.inf

    # (4) shuttle lateral translation (beams axial)
    V_lat = np.sqrt(4 * k_ax / (n_mov * kneg_f_perV2)) if n_mov > 0 else np.inf

    V_side = 1.0 / np.sqrt(1 / V_side_f**2 + 1 / V_rock**2
                           + 1 / V_spine**2 + 1 / V_lat**2)
    V_side_cal = V_SIDE_CAL * V_side
    r.update(V_side_f=V_side_f, V_rock=V_rock, V_spine=V_spine, V_lat=V_lat,
             V_side=V_side, V_side_cal=V_side_cal, K_theta=K_theta, d_op=d_op)
    r['V_PI'] = min(V_pi, V_side_cal)

    # ---- accelerations ----
    r['a_max_dc'] = k * g0 / m_static            # quasi-static, V=0 (worst)
    r['a_max_res'] = g0 * w0**2 / Q              # resonant AC, V=0 (worst)

    # ---- secondary mode estimates (rigidity checks) ----
    cbar = np.sqrt(E_SI / (12 * RHO))            # = 2466 m/s
    ell_bar = (L - Wm) / 2
    r['f_bar_tip'] = (3.516 / (2 * np.pi)) * cbar * Wm / ell_bar**2
    r['f_finger'] = (3.516 / (2 * np.pi)) * cbar * Wc / Lc**2
    r['f_beam_cc'] = (22.37 / (2 * np.pi)) * cbar * Wb / Lb**2

    # ---- comb fill ----
    end_w = 3 * Wc
    comb_w = n_mov * p + 2 * end_w               # movers at pitch p + end fingers
    r['comb_width'] = comb_w
    r['comb_fits'] = comb_w <= L
    r['if_label'] = f_label
    return r


def sensitivity_curve(rr, V):
    """Low-frequency sensitivity S(V) [rad/V] for evaluated design rr."""
    alpha = rr['dCdx'] / (2 * rr['k'])
    g = rr['g0'] + alpha * V**2
    S_phi = (2 * np.pi * rr['L'] / LAMBDA) * DNEFF * B_OPT * np.exp(-B_OPT * g)
    return S_phi * V * rr['dCdx'] / rr['k']


# ---------- reference designs ----------
P1_DIMS = dict(L=100e-6, g0=200e-9, Lb=75e-6, Wb=5e-6, Lm=40e-6, Wm=30e-6,
               h=1e-6, d0=19e-6, Lc=20e-6, Wc=1e-6, N=40)


def report(rr, name):
    f = lambda v, s=1e6: v * s
    print(f"\n===== {name} =====")
    print(f"k          = {rr['k']:8.3f} N/m   (FEM est {rr['k_fem_est']:.2f})")
    print(f"m_shuttle  = {rr['m_shuttle']*1e12:8.3f} pg, m_eff = {rr['m_eff']*1e12:.3f} pg")
    print(f"f_r        = {rr['f_r']/1e3:8.1f} kHz  (FEM est {rr['f_r_fem_est']/1e3:.1f})")
    print(f"Q          = {rr['Q']:8.1f}")
    print(f"dC/dx      = {rr['dCdx']:.3e} F/m")
    print(f"V_DC_MAX   = {rr['V_DC_MAX']:8.1f} V  (x at peak = {X_OPT*1e6:.3f} um)")
    print(f"V_PI axial = {rr['V_PI_axial']:8.1f} V  (x_PI = {rr['x_PI']*1e6:.3f} um, u = {rr['u_pi']*1e6:.3f} um)")
    print(f"V_side_f   = {rr['V_side_f']:8.1f} V, V_rock = {rr['V_rock']:.1f} V, "
          f"V_spine = {rr['V_spine']:.1f} V, V_lat = {rr['V_lat']:.1f} V")
    print(f"V_side     = {rr['V_side']:8.1f} V (calibrated {rr['V_side_cal']:.1f} V)")
    print(f"V_PI (min) = {rr['V_PI']:8.1f} V")
    print(f"S_max      = {rr['S_max']*1e3:8.2f} mrad/V (closed form {rr['S_max_closed']*1e3:.2f})")
    print(f"a_max DC   = {rr['a_max_dc']:.3e} m/s^2 ; a_max(res) = {rr['a_max_res']:.3e} m/s^2")
    print(f"f_bar_tip  = {rr['f_bar_tip']/1e6:6.2f} MHz, f_finger = {rr['f_finger']/1e6:.2f} MHz, f_beam = {rr['f_beam_cc']/1e6:.2f} MHz")
    print(f"comb width = {rr['comb_width']*1e6:.1f} um (fits: {rr['comb_fits']})")


if __name__ == '__main__':
    # ----- calibration vs Part 1 / Part 2 -----
    rr = evaluate(P1_DIMS)
    report(rr, 'Part 1/2 design (calibration)')
    print("\nCalibration targets: k=40.30 (P1), f_FEM=496 kHz, "
          "V_PI_axial=347 V (P1), V_side(FEM P2)=222 V")
    # sensitivity at P1's quoted bias for sanity
    S10 = sensitivity_curve(rr, np.array([10.0]))[0]
    print(f"S(10 V) = {S10*1e3:.3f} mrad/V (P1 quoted 0.40 mrad/V)")
