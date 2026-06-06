"""Regenerate all matplotlib figures from the exported FEM CSVs.

Run: ~/workspace/comsol-mcp/.venv/bin/python 04_make_plots.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

BASE = Path('/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p3')
DATA, FIGS = BASE / 'data', BASE / 'figures'

K_FEM = None  # filled from spring CSV
V_OP = 627.0
X_OPT = 1.0 / (2 * 0.3e6)


def save(fig, name):
    for ext in ('png', 'pdf'):
        fig.savefig(FIGS / f'{name}.{ext}', dpi=200)
    plt.close(fig)
    print('wrote', name)


# ---------------- step response (a.4 v) ----------------
f = DATA / 'p3a_q4v_step_response.csv'
if f.exists():
    d = np.genfromtxt(f, delimiter=',', names=True)
    t, uy = d['t_s'] * 1e6, d['uy_mass_m'] * 1e9
    xf = uy[-1]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 4))
    for a in (ax, ax2):
        a.plot(t, uy, lw=0.8, color='tab:blue')
        a.axhline(xf, color='k', lw=0.6)
        a.axhline(0.95 * xf, color='gray', ls='--', lw=0.6)
        a.axhline(1.05 * xf, color='gray', ls='--', lw=0.6)
        a.set_xlabel('t ($\\mu$s)')
        a.set_ylabel('$u_y$ mass (nm)')
        a.grid(alpha=0.3)
    out = np.abs(uy - xf) > 0.05 * abs(xf)
    ts = t[np.max(np.nonzero(out))] if out.any() else 0.0
    ax.axvline(ts, color='tab:red', ls='-.', lw=1, label=f'$t_s$(5%) = {ts:.1f} $\\mu$s')
    ax.legend(fontsize=8)
    ax.set_title('Force-step response (1 $\\mu$N), full window')
    ax2.set_xlim(0, 10)
    i10 = np.argmax(np.abs(uy) >= 0.1 * abs(xf))
    i90 = np.argmax(np.abs(uy) >= 0.9 * abs(xf))
    ax2.plot(t[[i10, i90]], uy[[i10, i90]], 'o', color='tab:red', ms=5,
             label=f'10-90% rise = {(t[i90]-t[i10])*1e3:.0f} ns')
    ax2.legend(fontsize=8)
    ax2.set_title('First microseconds (rise)')
    fig.tight_layout()
    save(fig, 'p3a_q4v_step_response')

# ---------------- x(V) sweep (a.4 iii) + FEM-true sensitivity ----------------
f = DATA / 'p3a_q4iii_vsweep.csv'
if f.exists():
    d = np.genfromtxt(f, delimiter=',', names=True)
    V, uy = d['Vd_V'], np.abs(d['uy_mass_m'])
    ks = np.genfromtxt(DATA / 'p3a_q4i_spring_constant.csv', delimiter=',',
                       names=True)
    k_fem = float(ks['k_mass'][-1])
    dCdx = 58 * 8.854187817e-12 * 200e-9 / 0.9e-6
    Vth = np.linspace(0, V[-1], 400)
    xth = dCdx * Vth**2 / (2 * 13.45) * 1e6     # analytic k
    xth_kfem = dCdx * Vth**2 / (2 * k_fem) * 1e6

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    ax.plot(V, uy * 1e6, 'o-', ms=2.5, lw=1, label='FEM (continuation sweep)')
    ax.plot(Vth, xth, 'k--', lw=1, label='analytic, $k$ = 13.45 N/m')
    ax.plot(Vth, xth_kfem, ':', color='tab:red', lw=1.2,
            label=f'analytic with $k_{{FEM}}$ = {k_fem:.2f} N/m')
    ax.axvline(V_OP, color='tab:red', ls='--', lw=0.8)
    ax.axhline(X_OPT * 1e6, color='tab:green', ls='--', lw=0.8,
               label='$x_{opt}=1/(2b)=1.667\\,\\mu$m')
    ax.set_xlabel('$V_{DC}$ (V)')
    ax.set_ylabel('|$u_y$| ($\\mu$m)')
    ax.set_title(f'Displacement vs bias (converged through {V[-1]:.0f} V)')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    # force-implied effective dC/dx: F = k_fem*x = 0.5 V^2 (dC/dx)_eff
    with np.errstate(divide='ignore', invalid='ignore'):
        dCdx_eff = 2 * k_fem * uy / V**2
    ax2.plot(V[1:], dCdx_eff[1:] * 1e12, 'o-', ms=2.5, lw=1, color='tab:purple',
             label='FEM force-implied (secant)')
    ax2.axhline(dCdx * 1e12, color='k', ls='--', lw=1,
                label='ideal comb $N\\varepsilon_0 t/h$')
    ax2.set_xlabel('$V_{DC}$ (V)')
    ax2.set_ylabel('$(dC/dx)_{eff}$ (pF/m)')
    ax2.set_title('Effective comb constant: shallow-engagement deficit')
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3)
    fig.tight_layout()
    save(fig, 'p3a_q4iii_vsweep')

    # FEM-true low-frequency sensitivity |H|(V) = S_phi(g0+x) * dx/dV
    LAM, DN, B = 1550e-9, 0.2, 0.3e6
    Sphi = (2 * np.pi * 100e-6 / LAM) * DN * B * np.exp(-B * (200e-9 + uy))
    dxdV = np.gradient(uy, V)
    H_fem = Sphi * dxdV * 1e3                  # mrad/V
    # analytic curve (ideal comb, analytic k)
    alpha = dCdx / (2 * 13.45)
    xa = alpha * Vth**2
    Ha = (2 * np.pi * 100e-6 / LAM) * DN * B * np.exp(-B * (200e-9 + xa)) \
        * 2 * alpha * Vth * 1e3
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot(Vth, Ha, 'k--', lw=1.2, label='analytic (ideal comb, $k_{ana}$)')
    ax.plot(V[1:], H_fem[1:], 'o-', ms=3, lw=1, color='tab:blue',
            label='FEM-derived: $S_\\phi(g_0+x)\\,dx/dV$')
    ax.axvline(V_OP, color='tab:red', ls='--', lw=0.8,
               label='design $V_{DC\\_MAX}$ = 627 V')
    imax = np.argmax(H_fem[1:]) + 1
    ax.plot(V[imax], H_fem[imax], '*', color='tab:red', ms=14, mec='k',
            label=f'FEM peak so far: {H_fem[imax]:.1f} mrad/V @ {V[imax]:.0f} V')
    ax.set_xlabel('$V_{DC}$ (V)')
    ax.set_ylabel('$|\\Delta\\phi/V_{AC}|_{f\\to0}$ (mrad/V)')
    ax.set_title('Low-frequency modulation sensitivity: analytic vs FEM-derived')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    save(fig, 'p3a_sensitivity_fem')
    print(f'FEM |H| at 627 V: {H_fem[np.argmin(np.abs(V-627))]:.1f} mrad/V; '
          f'at {V[-1]:.0f} V: {H_fem[-1]:.1f} mrad/V')

# ---------------- eigenfrequency vs bias ----------------
f = DATA / 'p3a_eig_vs_bias.csv'
if f.exists():
    d = np.genfromtxt(f, delimiter=',', names=True)
    V = d['Vd_V']
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    for i in range(1, 6):
        key = f'f{i}_Hz'
        if key in d.dtype.names:
            ax.plot(V, d[key] / 1e3, 'o-', ms=3, lw=1, label=f'mode {i}')
    ax.axvline(V_OP, color='tab:red', ls='--', lw=1,
               label=f'$V_{{DC\\_MAX}}$ = {V_OP:.0f} V')
    ax.set_xlabel('$V_{DC}$ (V)')
    ax.set_ylabel('eigenfrequency (kHz)')
    ax.set_title('Prestressed eigenfrequencies vs bias (electrostatic softening)')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    save(fig, 'p3a_eig_softening')

# ---------------- part b: modes vs thickness ----------------
f = DATA / 'p3b_modes_vs_t.csv'
if f.exists():
    d = np.genfromtxt(f, delimiter=',', names=True)
    d = d[d['f1_Hz'] > 1e3]          # drop degenerate eigensolves
    t = d['t_um']
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    for i in range(1, 7):
        fk, ok = d[f'f{i}_Hz'] / 1e3, d[f'oop{i}']
        ip = ok < 0.5
        ax.plot(t, fk, '-', color='lightgray', lw=0.7, zorder=1)
        ax.scatter(t[ip], fk[ip], c='tab:blue', s=18, zorder=2,
                   label='in-plane' if i == 1 else None)
        ax.scatter(t[~ip], fk[~ip], c='tab:orange', s=18, zorder=2,
                   label='out-of-plane' if i == 1 else None)
    ax.axvline(5.0, color='k', ls='--', lw=0.8, label='$t=W_b=5\\,\\mu$m')
    ax.axhline(500, color='tab:red', ls=':', lw=0.8, label='500 kHz spec')
    ax.set_xlabel('thickness $t$ ($\\mu$m)')
    ax.set_ylabel('eigenfrequency (kHz)')
    ax.set_yscale('log')
    ax.set_title('3D eigenfrequencies vs device thickness (Part-2 dims, no fingers)')
    ax.grid(alpha=0.3, which='both')
    ax.legend(fontsize=8)
    fig.tight_layout()
    save(fig, 'p3b_modes_vs_t')

print('all plots done')
