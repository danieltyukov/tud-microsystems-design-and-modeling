"""Generate report figures from the exported CSV data (matplotlib).
Re-runnable any time; skips plots whose data files are missing.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p2'
DATA = f'{BASE}/data'
FIGS = f'{BASE}/figures'

plt.rcParams.update({'font.size': 11, 'figure.dpi': 140,
                     'axes.grid': True, 'grid.alpha': 0.3})


def have(name):
    return os.path.exists(f'{DATA}/{name}')


# ---- Q3 step response ----
if have('q3_step_response.csv'):
    d = np.genfromtxt(f'{DATA}/q3_step_response.csv', delimiter=',', names=True)
    t = d['t_s'] * 1e6
    x = d['uy_mass_m'] * 1e9
    xf = x[-1]
    Q, f1 = 83.0, 496.26e3
    zw = 2 * np.pi * f1 / (2 * Q)
    ts_ana = np.log(20) / zw * 1e6
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(t, x, lw=0.7, label='FEM $u_y(t)$ (mass)')
    ax.plot(t, xf * (1 - np.exp(-zw * t * 1e-6)), 'r--', lw=1.2,
            label='envelope $x_f(1-e^{-\\zeta\\omega_1 t})$')
    ax.plot(t, xf * (1 + np.exp(-zw * t * 1e-6)), 'r--', lw=1.2)
    ax.axhline(xf * 1.05, color='gray', ls=':', lw=0.8)
    ax.axhline(xf * 0.95, color='gray', ls=':', lw=0.8)
    ax.axvline(158.4, color='g', ls='-.', lw=1,
               label='$t_s$(FEM, 5%) = 158.4 µs')
    ax.set_xlabel('t (µs)')
    ax.set_ylabel('$u_y$ (nm)')
    ax.set_title(f'Q3: step response, Rayleigh damping Q=83 '
                 f'(analytic $t_s$ = {ts_ana:.1f} µs)')
    ax.legend(loc='lower right')
    fig.tight_layout()
    fig.savefig(f'{FIGS}/q3_step_response.png')
    plt.close(fig)
    print('q3 plot done')

# ---- Q4a displacement vs voltage ----
if have('q4a_displacement_vs_voltage.csv'):
    d = np.genfromtxt(f'{DATA}/q4a_displacement_vs_voltage.csv',
                      delimiter=',', names=True)
    V = d['Vd_V']
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(V, -d['uy_mass_m'] * 1e9, 'o-', label='FEM (COMSOL)')
    ax.plot(V, d['x_P1_analytical_mdownwardneg'] * -1e9 if
            'x_P1_analytical_mdownwardneg' in d.dtype.names else
            0.5 * V**2 * 7.083e-11 / 40.2963 * 1e9, 's--',
            label='Part 1 analytical $x=\\frac{V^2}{2k}\\frac{N\\epsilon_0 t}{h}$')
    ax.set_xlabel('$V_d$ (V)')
    ax.set_ylabel('|displacement| (nm)')
    ax.set_title('Q4a: static displacement vs comb-drive voltage')
    ax.legend()
    fig.tight_layout()
    fig.savefig(f'{FIGS}/q4a_displacement_vs_voltage.png')
    plt.close(fig)
    print('q4a plot done')

# ---- Q4b pull-in curve ----
if have('q4b_pullin_curve.csv'):
    d = np.genfromtxt(f'{DATA}/q4b_pullin_curve.csv', delimiter=',', names=True)
    x = np.atleast_1d(-d['uy_mass_m']) * 1e9
    V = np.abs(np.atleast_1d(d['Vdc_V']))
    ipk = int(np.argmax(V))
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x, V, 'o-', label='FEM equilibrium curve V(x)')
    ax.plot(x[ipk], V[ipk], 'r*', ms=16,
            label=f'pull-in: {V[ipk]:.1f} V at {x[ipk]:.0f} nm')
    ax.axvline(404, color='gray', ls=':', label='P1: $x_{PI}$=404 nm, $V_{PI}$=347 V')
    ax.axhline(347.1, color='gray', ls=':')
    ax.set_xlabel('|displacement| (nm)')
    ax.set_ylabel('V (V)')
    ax.set_title('Q4b: displacement-controlled pull-in curve')
    ax.legend()
    fig.tight_layout()
    fig.savefig(f'{FIGS}/q4b_pullin_curve.png')
    plt.close(fig)
    print('q4b plot done')

# ---- Q4b voltage-sweep stable branch (continuation to divergence) ----
if have('q4b_vsweep_curve.csv'):
    d = np.genfromtxt(f'{DATA}/q4b_vsweep_curve.csv', delimiter=',', names=True)
    V = np.atleast_1d(d['Vd_V'])
    x = np.atleast_1d(-d['uy_mass_m']) * 1e9
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(V, x, '-', lw=1.5, label='FEM stable branch (continuation)')
    ax.plot(V[-1], x[-1], 'rx', ms=12, mew=2.5,
            label=f'last converged: {V[-1]:.0f} V, {x[-1]:.0f} nm')
    ax.axvline(347.1, color='gray', ls=':', label='P1 analytical $V_{PI}$ = 347 V')
    ax.set_xlabel('$V_d$ (V)')
    ax.set_ylabel('|displacement| (nm)')
    ax.set_title('Q4b: stable equilibrium branch up to pull-in')
    ax.legend()
    fig.tight_layout()
    fig.savefig(f'{FIGS}/q4b_vsweep_branch.png')
    plt.close(fig)
    print('q4b vsweep plot done')

# ---- Q4b eigenfrequency softening vs bias (pull-in mechanism) ----
if have('q4b_eig_vs_bias.csv'):
    d = np.genfromtxt(f'{DATA}/q4b_eig_vs_bias.csv', delimiter=',', names=True)
    V = d['Vd_V']
    fig, ax = plt.subplots(figsize=(8, 5))
    for col, lab, c in [('f1_Hz', 'mode 1: shuttle (axial)', 'tab:blue'),
                        ('f2_Hz', 'mode 2: rocking', 'tab:red'),
                        ('f3_Hz', 'mode 3: finger band', 'tab:green')]:
        ax.plot(V, (d[col] / d[col][0])**2, 'o-', color=c, label=lab)
    # f2^2 linear extrapolation to zero from the last two points
    f2n = (d['f2_Hz'] / d['f2_Hz'][0])**2
    sl = (f2n[-1] - f2n[-2]) / (V[-1] - V[-2])
    Vzero = V[-1] - f2n[-1] / sl
    ax.plot([V[-1], Vzero], [f2n[-1], 0], 'r:', lw=1.5)
    ax.plot(Vzero, 0, 'r*', ms=14,
            label=f'$f_2^2 \\to 0$ at {Vzero:.0f} V (FEM div. 222 V)')
    ax.axhline(0, color='k', lw=0.8)
    ax.set_xlabel('$V_d$ (V)')
    ax.set_ylabel('$(f/f_0)^2$  (normalized stiffness)')
    ax.set_title('Q4b mechanism: electrostatic softening of the lowest modes')
    ax.legend()
    fig.tight_layout()
    fig.savefig(f'{FIGS}/q4b_eig_softening.png')
    plt.close(fig)
    print('q4b eig plot done')

# ---- Q5 ----
if have('q5_max_acceleration.csv'):
    d = np.genfromtxt(f'{DATA}/q5_max_acceleration.csv', delimiter=',', names=True)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(d['Vd_V'], d['a_max_m_s2'] / 1e6, 'o-')
    # stagger labels: the 0/1/10 V points sit nearly on top of each other
    offsets = [(8, -14), (8, 4), (8, 14), (-12, 10)]
    for (vv, aa), off in zip(zip(d['Vd_V'], d['a_max_m_s2']), offsets):
        ax.annotate(f'{aa/9.81:,.0f} g', (vv, aa / 1e6),
                    textcoords='offset points', xytext=off, fontsize=9)
    ax.set_xlabel('$V_d$ (V)')
    ax.set_ylabel('$a_{max}$ (10$^6$ m/s$^2$)')
    ax.set_title('Q5: max acceleration before mass-waveguide contact')
    fig.tight_layout()
    fig.savefig(f'{FIGS}/q5_max_acceleration.png')
    plt.close(fig)
    print('q5 plot done')

# ---- Q6 voltage steps + comparison with Q3 ----
if have('q6_vstep_50V.csv') and have('q6_vstep_1V.csv'):
    d1 = np.genfromtxt(f'{DATA}/q6_vstep_1V.csv', delimiter=',', names=True)
    d50 = np.genfromtxt(f'{DATA}/q6_vstep_50V.csv', delimiter=',', names=True)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
    axs[0].plot(d1['t_s'] * 1e6, d1['uy_mass_m'] * 1e12, lw=0.7)
    axs[0].set_xlabel('t (µs)'); axs[0].set_ylabel('$u_y$ (pm)')
    axs[0].set_title('Q6: $V_d$ step 0 → 1 V')
    axs[1].plot(d50['t_s'] * 1e6, d50['uy_mass_m'] * 1e9, lw=0.7)
    axs[1].set_xlabel('t (µs)'); axs[1].set_ylabel('$u_y$ (nm)')
    axs[1].set_title('Q6: $V_d$ step 0 → 50 V')
    fig.tight_layout()
    fig.savefig(f'{FIGS}/q6_vstep_responses.png')
    plt.close(fig)

    # normalized overlay with Q3
    if have('q3_step_response.csv'):
        d3 = np.genfromtxt(f'{DATA}/q3_step_response.csv', delimiter=',', names=True)
        fig, ax = plt.subplots(figsize=(9, 5))
        for d, lab in [(d3, 'Q3 force step'), (d1, 'Q6 1 V step'),
                       (d50, 'Q6 50 V step')]:
            tt = d['t_s'] * 1e6
            uu = d['uy_mass_m'] / d['uy_mass_m'][-1]
            ax.plot(tt, uu, lw=0.7, label=lab)
        ax.set_xlabel('t (µs)')
        ax.set_ylabel('$u_y / u_{y,final}$')
        ax.set_title('Q3 vs Q6: normalized step responses')
        ax.legend()
        ax.set_xlim(0, 250)
        fig.tight_layout()
        fig.savefig(f'{FIGS}/q6_vs_q3_normalized.png')
        plt.close(fig)
    print('q6 plots done')

print('ALL PLOTS DONE')
