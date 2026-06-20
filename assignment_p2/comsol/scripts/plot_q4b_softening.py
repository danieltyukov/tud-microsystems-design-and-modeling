"""Regenerate the Q4b electrostatic-softening figure to match the reframed text:
(f/f0)^2 plotted against V^2, with the linear f^2 ~ (1 - V^2/V_inst^2) fits
extrapolated to zero -> the intercept IS the instability voltage. Mesh-symmetry
independent (eigenvalue method), unlike the parametric sweep.

Reads data/q4b_eig_vs_bias.csv (or *_fine.csv when re-run with finer mesh +
electrostatic thickness=1). Pure matplotlib, no COMSOL needed.

  python3 plot_q4b_softening.py [eig_csv]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p2'
csv = sys.argv[1] if len(sys.argv) > 1 else f'{BASE}/data/q4b_eig_vs_bias.csv'
d = np.genfromtxt(csv, delimiter=',', names=True)
V = d['Vd_V']
modes = [('f1_Hz', 'mode 1: shuttle (axial)', 'C0'),
         ('f2_Hz', 'mode 2: rocking', 'C3'),
         ('f3_Hz', 'mode 3: finger band', 'C2')]

fig, ax = plt.subplots(figsize=(7.4, 4.6))
V2 = V**2
for key, label, c in modes:
    f = d[key]
    y = (f / f[0])**2
    ax.plot(V2 / 1e3, y, 'o-', color=c, lw=1.6, ms=5, label=label)
    # linear fit of f^2 vs V^2 over the points nearest the fold (most accurate
    # extrapolation): last 4 bias points
    m = V >= V[-4]
    sl, it = np.polyfit(V2[m], f[m]**2, 1)
    if sl < 0:
        Vpi2 = -it / sl
        Vpi = np.sqrt(Vpi2)
        xs = np.linspace(V2.min(), Vpi2, 50)
        ax.plot(xs / 1e3, (sl * xs + it) / f[0]**2, '--', color=c, lw=1.0,
                alpha=0.7)
        ax.plot(Vpi2 / 1e3, 0, '*', color=c, ms=15,
                label=f'   $\\to V_{{\\rm inst}}\\approx{Vpi:.0f}$ V')

# overlay finer-mesh + thickness=1 rocking points to show mesh-convergence
import os
fine_csv = f'{BASE}/data/q4b_eig_vs_bias_fine.csv'
if os.path.exists(fine_csv):
    df = np.genfromtxt(fine_csv, delimiter=',', names=True)
    yf = (df['f2_Hz'] / df['f2_Hz'][0])**2
    ax.plot(df['Vd_V']**2 / 1e3, yf, 'kx', ms=8, mew=1.8,
            label='rocking, fine mesh (51k, $d{=}1$)')

ax.axhline(0, color='k', lw=0.8)
ax.set_xlabel(r'$V_d^{\,2}$  (kV$^2$)')
ax.set_ylabel(r'$(f/f_0)^2$   (normalized modal stiffness)')
ax.set_title('Q4b: electrostatic softening is linear in $V^2$ '
             r'($f^2\propto V_{\rm inst}^2-V^2$)')
ax.set_ylim(-0.05, 1.05)
ax.grid(alpha=0.3)
ax.legend(loc='upper right', fontsize=9, ncol=2)
fig.tight_layout()
out = f'{BASE}/report/figures/q4b_eig_softening.png'
fig.savefig(out, dpi=150)
fig.savefig(f'{BASE}/figures/q4b_eig_softening.png', dpi=150)
print('saved', out)
