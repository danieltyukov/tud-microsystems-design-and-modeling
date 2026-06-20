"""Regenerate the Q6 figures consistent with the corrected text. Uses a directly
resolved 1 V trace (q6_vstep_1V_fine.csv) if present; otherwise shows the
physical 1 V as the well-resolved 50 V trace scaled by 1/2500 (linearity, since
the device is far below pull-in) -- never the under-resolved first-submission
1 V trace, whose 75 us decay is a numerical-dissipation artefact.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p2'
D, F = f'{BASE}/data', f'{BASE}/report/figures'


def load(name):
    p = f'{D}/{name}'
    return np.genfromtxt(p, delimiter=',', names=True) if os.path.exists(p) else None


def settle(t, u):
    xf = u[-1]; out = np.abs(u - xf) > 0.05 * abs(xf)
    return (t[np.max(np.nonzero(out))] if out.any() else 0.0), xf


d50 = load('q6_vstep_50V_fine.csv') or load('q6_vstep_50V.csv')
d3 = load('q3_step_response.csv')
d1_real = load('q6_vstep_1V_fine.csv')   # only the re-run direct trace
t50, u50 = d50['t_s'], d50['uy_mass_m']

# 1 V curve: real direct trace if available, else 50 V / 2500 (linearity)
if d1_real is not None:
    t1, u1, lab1 = d1_real['t_s'], d1_real['uy_mass_m'], '0 $\\to$ 1 V (direct)'
else:
    t1, u1, lab1 = t50, u50 / 2500.0, '0 $\\to$ 1 V (= 50 V/2500, linearity)'

# side-by-side responses with settling markers
fig, axs = plt.subplots(1, 2, figsize=(12, 4.4))
for ax, t, u, sc, unit, lab in [(axs[0], t1, u1, 1e12, 'pm', lab1),
                                (axs[1], t50, u50, 1e9, 'nm', '0 $\\to$ 50 V')]:
    ts, xf = settle(t, u)
    ax.plot(t * 1e6, u * sc, lw=0.8, color='C0')
    ax.axhline(xf * sc, ls=':', color='0.4', lw=0.8)
    for k in (1.05, 0.95):
        ax.axhline(xf * sc * k, ls='--', color='0.75', lw=0.7)
    ax.axvline(ts * 1e6, ls=':', color='C3', lw=1.3)
    ax.text(ts * 1e6 + 4, ax.get_ylim()[1] * 0.85, f'$t_s$={ts*1e6:.0f} µs',
            color='C3', fontsize=9)
    ax.set_xlabel('t (µs)'); ax.set_ylabel(f'$u_y$ ({unit})'); ax.set_title(f'Q6: $V_d$ step {lab}')
fig.tight_layout(); fig.savefig(f'{F}/q6_vstep_responses.png', dpi=150); plt.close(fig)
print('saved q6_vstep_responses.png')

# normalized overlay: Q3 force step + 50 V voltage step (1 V coincides with 50 V)
fig, ax = plt.subplots(figsize=(9, 5))
for t, u, lab, c in [(d3['t_s'], d3['uy_mass_m'], 'Q3 force step', 'C2'),
                     (t50, u50, 'Q6 50 V step ($\\equiv$ 1 V, normalized)', 'C3')]:
    ax.plot(t * 1e6, u / u[-1], lw=0.9, label=lab, color=c)
ax.axhline(1, ls=':', color='0.5', lw=0.8)
for k in (1.05, 0.95):
    ax.axhline(k, ls='--', color='0.8', lw=0.6)
ax.set_xlabel('t (µs)'); ax.set_ylabel('$u_y / u_{y,\\rm final}$')
ax.set_title('Q3 vs Q6: normalized step responses (same $f_1$ ringing, same settling)')
ax.legend(); ax.set_xlim(0, 250)
fig.tight_layout(); fig.savefig(f'{F}/q6_vs_q3_normalized.png', dpi=150); plt.close(fig)
print('saved q6_vs_q3_normalized.png')
ts1, xf1 = settle(t1, u1)
print(f'1V shown: final={xf1*1e12:.3f} pm  settling={ts1*1e6:.1f} us  '
      f'({"direct" if d1_real is not None else "linearity"})')
