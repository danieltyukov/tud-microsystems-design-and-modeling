"""ET4260 Part 3a -- freeze the final design, export summary CSV + figures.

Final design C1 (from p3a_optimize grid search + rounding to clean values).
"""
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import p3a_model as M

BASE = '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p3'

FINAL = dict(L=100e-6, g0=200e-9, Lb=40e-6, Wb=1.85e-6, Lm=30e-6, Wm=10e-6,
             h=0.9e-6, d0=1.0e-6, Lc=5.5e-6, Wc=0.7e-6, N=58)

P2_FEM_VPI = 222.0   # Part-2 FEM pull-in of the old design (rocking)


def main():
    old = M.evaluate(M.P1_DIMS)
    new = M.evaluate(FINAL)
    M.report(old, 'Part 1/2 design')
    M.report(new, 'Part 3a final design')

    # ---------- summary CSV ----------
    rows = [('quantity', 'unit', 'Part1/2', 'Part3a')]
    take = [('k', 'N/m', 1), ('m_shuttle', 'pg', 1e12), ('m_eff', 'pg', 1e12),
            ('f_r', 'kHz', 1e-3), ('f_r_fem_est', 'kHz', 1e-3), ('Q', '-', 1),
            ('dCdx', 'F/m', 1), ('V_DC_MAX', 'V', 1), ('V_PI_axial', 'V', 1),
            ('V_side_cal', 'V', 1), ('V_PI', 'V', 1), ('x_PI', 'um', 1e6),
            ('u_pi', 'um', 1e6), ('S_max', 'mrad/V', 1e3),
            ('a_max_dc', 'm/s^2', 1), ('a_max_res', 'm/s^2', 1),
            ('f_bar_tip', 'MHz', 1e-6), ('f_finger', 'MHz', 1e-6),
            ('f_beam_cc', 'MHz', 1e-6)]
    for key, unit, sc in take:
        rows.append((key, unit, f"{old[key]*sc:.6g}", f"{new[key]*sc:.6g}"))
    with open(f'{BASE}/data/p3a_design_summary.csv', 'w', newline='') as fh:
        csv.writer(fh).writerows(rows)

    # dimension table CSV (final-value column of the assignment table)
    dim_rows = [('symbol', 'Part1/2 [um]', 'Part3a final [um]')]
    for key in ['L', 'g0', 'Lb', 'Wb', 'Lm', 'Wm', 'h', 'd0', 'Lc', 'Wc', 'N']:
        sc = 1 if key == 'N' else 1e6
        dim_rows.append((key, f"{M.P1_DIMS[key]*sc:g}", f"{FINAL[key]*sc:g}"))
    dim_rows.append(('3Wc', f"{3*M.P1_DIMS['Wc']*1e6:g}", f"{3*FINAL['Wc']*1e6:g}"))
    dim_rows.append(('t', '0.2', '0.2'))
    with open(f'{BASE}/data/p3a_dimensions.csv', 'w', newline='') as fh:
        csv.writer(fh).writerows(dim_rows)

    # ---------- figure: S(V) for both designs ----------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    V = np.linspace(0, 900, 2000)

    for ax, rr, vpi_fem, name, col in [
            (ax1, old, P2_FEM_VPI, 'Part 1/2 design', 'tab:gray'),
            (ax2, new, None, 'Part 3a design', 'tab:blue')]:
        S = M.sensitivity_curve(rr, V) * 1e3
        vpi = vpi_fem if vpi_fem else rr['V_PI']
        ok = V <= vpi
        ax.plot(V[ok], S[ok], color=col, lw=2, label='stable (below pull-in)')
        ax.plot(V[~ok], S[~ok], color=col, lw=1.2, ls=':',
                label='beyond pull-in (unreachable)')
        ax.axvline(rr['V_DC_MAX'], color='tab:red', ls='--', lw=1,
                   label=f"$V_{{DC\\_MAX}}$ = {rr['V_DC_MAX']:.0f} V")
        ax.axvline(vpi, color='k', ls='-.', lw=1,
                   label=f"$V_{{PI}}$ = {vpi:.0f} V")
        iop = np.argmin(np.abs(V - min(rr['V_DC_MAX'], vpi)))
        ax.plot(V[iop], S[iop], 'o', color='tab:red', ms=7, zorder=5)
        ax.set_xlabel('$V_{DC}$ (V)')
        ax.set_ylabel('$|\\Delta\\phi/V_{AC}|_{f\\to 0}$ (mrad/V)')
        ax.set_title(name)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc='upper right' if rr is old else 'upper left')
        ax.set_xlim(0, 900)
    ax1.set_ylim(0, 40)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(f'{BASE}/figures/p3a_sensitivity_vs_V.{ext}', dpi=200)
    plt.close(fig)

    # ---------- figure: travel condition (a.2) ----------
    # S(x) = Sphi(x) * dCdx * V(x) / k with V(x) = sqrt(2 k x / dCdx)
    fig, ax = plt.subplots(figsize=(7, 4.2))
    x = np.linspace(0, 4e-6, 800)
    for rr, col, name in [(old, 'tab:gray', 'Part 1/2'),
                          (new, 'tab:blue', 'Part 3a')]:
        Vx = np.sqrt(2 * rr['k'] * x / rr['dCdx'])
        Sx = (2*np.pi*rr['L']/M.LAMBDA)*M.DNEFF*M.B_OPT \
            * np.exp(-M.B_OPT*(rr['g0']+x)) * rr['dCdx'] * Vx / rr['k'] * 1e3
        # reachable travel: new design = axial tip-gap fold; old design =
        # displacement at the Part-2 FEM (rocking) pull-in voltage
        xpi = rr['x_PI'] if rr is new else \
            old['dCdx'] / (2 * old['k']) * P2_FEM_VPI**2
        ok = x <= xpi
        ax.plot(x[ok]*1e6, Sx[ok], color=col, lw=2, label=f'{name}: reachable')
        ax.plot(x[~ok]*1e6, Sx[~ok], color=col, lw=1.2, ls=':')
        ax.axvline(xpi*1e6, color=col, ls='-.', lw=1)
    ax.axvline(M.X_OPT*1e6, color='tab:red', ls='--', lw=1.2,
               label='$x_{opt}=1/(2b)=1.67\\,\\mu$m')
    ax.set_xlabel('DC displacement $x$ ($\\mu$m)')
    ax.set_ylabel('$|\\Delta\\phi/V_{AC}|_{f\\to 0}$ (mrad/V)')
    ax.set_title('Sensitivity vs DC displacement: the travel condition')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(f'{BASE}/figures/p3a_travel_condition.{ext}', dpi=200)
    plt.close(fig)
    print('\nwrote data/p3a_design_summary.csv, data/p3a_dimensions.csv,'
          ' figures/p3a_sensitivity_vs_V.*, figures/p3a_travel_condition.*')


if __name__ == '__main__':
    main()
