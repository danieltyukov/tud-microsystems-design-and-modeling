"""Q4b RE-RUN addressing reviewer feedback (final submission).

Reviewer (Xiaoxi Zhao) points on Q4b:
  - a parametric sweep is NOT a proper way to get pull-in (last-converged V is
    mesh- and step-dependent);
  - how can a *symmetric* structure show lateral pull-in? (suspect mesh asym.);
  - try electrostatic thickness = 1 -> same solution, far better convergence.

This script does three things, all with es thickness d = 1 (per the tip):
  (A) prestressed eigenfrequency softening on a FINER mesh (the proper, mesh-
      symmetry-independent instability test). f_i^2 -> 0 gives each mode's
      instability voltage; mode shape identifies axial pull-in vs side (rocking)
      instability. -> data/q4b_eig_vs_bias_fine.csv
  (B) the same softening on TWO mesh densities to show the eigenvalue result is
      mesh-CONVERGED (unlike the sweep).
  (C) a stationary voltage continuation sweep on the finer mesh, to show its
      last-converged voltage MOVES vs the coarse-mesh sweep -> demonstrates the
      reviewer's point that the sweep value is not the pull-in.
      -> data/q4b_vsweep_fine.csv

Run (license restored):
  LD_LIBRARY_PATH=/usr/local/comsol64/multiphysics/lib/glnxa64 \
  ~/workspace/comsol-mcp/.venv/bin/python 14_q4b_rerun.py
"""
import sys
sys.path.insert(0, '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p2/comsol/scripts')
import mph
import numpy as np
from jpype import JInt
from model_lib import MODELS, DATA, FIGS

MSC_FINE = 0.55      # finer than the 1.0 used in the first submission
client = mph.start(cores=6)
model = client.load(f'{MODELS}/phase_shifter_em.mph')
jm = model.java
comp = jm.component('comp1')
es = comp.physics('es')

# ---- reviewer tip: electrostatic out-of-plane thickness = 1 ----------------
# The Maxwell stress (Pa) driving the structure is independent of the ES
# thickness d; d only scales the charge/energy READOUT. Setting d = 1 leaves
# the displacement/force/instability unchanged but normalises the tiny
# (~2e-7-scaled) charge integrals, greatly improving Newton conditioning.
es.prop('d').set('d', '1')
print('es thickness set to 1 (was t_si); displacement solution unchanged')


def run_softening(msc, tag_csv):
    model.parameter('msc', f'{msc}')
    comp.mesh('mesh1').run()
    ne = comp.mesh('mesh1').getNumElem()
    print(f'\n--- softening sweep, msc={msc} ({ne} elements) ---')
    std = jm.study().create(f'stdEig_{int(msc*100)}')
    std.create('stat', 'Stationary')
    eig = std.create('eig', 'Eigenfrequency')
    eig.set('neigs', JInt(4))
    eig.set('shift', '1e4')
    rows, gev, dstag = [], None, None
    Vlist = [0, 50, 100, 150, 180, 200, 210, 215, 218, 220, 221, 222]
    for Vd in Vlist:
        model.parameter('Vd', f'{Vd}[V]')
        try:
            std.run()
        except Exception as e:
            print(f'  Vd={Vd}: solve failed ({str(e)[:60]})'); continue
        if gev is None:
            gev = jm.result().numerical().create(f'gevE_{int(msc*100)}', 'Global')
            gev.set('expr', ['freq'])
            for t in [str(x) for x in jm.result().dataset().tags()]:
                try:
                    gev.set('data', t); np.array(gev.getReal()[0], dtype=float)
                    dstag = t; break
                except Exception:
                    continue
        gev.set('data', dstag); gev.set('expr', ['freq'])
        f = np.array(gev.getReal()[0], dtype=float)
        rows.append([Vd] + list(f[:4]))
        print(f'  Vd={Vd:4d} V: ' + '  '.join(f'{x/1e3:8.1f} kHz' for x in f[:4]))
    arr = np.array(rows)
    np.savetxt(f'{DATA}/{tag_csv}', arr, delimiter=',',
               header='Vd_V,f1_Hz,f2_Hz,f3_Hz,f4_Hz', comments='')
    # f^2 ~ a*(Vpi^2 - V^2) extrapolation for each mode
    for col, name in [(1, 'axial(f1)'), (2, 'rocking(f2)')]:
        V = arr[:, 0]; y = arr[:, col]**2
        m = V >= 150
        sl, it = np.polyfit(V[m]**2, y[m], 1)
        vpi = np.sqrt(-it/sl) if sl < 0 else float('nan')
        print(f'  V_PI {name} = {vpi:.0f} V')
    return arr, dstag


# (A)+(B) finer-mesh softening, plus coarse re-confirm
arr_fine, dstag = run_softening(MSC_FINE, 'q4b_eig_vs_bias_fine.csv')
arr_coarse, _ = run_softening(1.0, 'q4b_eig_vs_bias_coarse.csv')

# critical-mode shape at last converged bias (finer mesh)
try:
    pg = jm.result().create('pgEigF', 'PlotGroup2D')
    pg.set('data', dstag); pg.set('looplevel', ['2'])  # 2nd eig = rocking
    surf = pg.create('s1', 'Surface'); surf.set('expr', 'solid.disp')
    dfm = surf.create('def1', 'Deform'); dfm.set('expr', ['u', 'v'])
    dfm.set('scaleactive', True); dfm.set('scale', '5e6')
    img = jm.result().export().create('imgEigF', 'Image')
    img.set('sourceobject', 'pgEigF')
    img.set('pngfilename', f'{FIGS}/q4b_critical_mode_fine.png')
    img.set('size', 'manualweb'); img.set('unit', 'px')
    img.set('height', '720'); img.set('width', '1280'); img.run()
    print('critical-mode image exported (fine mesh)')
except Exception as e:
    print('mode image failed:', str(e)[:100])

# (C) stationary continuation sweep on the finer mesh (demonstrates mesh-dep.)
model.parameter('msc', f'{MSC_FINE}')
comp.mesh('mesh1').run()
stdS = jm.study().create('stdSweepF')
ss = stdS.create('stat', 'Stationary')
ss.set('useparam', True)
ss.setIndex('pname', 'Vd', 0)
ss.setIndex('plistarr', 'range(0,5,260)', 0)
ss.set('pcontinuationmode', 'no')
rows = []
try:
    stdS.run()
    soltag = None
    for tg in [str(x) for x in jm.sol().tags()]:
        try:
            if str(jm.sol(tg).study()) == 'stdSweepF':
                soltag = tg
        except Exception:
            pass
    ds = jm.result().dataset().create('dsSwF', 'Solution')
    ds.set('solution', soltag)
    gev = jm.result().numerical().create('gevSwF', 'Global')
    gev.set('data', 'dsSwF'); gev.set('expr', ['uy_mass'])
    uy = np.array(gev.getReal()[0], dtype=float)
    gev.set('expr', ['Vd']); V = np.array(gev.getReal()[0], dtype=float)
    rows = np.column_stack([V, uy])
    np.savetxt(f'{DATA}/q4b_vsweep_fine.csv', rows, delimiter=',',
               header='Vd_V,uy_mass_m', comments='')
    print(f'fine-mesh sweep last converged Vd = {V.max():.0f} V '
          f'(coarse-mesh first submission: 222 V)')
except Exception as e:
    print('fine sweep diverged/failed at some V:', str(e)[:120])

print('DONE q4b rerun')
