"""Q6 final re-run (addresses reviewer #6): directly resolve the 1 V voltage-step
settling time. Root cause of the first submission's bad 1 V decay was NOT just
rtol -- it was generalized-alpha numerical dissipation swamping the physical
envelope at picometre amplitude. Fixes:
  - rho_inf = 0.98 (near-zero algorithmic damping; was 0.95/0.97),
  - rtol 1e-6 (tighter than 1e-5; 1e-7 makes the coupled DAE fail to converge),
  - 200 ns smoothed step (a fully sharp step stalls the ES+moving-mesh DAE; 200 ns
    is ~10% of the 2 us period, so the ~2x overshoot of a true step is preserved),
  - es thickness = 1 (teaching-team tip), msc 0.8 mesh.
Runs 1 V then 50 V; saves each immediately. -> q6_vstep_{1V,50V}_fine.csv
"""
import sys
sys.path.insert(0, '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p2/comsol/scripts')
import mph
import numpy as np
from model_lib import MODELS, DATA

F1_FEM = 496.26e3
client = mph.start(cores=6)
model = client.load(f'{MODELS}/phase_shifter_em.mph')
jm = model.java
comp = jm.component('comp1')
solid = comp.physics('solid')
es = comp.physics('es')
model.parameter('f1_fem', f'{F1_FEM}[Hz]')
model.parameter('msc', '0.8')
es.prop('d').set('d', '1')
comp.mesh('mesh1').run()
print('mesh:', comp.mesh('mesh1').getNumElem(), flush=True)

lemm = solid.feature('lemm1')
try:
    dmp = lemm.create('dmpF', 'Damping', 2)
    dmp.set('DampingType', 'RayleighDamping'); dmp.set('InputParameters', 'AlphaBeta')
    dmp.set('alpha_dM', 'alpha_dm'); dmp.set('beta_dK', 'beta_dk')
except Exception as e:
    print('damping:', str(e)[:50], flush=True)

step = jm.func().create('stepF', 'Step')
step.set('location', '0'); step.set('smooth', '2e-7')   # 200 ns
es.feature('pot1').set('V0', 'Vd*stepF(t[1/s])')

std = jm.study().create('stdF')
t = std.create('time', 'Transient')
t.set('tlist', 'range(0, 1e-7, 2e-4)')
t.set('rtolactive', True); t.set('rtol', '1e-6')
std.createAutoSequences('all')
soltag = None
for tg in [str(x) for x in jm.sol().tags()]:
    try:
        if str(jm.sol(tg).study()) == 'stdF':
            soltag = tg
    except Exception:
        pass
t1 = jm.sol(soltag).feature('t1')
for prop, val in [('storeudot', False), ('rhoinf', '0.98'), ('rtol', '1e-6'),
                  ('initialstepgenalphaactive', 'on'),
                  ('initialstepgenalpha', '1e-9')]:
    try:
        t1.set(prop, val)
    except Exception as e:
        print(f'{prop} skipped:', str(e)[:50], flush=True)

ds = jm.result().dataset().create('dsF', 'Solution'); ds.set('solution', soltag)
gev = jm.result().numerical().create('gevF', 'Global'); gev.set('data', 'dsF')


def geval(expr):
    gev.set('expr', [expr]); return np.array(gev.getReal()[0], dtype=float)


for Vlabel, Vval in [('1V', 1), ('50V', 50)]:
    model.parameter('Vd', f'{Vval}[V]')
    print(f'\n=== {Vlabel} transient ===', flush=True)
    try:
        std.runNoGen()
    except Exception as e:
        print(f'{Vlabel} FAILED:', str(e)[:300], flush=True); continue
    tt = geval('t'); uy = geval('uy_mass')
    np.savetxt(f'{DATA}/q6_vstep_{Vlabel}_fine.csv', np.column_stack([tt, uy]),
               delimiter=',', header='t_s,uy_mass_m', comments='')
    xf = uy[-1]; out = np.abs(uy - xf) > 0.05 * abs(xf)
    ts = tt[np.max(np.nonzero(out))] if out.any() else 0.0
    over = np.max(np.abs(uy)) / abs(xf) if xf else float('nan')
    unit = 1e12 if abs(xf) < 1e-10 else 1e9
    u = 'pm' if unit == 1e12 else 'nm'
    print(f'{Vlabel}: final={xf*unit:.3f} {u}  overshoot={over:.2f}x  '
          f'settling(5%)={ts*1e6:.2f} us  npts={len(tt)}', flush=True)
    jm.sol(soltag).clearSolutionData()

model.parameter('Vd', '0[V]')
print('DONE q6 final', flush=True)
