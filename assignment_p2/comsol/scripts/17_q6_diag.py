"""Q6 1V diagnostic: load once, try several (smooth, rtol) configs for the 1 V
coupled transient, printing the FULL error for each so we know what the solver
actually objects to. Baseline (1 us smooth, rtol 1e-5) replicates the working
first-submission run (09); the others escalate toward the reviewer's rtol 1e-7.
"""
import sys, traceback
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
    dmp = lemm.create('dmpD', 'Damping', 2)
    dmp.set('DampingType', 'RayleighDamping'); dmp.set('InputParameters', 'AlphaBeta')
    dmp.set('alpha_dM', 'alpha_dm'); dmp.set('beta_dK', 'beta_dk')
except Exception as e:
    print('damping:', str(e)[:50])

step = jm.func().create('stepD', 'Step')
step.set('location', '0'); step.set('smooth', '1e-6')
es.feature('pot1').set('V0', 'Vd*stepD(t[1/s])')

std = jm.study().create('stdD')
t = std.create('time', 'Transient')
t.set('tlist', 'range(0, 2e-7, 2e-4)')
t.set('rtolactive', True); t.set('rtol', '1e-5')
std.createAutoSequences('all')
soltag = None
for tg in [str(x) for x in jm.sol().tags()]:
    try:
        if str(jm.sol(tg).study()) == 'stdD':
            soltag = tg
    except Exception:
        pass
t1 = jm.sol(soltag).feature('t1')
t1.set('storeudot', False)

model.parameter('Vd', '1[V]')
ds = jm.result().dataset().create('dsD', 'Solution'); ds.set('solution', soltag)
gev = jm.result().numerical().create('gevD', 'Global'); gev.set('data', 'dsD')

CONFIGS = [('1us smooth, rtol 1e-5', '1e-6', '1e-5'),
           ('1us smooth, rtol 1e-7', '1e-6', '1e-7'),
           ('200ns smooth, rtol 1e-7', '2e-7', '1e-7')]
for label, sm, rt in CONFIGS:
    step.set('smooth', sm)
    t.set('rtol', rt); t1.set('rtol', rt)
    print(f'\n=== {label} ===', flush=True)
    try:
        std.runNoGen()
        gev.set('expr', ['t']); tt = np.array(gev.getReal()[0], dtype=float)
        gev.set('expr', ['uy_mass']); uy = np.array(gev.getReal()[0], dtype=float)
        xf = uy[-1]; out = np.abs(uy - xf) > 0.05 * abs(xf)
        ts = tt[np.max(np.nonzero(out))] if out.any() else 0.0
        over = np.max(np.abs(uy)) / abs(xf) if xf else float('nan')
        print(f'  OK: final={xf*1e12:.3f} pm  overshoot={over:.2f}x  '
              f'settling(5%)={ts*1e6:.2f} us  npts={len(tt)}', flush=True)
        np.savetxt(f'{DATA}/q6_vstep_1V_fine.csv', np.column_stack([tt, uy]),
                   delimiter=',', header='t_s,uy_mass_m', comments='')
        print('  -> saved q6_vstep_1V_fine.csv', flush=True)
        break
    except Exception as e:
        print('  FULL ERROR:', str(e)[:400], flush=True)
        try:
            print('  solver log tail:', str(jm.sol(soltag).getLastComputationInfo())[:300], flush=True)
        except Exception:
            pass

print('DONE diag', flush=True)
