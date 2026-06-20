"""Q6 1V direct re-run -- LIGHTWEIGHT to avoid overheating the laptop.
Only 2 cores (low thermal load), default mesh, 1 V only (50 V already on file).
Config that fixes the picometre decay artefact: rho_inf=0.98 (minimal
generalized-alpha dissipation) + rtol 1e-6 + 200 ns smoothed step + es
thickness=1. -> q6_vstep_1V_fine.csv
"""
import sys
sys.path.insert(0, '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p2/comsol/scripts')
import mph
import numpy as np
from model_lib import MODELS, DATA

F1_FEM = 496.26e3
client = mph.start(cores=2)            # 2 cores only -- gentle on the laptop
print('COMSOL up (2 cores):', client.version, flush=True)
model = client.load(f'{MODELS}/phase_shifter_em.mph')
jm = model.java
comp = jm.component('comp1')
solid = comp.physics('solid')
es = comp.physics('es')
model.parameter('f1_fem', f'{F1_FEM}[Hz]')
model.parameter('msc', '1.0')          # default (coarser) mesh -- lighter
es.prop('d').set('d', '1')
comp.mesh('mesh1').run()
print('mesh:', comp.mesh('mesh1').getNumElem(), flush=True)

lemm = solid.feature('lemm1')
try:
    dmp = lemm.create('dmpL', 'Damping', 2)
    dmp.set('DampingType', 'RayleighDamping'); dmp.set('InputParameters', 'AlphaBeta')
    dmp.set('alpha_dM', 'alpha_dm'); dmp.set('beta_dK', 'beta_dk')
except Exception as e:
    print('damping:', str(e)[:50], flush=True)

step = jm.func().create('stepL', 'Step')
step.set('location', '5e-7'); step.set('smooth', '1e-6')   # stable 1 us rise
es.feature('pot1').set('V0', 'Vd*stepL(t[1/s])')

std = jm.study().create('stdL')
t = std.create('time', 'Transient')
t.set('tlist', 'range(0, 2e-7, 2e-4)')   # 1000 pts -> lean storage
t.set('rtolactive', True); t.set('rtol', '1e-5')
std.createAutoSequences('all')
soltag = None
for tg in [str(x) for x in jm.sol().tags()]:
    try:
        if str(jm.sol(tg).study()) == 'stdL':
            soltag = tg
    except Exception:
        pass
t1 = jm.sol(soltag).feature('t1')
# rho_inf -> 1 removes the generalized-alpha numerical dissipation that made the
# pm-scale 1V decay too fast; 1 us smooth + rtol 1e-5 keeps the solve stable/light
for prop, val in [('storeudot', False), ('rhoinf', '0.995'), ('rtol', '1e-5')]:
    try:
        t1.set(prop, val)
    except Exception as e:
        print(f'{prop} skip:', str(e)[:40], flush=True)

model.parameter('Vd', '1[V]')
print('=== 1V transient (2 cores, rtol 1e-6) ===', flush=True)
std.runNoGen()
ds = jm.result().dataset().create('dsL', 'Solution'); ds.set('solution', soltag)
gev = jm.result().numerical().create('gevL', 'Global'); gev.set('data', 'dsL')
gev.set('expr', ['t']); tt = np.array(gev.getReal()[0], dtype=float)
gev.set('expr', ['uy_mass']); uy = np.array(gev.getReal()[0], dtype=float)
np.savetxt(f'{DATA}/q6_vstep_1V_fine.csv', np.column_stack([tt, uy]),
           delimiter=',', header='t_s,uy_mass_m', comments='')
xf = uy[-1]; out = np.abs(uy - xf) > 0.05 * abs(xf)
ts = tt[np.max(np.nonzero(out))] if out.any() else 0.0
over = np.max(np.abs(uy)) / abs(xf) if xf else float('nan')
print(f'1V: final={xf*1e12:.3f} pm  overshoot={over:.2f}x  settling(5%)={ts*1e6:.2f} us', flush=True)
print('DONE q6 light', flush=True)
