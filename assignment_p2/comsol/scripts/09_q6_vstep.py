"""Q6: time-domain step response to a voltage step (1 V and 50 V).

Smoothed step (1 us rise << 2 us oscillation period; sharp ns-steps stall the
coupled ES+moving-mesh DAE); Rayleigh damping as in Q3. Storage-lean: one
voltage per run (sequential), velocities not stored, solution cleared after
extraction - the coupled transient otherwise fills the disk (~4 GB/run).
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

# Rayleigh damping (same recipe as Q3)
lemm = solid.feature('lemm1')
try:
    dmp = lemm.create('dmpQ6', 'Damping', 2)
    dmp.set('DampingType', 'RayleighDamping')
    dmp.set('InputParameters', 'AlphaBeta')
    dmp.set('alpha_dM', 'alpha_dm')
    dmp.set('beta_dK', 'beta_dk')
except Exception as e:
    print('damping exists?', str(e)[:60])

# smoothed voltage step
step = jm.func().create('stepQ6', 'Step')
step.set('location', '1e-6')
step.set('smooth', '2e-6')
es.feature('pot1').set('V0', 'Vd*stepQ6(t[1/s])')

std = jm.study().create('stdQ6')
t = std.create('time', 'Transient')
t.set('tlist', 'range(0, 2e-7, 2e-4)')
t.set('rtolactive', True)
# rtol 1e-4 gave ~1.7x too-fast envelope decay (numerical dissipation);
# 1e-5 is viable with the smooth ramp + tiny initial step
t.set('rtol', '1e-5')

std.createAutoSequences('all')
soltag = None
for tg in [str(x) for x in jm.sol().tags()]:
    try:
        if str(jm.sol(tg).study()) == 'stdQ6':
            soltag = tg
    except Exception:
        pass
sol = jm.sol(soltag)
t1 = sol.feature('t1')
for prop, val in [('initialstepgenalphaactive', 'on'),
                  ('initialstepgenalpha', '1e-10'),
                  ('storeudot', False),
                  ('rhoinf', '0.95')]:   # less generalized-alpha dissipation
    try:
        t1.set(prop, val)
    except Exception as e:
        print(f'{prop} skipped:', str(e)[:60])

ds = jm.result().dataset().create('dsQ6', 'Solution')
ds.set('solution', soltag)
gev = jm.result().numerical().create('gevQ6', 'Global')
gev.set('data', 'dsQ6')


def geval(expr):
    gev.set('expr', [expr])
    return np.array(gev.getReal()[0], dtype=float)


for Vlabel, Vval in [('1V', 1), ('50V', 50)]:
    model.parameter('Vd', f'{Vval}[V]')
    std.runNoGen()
    tt = geval('t')
    uy = geval('uy_mass')
    np.savetxt(f'{DATA}/q6_vstep_{Vlabel}.csv',
               np.column_stack([tt, uy]), delimiter=',',
               header='t_s,uy_mass_m', comments='')
    xf = uy[-1]
    out = np.abs(uy - xf) > 0.05 * abs(xf)
    ts = tt[np.max(np.nonzero(out))] if out.any() else 0.0
    print(f'V={Vlabel}: final={xf*1e9:.4f} nm, settling(5%)={ts*1e6:.2f} us')

# keep the setup, drop the bulky solution, restore stationary-safe V0
sol.clearSolutionData()
model.parameter('Vd', '50[V]')
model.save(f'{MODELS}/phase_shifter_em_step.mph')
print('DONE Q6')
