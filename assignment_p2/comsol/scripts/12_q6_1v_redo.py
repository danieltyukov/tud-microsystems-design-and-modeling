"""Q6 redo, 1 V only: at pm amplitudes the rtol=1e-4 adaptive stepper takes
huge steps and generalized-alpha dissipation wrecks the envelope (apparent
settling 3x too fast - unphysical for a linear system). Tighten rtol + rhoinf.
Overwrites data/q6_vstep_1V.csv. 50 V trace from the first pass is kept.
"""
import sys
sys.path.insert(0, '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p2/comsol/scripts')
import mph
import numpy as np
from model_lib import MODELS, DATA

client = mph.start(cores=6)
model = client.load(f'{MODELS}/phase_shifter_em_step.mph')
jm = model.java

std = jm.study('stdQ6')
std.feature('time').set('rtol', '1e-5')

soltag = None
for tg in [str(x) for x in jm.sol().tags()]:
    try:
        if str(jm.sol(tg).study()) == 'stdQ6':
            soltag = tg
    except Exception:
        pass
print('sol:', soltag)
t1 = jm.sol(soltag).feature('t1')
for prop, val in [('rtol', '1e-5')]:   # rhoinf override made onset diverge
    try:
        t1.set(prop, val)
        print(f'{prop} set')
    except Exception as e:
        print(f'{prop} skipped:', str(e)[:60])

model.parameter('Vd', '1[V]')
std.runNoGen()
print('1V transient solved')

ds = jm.result().dataset().create('dsR', 'Solution')
ds.set('solution', soltag)
gev = jm.result().numerical().create('gevR', 'Global')
gev.set('data', 'dsR')


def geval(expr):
    gev.set('expr', [expr])
    return np.array(gev.getReal()[0], dtype=float)


tt = geval('t')
uy = geval('uy_mass')
np.savetxt(f'{DATA}/q6_vstep_1V.csv',
           np.column_stack([tt, uy]), delimiter=',',
           header='t_s,uy_mass_m', comments='')
xf = uy[-1]
out = np.abs(uy - xf) > 0.05 * abs(xf)
ts = tt[np.max(np.nonzero(out))] if out.any() else 0.0
print(f'1V redo: final={xf*1e12:.3f} pm, settling(5%)={ts*1e6:.2f} us')

jm.sol(soltag).clearSolutionData()
model.parameter('Vd', '50[V]')
model.save()
print('DONE 1V redo')
