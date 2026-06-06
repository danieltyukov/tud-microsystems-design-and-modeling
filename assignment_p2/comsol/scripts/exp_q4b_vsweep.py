"""Q4b (robust path): continuation voltage sweep 0 -> 450 V in one study.
Continuation from the previous parameter value keeps the solution on the
physical stable branch; the sweep fails just past pull-in -> last converged
voltage brackets V_PI. Partial parametric solutions remain accessible.
"""
import sys
sys.path.insert(0, '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p2/comsol/scripts')
import mph
import numpy as np
from model_lib import MODELS, DATA

client = mph.start(cores=6)
model = client.load(f'{MODELS}/phase_shifter_em.mph')
jm = model.java

std = jm.study().create('stdVs')
stat = std.create('stat', 'Stationary')
stat.set('useparam', 'on')
stat.set('pname', ['Vd'])
stat.set('plistarr', ['range(0,5,450)'])
stat.set('punit', ['V'])
stat.set('sweeptype', 'sparse')
stat.set('pcontinuationmode', 'last')
try:
    std.run()
    print('sweep completed full range (no pull-in below 450 V?)')
except Exception as e:
    print('sweep stopped (pull-in divergence expected):', str(e)[:150])

soltag = None
for t in [str(x) for x in jm.sol().tags()]:
    try:
        if str(jm.sol(t).study()) == 'stdVs':
            soltag = t
    except Exception:
        pass
print('solution:', soltag)
ds = jm.result().dataset().create('dsVs', 'Solution')
ds.set('solution', soltag)
gev = jm.result().numerical().create('gevVs', 'Global')
gev.set('data', 'dsVs')


def geval(expr):
    gev.set('expr', [expr])
    return np.array(gev.getReal()[0], dtype=float)


try:
    V = geval('Vd')
    uy = geval('uy_mass')
    np.savetxt(f'{DATA}/q4b_vsweep_curve.csv',
               np.column_stack([V, uy]), delimiter=',',
               header='Vd_V,uy_mass_m', comments='')
    for i in range(0, len(V), max(1, len(V)//30)):
        print(f'V={V[i]:6.1f}: uy={uy[i]*1e9:9.2f} nm')
    print(f'last converged: V={V[-1]:.1f} V, uy={uy[-1]*1e9:.2f} nm')
    print(f'-> V_PI between {V[-1]:.0f} and {V[-1]+5:.0f} V')
except Exception as e:
    print('extraction failed:', str(e)[:200])

# persist the study (with solutions) for GUI screenshots; reset params first
model.parameter('Vd', '50[V]')
model.save(f'{MODELS}/phase_shifter_em.mph')
print('model saved with stdVs')
print('DONE vsweep')
