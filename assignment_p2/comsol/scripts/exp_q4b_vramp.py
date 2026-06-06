"""Q4b insurance: bracket V_PI by ramping Vd until the stationary solve
diverges. Warm-starts each step from the previous solution (v1 self-init).
Adaptive bisection refines the bracket to ~1 V.
"""
import mph
import numpy as np

DATA = '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p2/data'
client = mph.start(cores=6)
model = client.load('/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p2/comsol/models/phase_shifter_em.mph')
jm = model.java

std = jm.study().create('stdVr')
std.create('stat', 'Stationary')
model.parameter('Vd', '150[V]')
std.run()
soltag = [str(t) for t in jm.sol().tags()][-1]
sol = jm.sol(soltag)
sol.feature('v1').set('initmethod', 'sol')
sol.feature('v1').set('initsol', soltag)

ds = jm.result().dataset().create('dsVr', 'Solution')
ds.set('solution', soltag)
gev = jm.result().numerical().create('gevVr', 'Global')
gev.set('data', 'dsVr')


def geval(expr):
    gev.set('expr', [expr])
    return float(gev.getReal()[0][0])


rows = [(150.0, geval('uy_mass'))]
print(f'V=150: uy={rows[0][1]*1e9:.2f} nm')

# coarse upward ramp, then bisection between last-good and first-fail
V, dV = 150.0, 25.0
V_good, V_fail = 150.0, None
while True:
    Vtry = V + dV
    model.parameter('Vd', f'{Vtry}[V]')
    try:
        sol.runAll()
        uy = geval('uy_mass')
        rows.append((Vtry, uy))
        print(f'V={Vtry:.2f}: uy={uy*1e9:.2f} nm')
        V_good, V = Vtry, Vtry
        if V_fail is not None:
            dV = (V_fail - V_good) / 2
    except Exception:
        print(f'V={Vtry:.2f}: DIVERGED')
        V_fail = Vtry
        dV = (V_fail - V_good) / 2
        # re-solve at last good V to restore a converged warm-start state
        model.parameter('Vd', f'{V_good}[V]')
        try:
            sol.runAll()
        except Exception:
            print('warm-state restore failed'); break
    if V_fail is not None and (V_fail - V_good) <= 2.0:
        break
    if Vtry > 800:
        break

arr = np.array(rows)
np.savetxt(f'{DATA}/q4b_vramp_bracket.csv', arr, delimiter=',',
           header='Vd_V,uy_mass_m', comments='')
print(f'\nBRACKET: V_PI between {V_good:.1f} V (converged) and '
      f'{V_fail:.1f} V (diverged)' if V_fail else 'no divergence below 800 V')
print('DONE vramp')
