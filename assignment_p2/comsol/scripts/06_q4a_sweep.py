"""Q4a: stationary displacement vs comb-drive voltage, Vd = 0..150 V.

Auxiliary sweep with continuation; records uy_mass, uy_bar and terminal
charge (capacitance). Output: data/q4a_displacement_vs_voltage.csv
"""
import sys
sys.path.insert(0, '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p2/comsol/scripts')
import mph
import numpy as np
from model_lib import MODELS, DATA

client = mph.start(cores=4)
model = client.load(f'{MODELS}/phase_shifter_em.mph')
jm = model.java

# clean any previous instance (script re-runnable)
for tag in ['stdQ4a']:
    try:
        jm.study().remove(tag)
    except Exception:
        pass

std = jm.study().create('stdQ4a')
stat = std.create('stat', 'Stationary')
stat.set('useparam', 'on')
stat.set('pname', ['Vd'])
stat.set('plistarr', ['range(0,10,150)'])
stat.set('punit', ['V'])
stat.set('sweeptype', 'sparse')
stat.set('pcontinuationmode', 'last')
std.run()
print('sweep solved')

soltag = None
for t in [str(x) for x in jm.sol().tags()]:
    try:
        if str(jm.sol(t).study()) == 'stdQ4a':
            soltag = t
    except Exception:
        pass
print('solution:', soltag)

ds = jm.result().dataset().create('dsQ4a', 'Solution')
ds.set('solution', soltag)
gev = jm.result().numerical().create('gevQ4a', 'Global')
gev.set('data', 'dsQ4a')


def geval(expr):
    gev.set('expr', [expr])
    return np.array(gev.getReal()[0], dtype=float)


Vd = geval('Vd')
uy = geval('uy_mass')
uyb = geval('uy_bar')
Qt = geval('es.Q0_1')
C = np.where(Vd > 0, np.abs(Qt) / np.maximum(Vd, 1e-12), np.nan)

# analytical P1 prediction: x = 0.5 V^2 (N eps0 t / h) / k_P1
k_P1 = 40.2963
dCdx_P1 = 40 * 8.854187817e-12 * 200e-9 / 1e-6
x_P1 = 0.5 * Vd**2 * dCdx_P1 / k_P1

out = np.column_stack([Vd, uy, uyb, C, -x_P1])
np.savetxt(f'{DATA}/q4a_displacement_vs_voltage.csv', out, delimiter=',',
           header='Vd_V,uy_mass_m,uy_bar_m,C_F,x_P1_analytical_m(downward=neg)',
           comments='')
for i in range(len(Vd)):
    print(f'Vd={Vd[i]:5.0f} V: uy={uy[i]*1e9:8.4f} nm  '
          f'(P1 ana {-x_P1[i]*1e9:8.4f} nm)  C={C[i]*1e15 if Vd[i]>0 else float("nan"):.4f} fF')

model.save(f'{MODELS}/phase_shifter_em.mph')
print('DONE Q4a')
