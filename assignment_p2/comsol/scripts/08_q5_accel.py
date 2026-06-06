"""Q5: maximum acceleration without mass-waveguide contact, Vd = 0,1,10,50 V.

Contact when the modulation bar top reaches the waveguide: uy_bar = +g0 (200 nm).
Body load rho*a_acc (+y, toward waveguide) on all suspended solids.
Linear extrapolation from two solves per Vd + verification solve at a_max.
"""
import sys
sys.path.insert(0, '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p2/comsol/scripts')
import mph
import numpy as np
import csv
from model_lib import MODELS, DATA

G0 = 200e-9

client = mph.start(cores=4)
model = client.load(f'{MODELS}/phase_shifter_em.mph')
jm = model.java
comp = jm.component('comp1')
solid = comp.physics('solid')

# inertial body load (frame acceleration -> +y force on suspended parts)
bl = solid.create('blQ5', 'BodyLoad', 2)
bl.selection().named('selSolids')
bl.set('forceReferenceVolume', ['0', 'solid.rho*a_acc', '0'])

std = jm.study().create('stdQ5')
std.create('stat', 'Stationary')

ds = None
gev = None


def solve_uy(Vd, a):
    global ds, gev
    model.parameter('Vd', f'{Vd}[V]')
    model.parameter('a_acc', f'{a}[m/s^2]')
    std.run()
    soltag = None
    for t in [str(x) for x in jm.sol().tags()]:
        try:
            if str(jm.sol(t).study()) == 'stdQ5':
                soltag = t
        except Exception:
            pass
    if ds is None:
        ds = jm.result().dataset().create('dsQ5', 'Solution')
        gev = jm.result().numerical().create('gevQ5', 'Global')
        gev.set('data', 'dsQ5')
    ds.set('solution', soltag)
    gev.set('expr', ['uy_bar'])
    return float(gev.getReal()[0][0])


rows = []
for Vd in [0, 1, 10, 50]:
    uy0 = solve_uy(Vd, 0.0)
    a_test = 1.0e6
    uy1 = solve_uy(Vd, a_test)
    a_max = a_test * (G0 - uy0) / (uy1 - uy0)
    # verification solve at a_max: uy_bar should hit +200 nm
    uyv = solve_uy(Vd, a_max)
    resid = (uyv - G0) / G0
    rows.append(dict(Vd_V=Vd, uy_rest_m=uy0, a_max_m_s2=a_max,
                     a_max_g=a_max / 9.81, verify_uy_m=uyv,
                     verify_resid=resid))
    print(f'Vd={Vd:3d} V: uy_rest={uy0*1e9:8.4f} nm, a_max={a_max:.4e} m/s^2 '
          f'({a_max/9.81:,.0f} g), verify resid={resid*100:.3f}%')

with open(f'{DATA}/q5_max_acceleration.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

# analytical reference with FEM k and figure-true mass
print('\nReference: a_max ~ (k*(g0+|x_dc|)+dF_tip)/m;'
      ' P1-style with k=36.87, m=3.54e-12: '
      f'{36.87*2e-7/3.5417e-12:.3e} m/s^2 at 0 V')

# CRITICAL: reset swept parameters before saving - leftover a_acc poisons
# every later study loaded from this file (body load stays active!)
model.parameter('a_acc', '0[m/s^2]')
model.parameter('Vd', '50[V]')
model.save(f'{MODELS}/phase_shifter_em.mph')
print('DONE Q5')
