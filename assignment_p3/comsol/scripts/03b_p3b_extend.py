"""Part 3b extension: full thickness sweep including t > 8 um, final t = 10 um.

The 3D FEM showed the naive t >= Wb = 5 um criterion fails: the lowest
out-of-plane mode is the bar/spine antiphase PITCH (see-saw) mode, lower than
the pure flap the beam formula describes. Crossover where the in-plane
actuation mode becomes the fundamental is near t ~ 8.5 um -> final t = 10 um.

Run: ~/workspace/comsol-mcp/.venv/bin/python 03b_p3b_extend.py
"""
import sys
sys.path.insert(0, '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p3/comsol/scripts')
import csv
import numpy as np
import mph
from jpype import JInt
from model_lib_p3 import MODELS, DATA, FIGS

T_SWEEP = [0.2, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 9.0, 10.0, 12.0]
T_FINAL = 10.0
E, RHO = 170e9, 2330.0
LB, WB = 75e-6, 5e-6

client = mph.start(cores=6)
model = client.load(f'{MODELS}/phase_shifter_p3b_3d.mph')
jm = model.java
print('3D model loaded')

# locate existing eigen dataset / numerical feature
gev = jm.result().numerical('gevF')
rows = []
for t_um in T_SWEEP:
    model.parameter('t_si', f'{t_um}[um]')
    jm.study('stdEig').run()
    gev.set('expr', ['freq'])
    freqs = np.array(gev.getReal()[0], dtype=float)
    gev.set('expr', ['oopfrac'])
    oop = np.array(gev.getReal()[0], dtype=float)
    rows.append([t_um] + [x for pair in zip(freqs[:6], oop[:6]) for x in pair])
    lab = '  '.join(f'{f/1e3:7.1f}k({o:.2f})' for f, o in zip(freqs[:6], oop[:6]))
    print(f't={t_um:4.1f} um: {lab}')

hdr = 't_um,' + ','.join(f'f{i}_Hz,oop{i}' for i in range(1, 7))
np.savetxt(f'{DATA}/p3b_modes_vs_t.csv', np.array(rows), delimiter=',',
           header=hdr, comments='')

# ---- final t: mode maps ----
model.parameter('t_si', f'{T_FINAL}[um]')
jm.study('stdEig').run()
gev.set('expr', ['freq'])
freqs = np.array(gev.getReal()[0], dtype=float)
gev.set('expr', ['oopfrac'])
oop = np.array(gev.getReal()[0], dtype=float)
print('final t modes:', ', '.join(f'{f/1e3:.1f} kHz (oop {o:.2f})'
                                  for f, o in zip(freqs[:5], oop[:5])))

pg = jm.result('pgM')
img = jm.result().export('imgM')
import glob, os
for f in glob.glob(f'{FIGS}/p3b_mode*.png'):
    os.remove(f)
for i in range(5):
    pg.set('looplevel', [str(i + 1)])
    pg.run()
    fn = f'{FIGS}/p3b_mode{i+1}_{freqs[i]/1e3:.0f}kHz.png'
    img.set('pngfilename', fn)
    img.run()
    print('exported', fn)

# ---- static spring constants at final t ----
solid = jm.component('comp1').physics('solid')
bl = solid.feature('blK')
std2 = jm.study('stdK')
ds2 = jm.result().dataset('dsK')
gev2 = jm.result().numerical('gevK')


def solve_k(direction):
    fexpr = ['0', '0', '0']
    fexpr[direction] = 'F_test/(100[um]*30[um]*t_si)'
    bl.set('forceReferenceVolume', fexpr)
    std2.run()
    soltag2 = None
    for tg in [str(x) for x in jm.sol().tags()]:
        try:
            if str(jm.sol(tg).study()) == 'stdK':
                soltag2 = tg
        except Exception:
            pass
    ds2.set('solution', soltag2)
    gev2.set('expr', ['uy_plate' if direction == 1 else 'uz_plate'])
    return 1e-6 / abs(float(gev2.getReal()[0][0]))


k_ip = solve_k(1)
k_oop = solve_k(2)
t_m = T_FINAL * 1e-6
k_ip_ana = 4 * E * t_m * WB**3 / LB**3
k_oop_ana = 4 * E * t_m**3 * WB / LB**3
print(f'k in-plane : FEM {k_ip:8.2f} N/m, analytic {k_ip_ana:8.2f} N/m')
print(f'k out-plane: FEM {k_oop:8.2f} N/m, analytic {k_oop_ana:8.2f} N/m')
with open(f'{DATA}/p3b_spring_constants.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['t_um', 'k_ip_FEM', 'k_ip_analytic', 'k_oop_FEM',
                'k_oop_analytic'])
    w.writerow([T_FINAL, k_ip, k_ip_ana, k_oop, k_oop_ana])

for t in [str(x) for x in jm.sol().tags()]:
    try:
        jm.sol(t).clearSolutionData()
    except Exception:
        pass
model.parameter('t_si', f'{T_FINAL}[um]')
model.save(f'{MODELS}/phase_shifter_p3b_3d.mph')
print('DONE Part 3b extension')
