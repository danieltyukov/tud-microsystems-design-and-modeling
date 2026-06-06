"""Part 3a.4 (i),(ii),(v): build mechanics model, spring constant with mesh
convergence, eigenfrequencies + 5 mode maps, Rayleigh-damped force-step
transient (rise time + settling).

Run: ~/workspace/comsol-mcp/.venv/bin/python 01_mech_all.py
"""
import sys
sys.path.insert(0, '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p3/comsol/scripts')
import csv
import numpy as np
import mph
from jpype import JInt
import model_lib_p3 as lib
from model_lib_p3 import MODELS, DATA, FIGS

Q_TARGET = 80.7          # analytic damping model (air + clamping) at f_r

client = mph.start(cores=6)
model = lib.build_mech_model(client)
jm = model.java
comp = jm.component('comp1')
solid = comp.physics('solid')
print('model built')

# ---------------- (i) spring constant -----------------
bl = solid.create('blQk', 'BodyLoad', 2)
bl.selection().named('geom1_cselMass_dom')
bl.set('forceReferenceVolume', ['0', '-F_test/(L*Wm*t_si)', '0'])

std = jm.study().create('stdK')
std.create('stat', 'Stationary')

rows = []
for msc in [2.0, 1.0, 0.5]:
    model.parameter('msc', str(msc))
    comp.mesh('mesh1').run()
    nelem = comp.mesh('mesh1').stat().getNumElem()
    jm.study('stdK').run()
    uy_mass = float(model.evaluate('uy_mass', 'm'))
    uy_bar = float(model.evaluate('uy_bar', 'm'))
    m_sh = float(model.evaluate('m_shuttle', 'kg'))
    k_mass = 1e-6 / abs(uy_mass)
    rows.append(dict(msc=msc, nelem=int(nelem), uy_mass=uy_mass,
                     uy_bar=uy_bar, k_mass=k_mass,
                     k_bar=1e-6 / abs(uy_bar), m_shuttle=m_sh))
    print(f'msc={msc}: nelem={nelem}, uy={uy_mass:.4e} m, k={k_mass:.4f} N/m, '
          f'm_shuttle={m_sh:.4e} kg')

with open(f'{DATA}/p3a_q4i_spring_constant.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
k_ana = 4 * 170e9 * 200e-9 * (1.85e-6)**3 / (40e-6)**3
print(f'analytic k = {k_ana:.3f} N/m; FEM (finest) = {rows[-1]["k_mass"]:.3f} '
      f'({(rows[-1]["k_mass"]/k_ana-1)*100:+.1f} %)')
# disable the test load for subsequent studies
bl.active(False)

# ---------------- (ii) eigenfrequencies + maps -----------------
std2 = jm.study().create('stdEig')
eig = std2.create('eig', 'Eigenfrequency')
eig.set('neigsactive', True)
eig.set('neigs', JInt(8))
eig.set('shift', '1[kHz]')

model.parameter('msc', '1')
comp.mesh('mesh1').run()
jm.study('stdEig').run()

soltag = None
for t in [str(x) for x in jm.sol().tags()]:
    try:
        if str(jm.sol(t).study()) == 'stdEig':
            soltag = t
    except Exception:
        pass
ds = jm.result().dataset().create('dsEig', 'Solution')
ds.set('solution', soltag)
gev = jm.result().numerical().create('gevF', 'Global')
gev.set('data', 'dsEig')
gev.set('expr', ['freq'])


def get_freqs():
    vals = gev.getReal()
    return [float(vals[0][j]) for j in range(len(vals[0]))]


res = {1.0: get_freqs()}
print('msc=1.0: f =', ', '.join(f'{f/1e3:.2f}k' for f in res[1.0]))
model.parameter('msc', '0.5')
comp.mesh('mesh1').run()
jm.study('stdEig').run()
res[0.5] = get_freqs()
print('msc=0.5: f =', ', '.join(f'{f/1e3:.2f}k' for f in res[0.5]))

with open(f'{DATA}/p3a_q4ii_eigenfrequencies.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['mode', 'f_msc1_Hz', 'f_msc05_Hz', 'rel_change'])
    for i, (a, b) in enumerate(zip(res[1.0], res[0.5]), start=1):
        w.writerow([i, a, b, abs(b - a) / b])

pg = jm.result().create('pgMode', 'PlotGroup2D')
pg.set('data', 'dsEig')
surf = pg.create('surf1', 'Surface')
surf.set('expr', 'solid.disp')
surf.set('unit', 'nm')
defo = surf.create('def1', 'Deform')
defo.set('expr', ['u', 'v'])
defo.set('scaleactive', False)
exp = jm.result().export().create('imgMode', 'pgMode', 'Image')
exp.set('size', 'manualweb')
exp.set('width', '1400')
exp.set('height', '900')
for i in range(5):
    pg.set('solnum', JInt(i + 1))
    pg.run()
    fname = f'{FIGS}/p3a_mode{i+1}_{res[0.5][i]/1e3:.0f}kHz.png'
    exp.set('pngfilename', fname)
    exp.run()
    print('exported', fname)

# mesh picture
pgm = jm.result().create('pgMesh', 'PlotGroup2D')
pgm.set('data', 'dsEig')
mp = pgm.create('mesh1', 'Mesh')
pgm.run()
expm = jm.result().export().create('imgMesh', 'pgMesh', 'Image')
expm.set('size', 'manualweb')
expm.set('width', '1400')
expm.set('height', '900')
expm.set('pngfilename', f'{FIGS}/p3a_model_mesh.png')
expm.run()

# ---------------- (v) step response -----------------
f1 = res[0.5][0]
model.parameter('f1_fem', f'{f1}[Hz]')
model.parameter('Q_target', str(Q_TARGET))
model.parameter('msc', '1')
comp.mesh('mesh1').run()

lemm = solid.feature('lemm1')
dmp = lemm.create('dmpStep', 'Damping', 2)
dmp.set('DampingType', 'RayleighDamping')
dmp.set('InputParameters', 'AlphaBeta')
dmp.set('alpha_dM', 'alpha_dm')
dmp.set('beta_dK', 'beta_dk')
bl.active(True)   # re-enable the 1 uN load as the step input

std3 = jm.study().create('stdStep')
tstep = std3.create('time', 'Transient')
tstep.set('tlist', 'range(0, 5e-8, 3e-4)')
tstep.set('rtolactive', True)
tstep.set('rtol', '1e-5')
std3.run()
print('transient solved')

soltag = None
for t in [str(x) for x in jm.sol().tags()]:
    try:
        if str(jm.sol(t).study()) == 'stdStep':
            soltag = t
    except Exception:
        pass
ds3 = jm.result().dataset().create('dsStep', 'Solution')
ds3.set('solution', soltag)
gev3 = jm.result().numerical().create('gevT', 'Global')
gev3.set('data', 'dsStep')
gev3.set('expr', ['uy_mass'])
uy = np.array(gev3.getReal(), dtype=float)[0]
gev3.set('expr', ['t'])
tt = np.array(gev3.getReal(), dtype=float)[0]
np.savetxt(f'{DATA}/p3a_q4v_step_response.csv',
           np.column_stack([tt, uy]), delimiter=',',
           header='t_s,uy_mass_m', comments='')

xf = uy[-1]
out = np.abs(uy - xf) > 0.05 * abs(xf)
ts = tt[np.max(np.nonzero(out))] if out.any() else 0.0
# 10-90% rise time of the first transition
i10 = np.argmax(np.abs(uy) >= 0.1 * abs(xf))
i90 = np.argmax(np.abs(uy) >= 0.9 * abs(xf))
zeta = 1 / (2 * Q_TARGET)
w1 = 2 * np.pi * f1
print(f'final uy = {xf*1e9:.3f} nm (static F/k = {1e-6/rows[-1]["k_mass"]*1e9:.3f} nm)')
print(f'rise time 10-90%   = {(tt[i90]-tt[i10])*1e6:.3f} us')
print(f'settling time (5%) = {ts*1e6:.2f} us '
      f'(analytic ln20/(zeta w1) = {np.log(20)/(zeta*w1)*1e6:.2f} us)')

# reset swept parameters before saving (P2 gotcha: parameter contamination)
model.parameter('msc', '1')
model.save(f'{MODELS}/phase_shifter_p3_mech.mph')
print('DONE mech (i),(ii),(v)')
