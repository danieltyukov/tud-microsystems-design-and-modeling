"""Q3: Rayleigh damping for Q=83 + time-domain step response (force input).

beta_dK = 1/(2*pi*f1*Q) gives zeta = 1/(2Q) at the fundamental mode
(stiffness-proportional damping; higher modes are damped more -> clean response).
Step input = constant body load F_test starting from rest.
"""
import sys
sys.path.insert(0, '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p2/comsol/scripts')
import mph
import numpy as np
from jpype import JInt
from model_lib import MODELS, DATA

F1_FEM = 496.26e3   # from Q2 (msc=0.5)
Q = 83.0

client = mph.start(cores=4)
model = client.load(f'{MODELS}/phase_shifter_mech.mph')
jm = model.java
comp = jm.component('comp1')
solid = comp.physics('solid')

model.parameter('f1_fem', f'{F1_FEM}[Hz]')
model.parameter('msc', '1')      # normal mesh for transient
comp.mesh('mesh1').run()

# ---- Rayleigh damping subnode on the linear elastic material ----
lemm = solid.feature('lemm1')
dmp = lemm.create('dmpQ3', 'Damping', 2)
node = None
for c in (model/'physics'/'Solid Mechanics').children():
    if c.type() == 'LinearElasticModel':
        for cc in c.children():
            if 'amping' in cc.name():
                node = cc
print('damping node:', node)
print({k: str(v)[:60] for k, v in node.properties().items()})
dmp.set('DampingType', 'RayleighDamping')
dmp.set('InputParameters', 'AlphaBeta')
dmp.set('alpha_dM', 'alpha_dm')
dmp.set('beta_dK', 'beta_dk')
print('Rayleigh set: alpha=alpha_dm(0), beta=beta_dk=1/(2*pi*f1_fem*Q_target)')

# ---- transient study ----
std = jm.study().create('stdQ3')
t = std.create('time', 'Transient')
t.set('tlist', 'range(0, 1e-7, 4e-4)')
t.set('rtolactive', True)
t.set('rtol', '1e-5')
std.run()
print('transient solved')

# find solution tag, evaluate uy_mass(t) via efficient Global eval
soltag = None
for tg in [str(x) for x in jm.sol().tags()]:
    try:
        if str(jm.sol(tg).study()) == 'stdQ3':
            soltag = tg
    except Exception:
        pass
print('solution tag:', soltag)
ds = jm.result().dataset().create('dsQ3', 'Solution')
ds.set('solution', soltag)
gev = jm.result().numerical().create('gevT', 'Global')
gev.set('data', 'dsQ3')
gev.set('expr', ['uy_mass'])
uy = np.array(gev.getReal(), dtype=float)[0]
gev.set('expr', ['t'])
tt = np.array(gev.getReal(), dtype=float)[0]
print('points:', tt.shape, uy.shape, 'final uy:', uy[-1])

np.savetxt(f'{DATA}/q3_step_response.csv',
           np.column_stack([tt, uy]), delimiter=',',
           header='t_s,uy_mass_m', comments='')

# settling time (5% band around final value)
xf = uy[-1]
out = np.abs(uy - xf) > 0.05 * abs(xf)
ts = tt[np.max(np.nonzero(out))] if out.any() else 0.0
zeta = 1 / (2 * Q)
w1 = 2 * np.pi * F1_FEM
print(f'FEM settling time (5%): {ts*1e6:.2f} us')
print(f'analytic ln(20)/(zeta*w1) with FEM f1: {np.log(20)/(zeta*w1)*1e6:.2f} us')
print(f'P1 prediction (their Q=96.8, w0=4.16e6): 139.4 us')
print(f'final displacement: {xf*1e9:.3f} nm  (static F/k = {1e-6/36.87*1e9:.3f} nm)')

model.save(f'{MODELS}/phase_shifter_mech.mph')
print('DONE Q3')
