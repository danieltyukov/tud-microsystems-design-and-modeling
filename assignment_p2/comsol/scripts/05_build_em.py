"""Build the electromechanical model (Q4-Q6) and save it.
Quick verification solve at Vd=50V: compare displacement against analytical
x = 0.5*V^2*(N*eps0*t/h)/k_FEM and capacitance readout vs N*eps0*t*d0/h.
"""
import sys
sys.path.insert(0, '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p2/comsol/scripts')
import mph
from model_lib import build_em_model, MODELS, FIGS

client = mph.start(cores=4)
model = build_em_model(client)
jm = model.java
comp = jm.component('comp1')

# sanity: list selections sizes
for sel in ['selSolids', 'selTerminal', 'selAirSidesU', 'geom1_cselAir_dom']:
    try:
        ents = comp.selection(sel).entities()
        print(f'{sel}: {len(list(ents))} entities')
    except Exception as e:
        print(f'{sel}: ERROR {e}')

nelem = comp.mesh('mesh1').stat().getNumElem()
print('mesh elements:', nelem)

model.save(f'{MODELS}/phase_shifter_em.mph')
print('SAVED phase_shifter_em.mph')

# ---- verification stationary solve at Vd=50 V ----
model.parameter('Vd', '50[V]')
std = jm.study().create('stdV')
stat = std.create('stat', 'Stationary')
std.run()
print('verification solve OK')

soltag = [str(t) for t in jm.sol().tags()][-1]
ds = jm.result().dataset().create('dsV', 'Solution')
ds.set('solution', soltag)
gev = jm.result().numerical().create('gevV', 'Global')
gev.set('data', 'dsV')
gev.set('expr', ['uy_mass', 'es.Q0_1'])
vals = gev.getReal()
uy = float(vals[0][0])
Q_term = float(vals[1][0])
print(f'uy_mass(50V) = {uy*1e9:.4f} nm')
print(f'analytical x(50V) with k_FEM=36.87: '
      f'{0.5*50**2*7.083e-11/36.87*1e9:.4f} nm')
C = abs(Q_term) / 50.0
C_ana = 40 * 8.854e-12 * 200e-9 * 19e-6 / 1e-6
print(f'C from terminal charge = {C*1e15:.3f} fF (parallel-plate analytical {C_ana*1e15:.3f} fF)')

model.save(f'{MODELS}/phase_shifter_em.mph')
print('DONE EM build')
