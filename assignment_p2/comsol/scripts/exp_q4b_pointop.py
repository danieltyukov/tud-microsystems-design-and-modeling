"""Experiment: does a GE referencing a POINT operator (instead of aveMass)
assemble cleanly? If yes -> use it for the pull-in continuation."""
import mph
from jpype import JInt

client = mph.start(cores=6)
model = client.load('/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p2/comsol/models/phase_shifter_em.mph')
jm = model.java
comp = jm.component('comp1')
es = comp.physics('es')
model.parameter('vset', '-0.025[um]')

# point selection: a corner vertex of the mass plate
s = comp.selection().create('selPt', 'Box')
s.set('entitydim', JInt(0))
s.set('xmin', -50.01e-6); s.set('xmax', -49.99e-6)
s.set('ymin', -100.01e-6); s.set('ymax', -99.99e-6)
s.set('condition', 'inside')

pt = comp.cpl().create('ptMass', 'Integration')
pt.selection().geom('geom1', 0)
pt.selection().named('selPt')
pt.set('opname', 'ptMass')

# verify the selection grabbed exactly one vertex
n = len(comp.selection('selPt').entities())
print('selPt vertices:', n)

ge = es.create('ge1', 'GlobalEquations', -1)
ge.set('DependentVariableQuantity', 'electricpotential')
ge.set('SourceTermQuantity', 'dimensionless')
ge.set('name', ['Vdc'])
ge.set('equation', ['(Vdc-120[V])/1[V]'])
ge.set('initialValueU', ['120[V]'])
es.feature('pot1').set('V0', 'Vdc')

std = jm.study().create('stdT')
std.create('stat', 'Stationary')
std.run()
print('baseline OK')
soltag = [str(t) for t in jm.sol().tags()][-1]
sol = jm.sol(soltag)
sol.feature('v1').set('initmethod', 'sol')
sol.feature('v1').set('initsol', soltag)

ge.set('equation', ['(ptMass(v)-vset)/1[um]'])
try:
    sol.runAll()
    ds = jm.result().dataset().create('dsT', 'Solution')
    ds.set('solution', soltag)
    gev = jm.result().numerical().create('gevT', 'Global')
    gev.set('data', 'dsT')
    for e in ['Vdc', 'uy_mass', 'ptMass(v)']:
        gev.set('expr', [e])
        print(e, '=', float(gev.getReal()[0][0]))
    print('POINT-OP GE STEP OK')
except Exception as e:
    print('POINT-OP GE FAILED:', str(e)[:500])
