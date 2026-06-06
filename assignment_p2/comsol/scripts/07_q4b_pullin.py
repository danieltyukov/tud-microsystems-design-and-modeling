"""Q4b: pull-in voltage via displacement-controlled continuation (v4).

Exact replica of the COMSOL biased_resonator_2d_pull_in tutorial structure:
- POINT PROBE (material frame) measuring v/1[um] on a mass vertex -> global
  variable vmid + helper DOF; the Global Equation references the probe
  variable (GE equation: vrel - vmid).
- Baseline study solves with the GE DISABLED (disabledphysics) at parameter
  Vdc = 120 V -> warm nonzero field (at E=0, dF/dV=0 is singular).
- Pull-in study enables the GE, warm-starts from baseline, sweeps vrel.
- Solver: SEGREGATED in 3 groups [V + probe DOF + ODE] / [u] / [spatial] -
  keeps the GE row out of the moving-mesh block (fully-coupled assembly
  of GE x hyperelastic mesh smoothing produces NaN).
"""
import sys
sys.path.insert(0, '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p2/comsol/scripts')
import mph
import numpy as np
from jpype import JInt
from model_lib import MODELS, DATA

client = mph.start(cores=6)
model = client.load(f'{MODELS}/phase_shifter_em.mph')
jm = model.java
comp = jm.component('comp1')
es = comp.physics('es')

model.parameter('Vdc', '120[V]')   # baseline bias (GE shadows it when active)
model.parameter('vrel', '-0.025')  # prescribed displacement in um (dimensionless)

# --- point probe on a mass-plate corner vertex (material frame) ---
s = comp.selection().create('selPt', 'Box')
s.set('entitydim', JInt(0))
s.set('xmin', -50.01e-6); s.set('xmax', -49.99e-6)
s.set('ymin', -100.01e-6); s.set('ymax', -99.99e-6)
s.set('condition', 'inside')
prb = comp.probe().create('point1', 'Point')
prb.selection().named('selPt')
prb.set('expr', 'v/1[um]')
prb.set('probename', 'vmid')
prb.set('frame', 'material')

# --- global equation (tutorial-style) ---
ge = es.create('ge1', 'GlobalEquations', -1)
ge.set('DependentVariableQuantity', 'electricpotential')
ge.set('SourceTermQuantity', 'dimensionless')
ge.set('name', ['Vdc'])
ge.set('equation', ['vrel-vmid'])
ge.set('initialValueU', ['0'])
ge.set('description', ['DC bias voltage'])
es.feature('pot1').set('V0', 'Vdc')

# --- baseline study: GE disabled, fixed 120 V ---
std1 = jm.study().create('stdQ4bInit')
st1 = std1.create('stat', 'Stationary')
st1.set('useadvanceddisable', 'on')   # REQUIRED for disabledphysics to apply
st1.set('disabledphysics', ['es/ge1'])
std1.run()
print('baseline solved (GE disabled, Vdc=120 V)')

# --- pull-in study: GE active, warm start, sweep vrel ---
std2 = jm.study().create('stdQ4b')
st2 = std2.create('stat', 'Stationary')
st2.set('useinitsol', 'on')
st2.set('initmethod', 'sol')
st2.set('initstudy', 'stdQ4bInit')
st2.set('solnum', 'auto')
st2.set('useparam', 'on')
st2.set('pname', ['vrel'])
st2.set('plistarr', ['range(-0.025,-0.025,-0.9)'])
st2.set('punit', [''])
st2.set('sweeptype', 'sparse')
st2.set('pcontinuationmode', 'last')

# generate the default solver sequence WITHOUT running, then restructure
try:
    std2.createAutoSequences('all')
except Exception:
    try:
        jm.sol().createAutoSequence('stdQ4b')
    except Exception as e:
        print('auto-seq generation issue:', str(e)[:120])

soltag2 = None
for t in [str(x) for x in jm.sol().tags()]:
    try:
        if str(jm.sol(t).study()) == 'stdQ4b':
            soltag2 = t
    except Exception:
        pass
print('pull-in solution tag:', soltag2)
sol2 = jm.sol(soltag2)
s1 = sol2.feature('s1')

# list variables to find the probe helper DOF tag
vtags = [str(t) for t in sol2.feature('v1').feature().tags()]
print('solver variables:', vtags)
ode_vars = [v for v in vtags if 'ODE' in v or v.endswith('_D')]
es_group = ['comp1_V'] + ode_vars
solid_group = [v for v in vtags if v in ('comp1_u', 'comp1_solid_wZ')]
mesh_group = [v for v in vtags if 'spatial' in v]
print('groups:', es_group, solid_group, mesh_group)

# restructure the nonlinear solver: drop FullyCoupled, configure the
# Segregated attribute (reuse the auto-generated one) with tutorial groups
se = None
for ft in [str(t) for t in s1.feature().tags()]:
    ftype = str(s1.feature(ft).getType())
    if ftype == 'FullyCoupled':
        s1.feature().remove(ft)
    elif ftype == 'Segregated':
        se = s1.feature(ft)
if se is None:
    se = s1.create('seQ4b', 'Segregated')
for ft in [str(t) for t in se.feature().tags()]:
    try:
        se.feature().remove(ft)
    except Exception:
        print('could not remove se step', ft)
for tag, grp in [('ssA', es_group), ('ssB', solid_group), ('ssC', mesh_group)]:
    ss = se.create(tag, 'SegregatedStep')
    ss.set('segvar', grp)
    ss.set('linsolver', 'dDef')

std2.runNoGen()
print('pull-in sweep solved')

# --- extract V(x) ---
dstag = None
for t in [str(x) for x in jm.result().dataset().tags()]:
    try:
        if str(jm.result().dataset(t).getString('solution')) == soltag2:
            dstag = t
    except Exception:
        pass
if dstag is None:
    ds = jm.result().dataset().create('dsQ4b', 'Solution')
    ds.set('solution', soltag2)
    dstag = 'dsQ4b'
gev = jm.result().numerical().create('gevQ4b', 'Global')
gev.set('data', dstag)


def geval(expr):
    gev.set('expr', [expr])
    return np.array(gev.getReal()[0], dtype=float)


Vdc = geval('Vdc')
uy = geval('uy_mass')
np.savetxt(f'{DATA}/q4b_pullin_curve.csv',
           np.column_stack([uy, Vdc]), delimiter=',',
           header='uy_mass_m,Vdc_V', comments='')
ipk = int(np.argmax(np.abs(Vdc)))
print('\n  x (nm)    V (V)')
for i in range(len(uy)):
    mark = '  <-- PULL-IN' if i == ipk else ''
    print(f'{uy[i]*1e9:9.2f}  {abs(Vdc[i]):8.2f}{mark}')
print(f'\nV_PI(FEM) = {abs(Vdc[ipk]):.2f} V at x = {uy[ipk]*1e9:.1f} nm')
print('P1 analytical: V_PI = 347.1 V at x = -404 nm (tip-gap model)')

model.save(f'{MODELS}/phase_shifter_em_pullin.mph')
print('DONE Q4b')
