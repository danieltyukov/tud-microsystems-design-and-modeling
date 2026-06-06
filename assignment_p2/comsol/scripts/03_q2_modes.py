"""Q2: eigenfrequency analysis - first 5 modes, 2D modal displacement maps.

Solves at msc=1 and msc=0.5 (convergence check on f1), exports mode-shape PNGs.
"""
import sys
sys.path.insert(0, '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p2/comsol/scripts')
import mph
import csv
from jpype import JInt
from model_lib import MODELS, DATA, FIGS

client = mph.start(cores=4)
model = client.load(f'{MODELS}/phase_shifter_mech.mph')
jm = model.java
comp = jm.component('comp1')

# eigenfrequency study
std = jm.study().create('stdQ2')
eig = std.create('eig', 'Eigenfrequency')
eig.set('neigsactive', True)
eig.set('neigs', JInt(6))
eig.set('shift', '1[kHz]')

# run once to discover the solution tag, set up dataset + global eval
model.parameter('msc', '1')
comp.mesh('mesh1').run()
jm.study('stdQ2').run()

soltag = None
for t in [str(x) for x in jm.sol().tags()]:
    try:
        if str(jm.sol(t).study()) == 'stdQ2':
            soltag = t
    except Exception:
        pass
print('solution tag for stdQ2:', soltag)

ds = jm.result().dataset().create('dsQ2', 'Solution')
ds.set('solution', soltag)
gev = jm.result().numerical().create('gevF', 'Global')
gev.set('data', 'dsQ2')
gev.set('expr', ['freq'])


def get_freqs():
    vals = gev.getReal()
    return [float(vals[0][j]) for j in range(len(vals[0]))]


results = {}
results[1.0] = get_freqs()
print('msc=1.0: f =', ', '.join(f'{f/1e3:.2f} kHz' for f in results[1.0]))

model.parameter('msc', '0.5')
comp.mesh('mesh1').run()
jm.study('stdQ2').run()
results[0.5] = get_freqs()
print('msc=0.5: f =', ', '.join(f'{f/1e3:.2f} kHz' for f in results[0.5]))

with open(f'{DATA}/q2_eigenfrequencies.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['mode', 'f_msc1_Hz', 'f_msc05_Hz', 'rel_change'])
    for i, (f1v, f2v) in enumerate(zip(results[1.0], results[0.5]), start=1):
        w.writerow([i, f1v, f2v, abs(f2v - f1v) / f2v])

# --- export mode shape images (finest mesh solution is current) ---
pg = jm.result().create('pgMode', 'PlotGroup2D')
pg.set('data', 'dsQ2')
surf = pg.create('surf1', 'Surface')
surf.set('expr', 'solid.disp')
surf.set('unit', 'nm')
defo = surf.create('def1', 'Deform')
defo.set('expr', ['u', 'v'])
defo.set('scaleactive', False)
exp = jm.result().export().create('imgMode', 'pgMode', 'Image')
exp.set('size', 'manualweb')
exp.set('width', '1400')
exp.set('height', '800')

flist = results[0.5]
for i in range(min(5, len(flist))):
    pg.set('solnum', JInt(i + 1))
    pg.run()
    fname = f'{FIGS}/q2_mode{i+1}_{flist[i]/1e3:.0f}kHz.png'
    exp.set('pngfilename', fname)
    exp.run()
    print('exported', fname)

model.save(f'{MODELS}/phase_shifter_mech.mph')
print('DONE Q2')
