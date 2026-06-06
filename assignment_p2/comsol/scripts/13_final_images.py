"""Final report images:
1. q4b_critical_mode.png - the ROCKING mode (mode 2) at 220 V bias, the mode
   whose stiffness collapses at pull-in (re-export; first pass showed mode 1).
2. model_mech_mesh.png - mesh view of the mechanics model (first export was
   blank; use a Mesh dataset, works without stored solutions).
"""
import sys
sys.path.insert(0, '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p2/comsol/scripts')
import mph
from jpype import JInt
from model_lib import MODELS, FIGS

client = mph.start(cores=6)

# ---- 1) rocking mode at 220 V ----
model = client.load(f'{MODELS}/phase_shifter_em.mph')
jm = model.java
model.parameter('Vd', '220[V]')
std = jm.study().create('stdEig2')
std.create('stat', 'Stationary')
eig = std.create('eig', 'Eigenfrequency')
eig.set('neigs', JInt(3))
eig.set('shift', '1e4')
std.run()
print('biased eigen solved at 220 V')

# eigen dataset = the one that evaluates freq
import numpy as np
gev = jm.result().numerical().create('gevF', 'Global')
gev.set('expr', ['freq'])
dstag = None
for t in [str(x) for x in jm.result().dataset().tags()]:
    try:
        gev.set('data', t)
        v = np.array(gev.getReal()[0], dtype=float)
        if len(v) >= 2:
            dstag = t
            print('eigen dataset:', t, 'freqs:', np.round(v / 1e3, 1), 'kHz')
    except Exception:
        continue

pg = jm.result().create('pgM2', 'PlotGroup2D')
pg.set('data', dstag)
pg.set('looplevel', ['2'])          # mode 2 = rocking
surf = pg.create('s1', 'Surface')
surf.set('expr', 'solid.disp')
dfm = surf.create('def1', 'Deform')
dfm.set('expr', ['u', 'v'])         # auto scale
img = jm.result().export().create('imgM2', 'Image')
img.set('sourceobject', 'pgM2')
img.set('pngfilename', f'{FIGS}/q4b_critical_mode.png')
img.set('size', 'manualweb')
img.set('unit', 'px')
img.set('height', '720')
img.set('width', '1280')
img.run()
print('rocking mode image exported')
client.remove(model)

# ---- 2) mech mesh view ----
model = client.load(f'{MODELS}/phase_shifter_mech.mph')
jm = model.java
jm.component('comp1').mesh('mesh1').run()
dsm = jm.result().dataset().create('dsMesh', 'Mesh')
dsm.set('mesh', 'mesh1')
pg = jm.result().create('pgMesh', 'PlotGroup2D')
pg.set('data', 'dsMesh')
pg.create('m1', 'Mesh')
img = jm.result().export().create('imgMesh', 'Image')
img.set('sourceobject', 'pgMesh')
img.set('pngfilename', f'{FIGS}/model_mech_mesh.png')
img.set('size', 'manualweb')
img.set('unit', 'px')
img.set('height', '720')
img.set('width', '1280')
img.run()
print('mesh image exported')
print('DONE final images')
