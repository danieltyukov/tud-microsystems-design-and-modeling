"""Q4b verification: prestressed eigenfrequencies vs bias voltage.

For each Vd: stationary solve (warm continuation), then eigenfrequency
linearized about that state (includes electrostatic softening via the EM
coupling). The voltage where f1^2 -> 0 is the true instability; the mode
shape identifies the mechanism (axial shuttle pull-in vs comb-finger side
pull-in). Exports mode image at the highest converged bias.
"""
import sys
sys.path.insert(0, '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p2/comsol/scripts')
import mph
import numpy as np
from jpype import JInt
from model_lib import MODELS, DATA, FIGS

client = mph.start(cores=6)
model = client.load(f'{MODELS}/phase_shifter_em.mph')
jm = model.java

std = jm.study().create('stdEig')
stat = std.create('stat', 'Stationary')
eig = std.create('eig', 'Eigenfrequency')
eig.set('neigs', JInt(4))
eig.set('shift', '1e4')

rows = []
gev = None
dstag_ok = None
for Vd in [0, 50, 100, 150, 180, 200, 210, 215, 220]:
    model.parameter('Vd', f'{Vd}[V]')
    try:
        std.run()
    except Exception as e:
        print(f'Vd={Vd}: solve failed ({str(e)[:80]}) - skipping')
        continue
    # locate the auto-created dataset that carries the eigen solution:
    # with a stat+eig study only that dataset evaluates 'freq'
    if gev is None:
        gev = jm.result().numerical().create('gevEig', 'Global')
        gev.set('expr', ['freq'])
        for t in [str(x) for x in jm.result().dataset().tags()]:
            try:
                gev.set('data', t)
                np.array(gev.getReal()[0], dtype=float)
                dstag_ok = t
                break
            except Exception:
                continue
        print('eigen dataset:', dstag_ok)
    gev.set('data', dstag_ok)
    gev.set('expr', ['freq'])
    freqs = np.array(gev.getReal()[0], dtype=float)
    rows.append([Vd] + list(freqs[:4]))
    print(f'Vd={Vd:4d} V: f = ' + '  '.join(f'{f/1e3:9.2f} kHz' for f in freqs[:4]))

arr = np.array(rows)
np.savetxt(f'{DATA}/q4b_eig_vs_bias.csv', arr, delimiter=',',
           header='Vd_V,f1_Hz,f2_Hz,f3_Hz,f4_Hz', comments='')

# mode-1 shape image at the highest converged bias (mechanism identification)
try:
    pg = jm.result().create('pgEig', 'PlotGroup2D')
    pg.set('data', dstag_ok)
    pg.set('looplevel', ['1'])
    surf = pg.create('s1', 'Surface')
    surf.set('expr', 'solid.disp')
    dfm = surf.create('def1', 'Deform')
    dfm.set('expr', ['u', 'v'])
    dfm.set('scaleactive', True)
    dfm.set('scale', '5e6')
    img = jm.result().export().create('imgEig', 'Image')
    img.set('sourceobject', 'pgEig')
    img.set('pngfilename', f'{FIGS}/q4b_critical_mode.png')
    img.set('size', 'manualweb')
    img.set('unit', 'px')
    img.set('height', '720')
    img.set('width', '1280')
    img.run()
    print('mode image exported')
except Exception as e:
    print('image export failed:', str(e)[:100])

# quadratic extrapolation of f1^2 -> 0 for V_PI estimate
if len(arr) >= 3:
    V = arr[:, 0]
    f1sq = arr[:, 1]**2
    c = np.polyfit(V[-4:]**2, f1sq[-4:], 1)   # f1^2 ~ a + b V^2
    Vpi = np.sqrt(-c[1] / c[0]) if c[0] < 0 else float('nan')
    print(f'f1^2 -> 0 extrapolation: V_PI ~ {Vpi:.1f} V')
print('DONE eigcheck')
