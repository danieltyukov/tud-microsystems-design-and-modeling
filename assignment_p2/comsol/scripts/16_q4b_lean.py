"""Q4b lean re-run (monitored): prestressed-eigenvalue softening on a finer mesh
with electrostatic thickness = 1, per the teaching-team hints. Focused voltage
points near the instability; compares against the existing June-6 coarse data
(which used ES thickness = 200 nm) to expose the spurious-fringing effect the
teaching team described. Writes q4b_eig_vs_bias_fine.csv + mode image.

Unbuffered (run with python -u) for live monitoring.
"""
import sys
sys.path.insert(0, '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p2/comsol/scripts')
import mph
import numpy as np
from jpype import JInt
from model_lib import MODELS, DATA, FIGS   # model_lib sets mph.option('classkit', True)

MSC = 0.7
VLIST = [0, 100, 150, 180, 200, 210, 215, 218, 220, 222]

print('starting COMSOL (classkit)...', flush=True)
client = mph.start(cores=6)
print('COMSOL up:', client.version, flush=True)
model = client.load(f'{MODELS}/phase_shifter_em.mph')
print('model loaded', flush=True)
jm = model.java
comp = jm.component('comp1')
es = comp.physics('es')

# teaching-team hint: ES out-of-plane thickness = 1 (200 nm pulls in spurious
# fringing in 6.4); mechanical thickness stays t_si.
es.prop('d').set('d', '1')
print('es thickness = 1', flush=True)

model.parameter('msc', f'{MSC}')
comp.mesh('mesh1').run()
ne = comp.mesh('mesh1').getNumElem()
print(f'fine mesh msc={MSC}: {ne} elements', flush=True)

std = jm.study().create('stdEigF')
std.create('stat', 'Stationary')
eig = std.create('eig', 'Eigenfrequency')
eig.set('neigs', JInt(4))
eig.set('shift', '1e4')

rows, gev, dstag = [], None, None
for Vd in VLIST:
    model.parameter('Vd', f'{Vd}[V]')
    try:
        std.run()
    except Exception as e:
        print(f'  Vd={Vd}: solve failed ({str(e)[:70]})', flush=True)
        continue
    if gev is None:
        gev = jm.result().numerical().create('gevEF', 'Global')
        gev.set('expr', ['freq'])
        for t in [str(x) for x in jm.result().dataset().tags()]:
            try:
                gev.set('data', t); np.array(gev.getReal()[0], dtype=float)
                dstag = t; break
            except Exception:
                continue
        print('eigen dataset:', dstag, flush=True)
    gev.set('data', dstag); gev.set('expr', ['freq'])
    f = np.array(gev.getReal()[0], dtype=float)
    rows.append([Vd] + list(f[:4]))
    print(f'  Vd={Vd:4d} V: ' + '  '.join(f'{x/1e3:8.1f} kHz' for x in f[:4]), flush=True)

arr = np.array(rows)
np.savetxt(f'{DATA}/q4b_eig_vs_bias_fine.csv', arr, delimiter=',',
           header='Vd_V,f1_Hz,f2_Hz,f3_Hz,f4_Hz', comments='')
print(f'wrote q4b_eig_vs_bias_fine.csv ({len(arr)} rows)', flush=True)


def vinst(col):
    V = arr[:, 0]; y = arr[:, col]**2
    m = V >= V[-4]
    sl, it = np.polyfit(V[m]**2, y[m], 1)
    return np.sqrt(-it/sl) if sl < 0 else float('nan')


print(f'\nFINE (thickness=1): V_inst axial={vinst(1):.0f} V  '
      f'rocking={vinst(2):.0f} V  finger={vinst(3):.0f} V', flush=True)

# compare with existing coarse 200nm data
try:
    old = np.genfromtxt(f'{DATA}/q4b_eig_vs_bias.csv', delimiter=',', names=True)
    Vo = old['Vd_V']
    def vinst_old(key):
        y = old[key]**2; m = Vo >= Vo[-4]
        sl, it = np.polyfit(Vo[m]**2, y[m], 1)
        return np.sqrt(-it/sl) if sl < 0 else float('nan')
    print(f'COARSE (thickness=200nm, Jun-6): V_inst axial={vinst_old("f1_Hz"):.0f} V  '
          f'rocking={vinst_old("f2_Hz"):.0f} V  finger={vinst_old("f3_Hz"):.0f} V', flush=True)
    print(f'f2(0V): fine={arr[0,2]/1e3:.0f} kHz vs coarse={old["f2_Hz"][0]/1e3:.0f} kHz', flush=True)
except Exception as e:
    print('coarse compare failed:', str(e)[:80], flush=True)

# critical (rocking) mode shape at last converged bias
try:
    pg = jm.result().create('pgEF', 'PlotGroup2D')
    pg.set('data', dstag); pg.set('looplevel', ['2'])
    surf = pg.create('s1', 'Surface'); surf.set('expr', 'solid.disp')
    dfm = surf.create('def1', 'Deform'); dfm.set('expr', ['u', 'v'])
    dfm.set('scaleactive', True); dfm.set('scale', '5e6')
    img = jm.result().export().create('imgEF', 'Image')
    img.set('sourceobject', 'pgEF')
    img.set('pngfilename', f'{FIGS}/q4b_critical_mode.png')
    img.set('size', 'manualweb'); img.set('unit', 'px')
    img.set('height', '720'); img.set('width', '1280'); img.run()
    print('mode image exported', flush=True)
except Exception as e:
    print('mode image failed:', str(e)[:100], flush=True)

print('DONE q4b lean', flush=True)
