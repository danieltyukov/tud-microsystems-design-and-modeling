"""Prestressed eigenfrequencies vs bias, robust version.

For each target bias the stationary step runs a continuation ramp
0 -> Vt (50 V steps, continuation from last parameter), so the eigenfrequency
step linearizes about a state reached on the physical branch -- no cold-start
failures (Part-2 gotcha: cold high-V EM solves diverge or converge to garbage).

Also re-exports the potential field at the operating point and saves the EM
model with solutions cleared.

Run: ~/workspace/comsol-mcp/.venv/bin/python 06_eigcheck_warm.py
"""
import sys
sys.path.insert(0, '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p3/comsol/scripts')
import numpy as np
import mph
from jpype import JInt
from model_lib_p3 import MODELS, DATA, FIGS

TARGETS = [0, 100, 200, 300, 400, 500, 600, 627, 700, 760]

client = mph.start(cores=6)
model = client.load(f'{MODELS}/phase_shifter_p3_em.mph')
jm = model.java
print('EM model loaded')

std = jm.study().create('stdPre2')
stat = std.create('stat', 'Stationary')
stat.set('useparam', 'on')
stat.set('pname', ['Vd'])
stat.set('punit', ['V'])
stat.set('sweeptype', 'sparse')
stat.set('pcontinuationmode', 'last')
eig = std.create('eig', 'Eigenfrequency')
eig.set('neigsactive', True)
eig.set('neigs', JInt(5))
eig.set('shift', '1e4')

rows = []
gev = None
dstag_ok = None
for Vt in TARGETS:
    plist = f'range(0,50,{Vt}) {Vt}' if Vt > 0 else '0'
    stat.set('plistarr', [plist])
    try:
        std.run()
    except Exception as e:
        print(f'Vt={Vt}: failed ({str(e)[:90]})')
        break
    if gev is None:
        gev = jm.result().numerical().create('gevPre2', 'Global')
        gev.set('expr', ['freq'])
        for t in [str(x) for x in jm.result().dataset().tags()]:
            try:
                gev.set('data', t)
                arr = np.array(gev.getReal()[0], dtype=float)
                if len(arr) >= 5 and arr[0] > 1e4:
                    dstag_ok = t
                    break
            except Exception:
                continue
        print('eigen dataset:', dstag_ok)
    gev.set('data', dstag_ok)
    gev.set('expr', ['freq'])
    freqs = np.array(gev.getReal()[0], dtype=float)
    rows.append([Vt] + list(freqs[:5]))
    print(f'Vd={Vt:4d} V: f = ' + '  '.join(f'{f/1e3:9.2f}k' for f in freqs[:5]))

arr = np.array(rows)
np.savetxt(f'{DATA}/p3a_eig_vs_bias.csv', arr, delimiter=',',
           header='Vd_V,f1_Hz,f2_Hz,f3_Hz,f4_Hz,f5_Hz', comments='')
if len(arr) >= 4:
    V = arr[:, 0]
    for col, nm in [(1, 'f1'), (2, 'f2')]:
        fsq = arr[:, col] ** 2
        c = np.polyfit(V[-4:] ** 2, fsq[-4:], 1)
        if c[0] < 0:
            print(f'{nm}^2 -> 0 extrapolation: V_PI ~ {np.sqrt(-c[1]/c[0]):.1f} V')

# softened mode-2 shape at the last bias
try:
    pg2 = jm.result().create('pgEig2', 'PlotGroup2D')
    pg2.set('data', dstag_ok)
    pg2.set('looplevel', ['2'])
    s1 = pg2.create('s1', 'Surface')
    s1.set('expr', 'solid.disp')
    dfm = s1.create('def1', 'Deform')
    dfm.set('expr', ['u', 'v'])
    pg2.run()
    img2 = jm.result().export().create('imgEig2', 'Image')
    img2.set('sourceobject', 'pgEig2')
    img2.set('pngfilename', f'{FIGS}/p3a_softening_mode2.png')
    img2.set('size', 'manualweb')
    img2.set('unit', 'px')
    img2.set('height', '800')
    img2.set('width', '1400')
    img2.run()
    print('mode-2 image exported')
except Exception as e:
    print('mode-2 image failed:', str(e)[:100])

# potential field at the operating point (fresh single ramp to 627 V)
try:
    std2 = jm.study().create('stdOp')
    stat2 = std2.create('stat', 'Stationary')
    stat2.set('useparam', 'on')
    stat2.set('pname', ['Vd'])
    stat2.set('punit', ['V'])
    stat2.set('plistarr', ['range(0,50,600) 627'])
    stat2.set('sweeptype', 'sparse')
    stat2.set('pcontinuationmode', 'last')
    std2.run()
    soltag = None
    for t in [str(x) for x in jm.sol().tags()]:
        try:
            if str(jm.sol(t).study()) == 'stdOp':
                soltag = t
        except Exception:
            pass
    ds = jm.result().dataset().create('dsOp', 'Solution')
    ds.set('solution', soltag)
    pg = jm.result().create('pgPot2', 'PlotGroup2D')
    pg.set('data', 'dsOp')
    surf = pg.create('s1', 'Surface')
    surf.set('expr', 'V')
    pg.run()
    img = jm.result().export().create('imgPot2', 'Image')
    img.set('sourceobject', 'pgPot2')
    img.set('pngfilename', f'{FIGS}/p3a_potential_op.png')
    img.set('size', 'manualweb')
    img.set('unit', 'px')
    img.set('height', '800')
    img.set('width', '1400')
    img.run()
    print('potential image exported')
except Exception as e:
    print('potential image failed:', str(e)[:140])

# clean + save
for t in [str(x) for x in jm.sol().tags()]:
    try:
        jm.sol(t).clearSolutionData()
    except Exception:
        pass
model.parameter('Vd', '0[V]')
model.save(f'{MODELS}/phase_shifter_p3_em.mph')
print('DONE eigcheck warm')
