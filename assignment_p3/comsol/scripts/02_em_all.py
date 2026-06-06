"""Part 3a.4 (iii),(iv) + pull-in verification: EM model studies.

 1. x(V) continuation sweep 0 -> 760 V (validates V_DC_MAX operating point,
    x(627 V) ~ 1/(2b) = 1.667 um, and brackets pull-in by divergence).
    Also records terminal capacitance C11(V).
 2. prestressed eigenfrequencies vs bias (warm restart): which mode softens,
    f^2 -> 0 extrapolation = V_PI and its character (axial vs side).
 3. max DC acceleration before bar-waveguide contact at Vd = 0.

Run: ~/workspace/comsol-mcp/.venv/bin/python 02_em_all.py
"""
import sys
sys.path.insert(0, '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p3/comsol/scripts')
import csv
import numpy as np
import mph
from jpype import JInt
import model_lib_p3 as lib
from model_lib_p3 import MODELS, DATA, FIGS

G0 = 200e-9
V_OP = 627.0

client = mph.start(cores=6)
model = client.load(f'{MODELS}/phase_shifter_p3_em.mph')
jm = model.java
comp = jm.component('comp1')
print('EM model loaded')

# ---------------- 1) x(V) continuation sweep -----------------
std = jm.study().create('stdVs')
stat = std.create('stat', 'Stationary')
stat.set('useparam', 'on')
stat.set('pname', ['Vd'])
stat.set('plistarr', ['range(0,10,620) 627 range(640,10,760)'])
stat.set('punit', ['V'])
stat.set('sweeptype', 'sparse')
stat.set('pcontinuationmode', 'last')
try:
    std.run()
    print('sweep completed full range (no pull-in below 760 V)')
except Exception as e:
    print('sweep stopped (divergence expected):', str(e)[:150])

soltag = None
for t in [str(x) for x in jm.sol().tags()]:
    try:
        if str(jm.sol(t).study()) == 'stdVs':
            soltag = t
    except Exception:
        pass
ds = jm.result().dataset().create('dsVs', 'Solution')
ds.set('solution', soltag)
gev = jm.result().numerical().create('gevVs', 'Global')
gev.set('data', 'dsVs')


def geval(expr):
    gev.set('expr', [expr])
    return np.array(gev.getReal()[0], dtype=float)


V = geval('Vd')
uy = geval('uy_mass')
uyb = geval('uy_bar')
C11 = geval('es.C11')
np.savetxt(f'{DATA}/p3a_q4iii_vsweep.csv',
           np.column_stack([V, uy, uyb, C11]), delimiter=',',
           header='Vd_V,uy_mass_m,uy_bar_m,C11_F', comments='')
for i in range(0, len(V), max(1, len(V) // 25)):
    print(f'V={V[i]:6.1f}: uy={uy[i]*1e9:10.2f} nm  C={C11[i]*1e15:8.4f} fF')
print(f'last converged: V={V[-1]:.1f} V, uy={uy[-1]*1e9:.1f} nm')
iop = np.argmin(np.abs(V - V_OP))
if abs(V[iop] - V_OP) < 1:
    print(f'x({V_OP:.0f} V) FEM = {abs(uy[iop])*1e6:.4f} um '
          f'(analytic 1/(2b) = 1.667 um)')

# potential field image near the operating point
try:
    lvl = int(iop + 1)
    pg = jm.result().create('pgPot', 'PlotGroup2D')
    pg.set('data', 'dsVs')
    pg.set('looplevel', [str(lvl)])
    surf = pg.create('s1', 'Surface')
    surf.set('expr', 'es.V')
    pg.run()
    img = jm.result().export().create('imgPot', 'Image')
    img.set('sourceobject', 'pgPot')
    img.set('pngfilename', f'{FIGS}/p3a_potential_op.png')
    img.set('size', 'manualweb')
    img.set('unit', 'px')
    img.set('height', '800')
    img.set('width', '1400')
    img.run()
    print('potential image exported')
except Exception as e:
    print('potential image failed:', str(e)[:120])

# ---------------- 2) prestressed eigenfrequencies vs bias -----------------
std2 = jm.study().create('stdPre')
stat2 = std2.create('stat', 'Stationary')
eig2 = std2.create('eig', 'Eigenfrequency')
eig2.set('neigsactive', True)
eig2.set('neigs', JInt(5))
eig2.set('shift', '1e4')

rows = []
gev2 = None
dstag_ok = None
warm = False
for Vd in [0, 100, 200, 300, 400, 450, 500, 550, 600, 627, 650, 675, 700, 720]:
    model.parameter('Vd', f'{Vd}[V]')
    try:
        std2.run()
    except Exception as e:
        print(f'Vd={Vd}: solve failed ({str(e)[:90]})')
        break
    if not warm:
        # warm restart: continue each next bias from the previous solution
        s2tag = None
        for t in [str(x) for x in jm.sol().tags()]:
            try:
                if str(jm.sol(t).study()) == 'stdPre':
                    s2tag = t
            except Exception:
                pass
        try:
            jm.sol(s2tag).feature('v1').set('initmethod', 'sol')
            jm.sol(s2tag).feature('v1').set('initsol', s2tag)
            warm = True
            print('warm restart enabled on', s2tag)
        except Exception as e:
            print('warm restart setup failed:', str(e)[:100])
    if gev2 is None:
        gev2 = jm.result().numerical().create('gevPre', 'Global')
        gev2.set('expr', ['freq'])
        for t in [str(x) for x in jm.result().dataset().tags()]:
            try:
                gev2.set('data', t)
                arrtest = np.array(gev2.getReal()[0], dtype=float)
                if len(arrtest) >= 5:
                    dstag_ok = t
                    break
            except Exception:
                continue
        print('eigen dataset:', dstag_ok)
    gev2.set('data', dstag_ok)
    gev2.set('expr', ['freq'])
    freqs = np.array(gev2.getReal()[0], dtype=float)
    rows.append([Vd] + list(freqs[:5]))
    print(f'Vd={Vd:4d} V: f = ' + '  '.join(f'{f/1e3:9.2f}k' for f in freqs[:5]))

arr = np.array(rows)
np.savetxt(f'{DATA}/p3a_eig_vs_bias.csv', arr, delimiter=',',
           header='Vd_V,f1_Hz,f2_Hz,f3_Hz,f4_Hz,f5_Hz', comments='')

# critical-mode image at highest converged bias
try:
    pg2 = jm.result().create('pgEig', 'PlotGroup2D')
    pg2.set('data', dstag_ok)
    pg2.set('looplevel', ['1'])
    s1 = pg2.create('s1', 'Surface')
    s1.set('expr', 'solid.disp')
    dfm = s1.create('def1', 'Deform')
    dfm.set('expr', ['u', 'v'])
    dfm.set('scaleactive', True)
    dfm.set('scale', '2e6')
    pg2.run()
    img2 = jm.result().export().create('imgEig', 'Image')
    img2.set('sourceobject', 'pgEig')
    img2.set('pngfilename', f'{FIGS}/p3a_critical_mode.png')
    img2.set('size', 'manualweb')
    img2.set('unit', 'px')
    img2.set('height', '800')
    img2.set('width', '1400')
    img2.run()
    print('critical mode image exported')
except Exception as e:
    print('eig image failed:', str(e)[:120])

# f1^2 -> 0 extrapolation
if len(arr) >= 4:
    Vv = arr[:, 0]
    for col, nm in [(1, 'f1'), (2, 'f2')]:
        fsq = arr[:, col] ** 2
        c = np.polyfit(Vv[-4:] ** 2, fsq[-4:], 1)
        if c[0] < 0:
            print(f'{nm}^2 -> 0 extrapolation: V_PI ~ {np.sqrt(-c[1]/c[0]):.1f} V')

# ---------------- 3) max DC acceleration (Vd = 0) -----------------
model.parameter('Vd', '0[V]')
solid = comp.physics('solid')
bl = solid.create('blAcc', 'BodyLoad', 2)
bl.selection().named('selSolids')
bl.set('forceReferenceVolume', ['0', 'solid.rho*a_acc', '0'])

std3 = jm.study().create('stdAcc')
std3.create('stat', 'Stationary')
ds3 = None
gev3 = None


def solve_uy_bar(a):
    global ds3, gev3
    model.parameter('a_acc', f'{a}[m/s^2]')
    jm.study('stdAcc').run()
    soltag3 = None
    for t in [str(x) for x in jm.sol().tags()]:
        try:
            if str(jm.sol(t).study()) == 'stdAcc':
                soltag3 = t
        except Exception:
            pass
    if ds3 is None:
        ds3 = jm.result().dataset().create('dsAcc', 'Solution')
        gev3 = jm.result().numerical().create('gevAcc', 'Global')
        gev3.set('data', 'dsAcc')
    ds3.set('solution', soltag3)
    gev3.set('expr', ['uy_bar'])
    return float(gev3.getReal()[0][0])


uy0 = solve_uy_bar(0.0)
a_test = 1.0e6
uy1 = solve_uy_bar(a_test)
a_max = a_test * (G0 - uy0) / (uy1 - uy0)
uyv = solve_uy_bar(a_max)
resid = (uyv - G0) / G0
print(f'a_max(Vd=0) = {a_max:.4e} m/s^2 ({a_max/9.81:,.0f} g), '
      f'verify uy_bar = {uyv*1e9:.2f} nm (resid {resid*100:.3f} %)')
with open(f'{DATA}/p3a_q4iv_max_acceleration.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['Vd_V', 'uy_rest_m', 'a_max_m_s2', 'a_max_g', 'verify_resid'])
    w.writerow([0, uy0, a_max, a_max / 9.81, resid])

# reset swept parameters before saving (P2 contamination gotcha)
model.parameter('a_acc', '0[m/s^2]')
model.parameter('Vd', '0[V]')
bl.active(False)
model.save(f'{MODELS}/phase_shifter_p3_em.mph')
print('DONE EM (iii),(iv) + pull-in verification')
