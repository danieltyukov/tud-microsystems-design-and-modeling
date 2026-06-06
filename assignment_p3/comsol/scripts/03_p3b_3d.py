"""Part 3b: 3D model, Part-2 dimensions, comb fingers removed.

 b.1 background: out-of-plane spring constant ratio (t/Wb)^2 (analytic).
 b.2: sweep thickness t -> eigenfrequencies + in-plane/out-of-plane character;
      find t where the fundamental mode is the in-plane actuation mode;
      at the final t: first 5 3D mode maps + static k (in-plane and
      out-of-plane) vs analytic.

Geometry (um, z = extrusion):
  bar 100x30 y[-30,0]; column 30x40 y[-70,-30]; plate 100x30 y[-100,-70];
  beams 75x5 (x to +-90, anchors at x=+-90), slots 1 um as in Part 2.

Run: ~/workspace/comsol-mcp/.venv/bin/python 03_p3b_3d.py
"""
import sys
sys.path.insert(0, '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p3/comsol/scripts')
import csv
import numpy as np
import mph
from jpype import JInt
from model_lib_p3 import MODELS, DATA, FIGS

T_SWEEP = [0.2, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0]  # um
T_FINAL = 6.0
E, RHO = 170e9, 2330.0
LB, WB = 75e-6, 5e-6

client = mph.start(cores=6)
model = client.create('phase_shifter_p3b_3d')
jm = model.java

model.parameter('t_si', '6[um]')
model.parameter('F_test', '1[uN]')

comp = jm.component().create('comp1', True)
geom = comp.geom().create('geom1', 3)
geom.lengthUnit('m')

wp = geom.create('wp1', 'WorkPlane')
wp.set('unite', True)
wpg = wp.geom()


def rect(tag, x, y, w, h):
    r = wpg.create(tag, 'Rectangle')
    r.set('pos', [f'{x}[um]', f'{y}[um]'])
    r.set('size', [f'{w}[um]', f'{h}[um]'])


rect('bar', -50, -30, 100, 30)
rect('col', -15, -70, 30, 40)
rect('mas', -50, -100, 100, 30)
rect('bUL', -90, -36, 75, 5)
rect('bUR', 15, -36, 75, 5)
rect('bLL', -90, -69, 75, 5)
rect('bLR', 15, -69, 75, 5)

ext = geom.create('ext1', 'Extrude')
ext.set('extrudefrom', 'workplane')
ext.selection('input').set('wp1')
ext.setIndex('distance', 't_si', 0)
geom.run('fin')
print('3D geometry built')

# ---- selections ----
from jpype import JInt as JI


def box3(tag, x1, x2, y1, y2, z1, z2, dim):
    s = comp.selection().create(tag, 'Box')
    s.set('entitydim', JI(dim))
    s.set('xmin', x1 * 1e-6); s.set('xmax', x2 * 1e-6)
    s.set('ymin', y1 * 1e-6); s.set('ymax', y2 * 1e-6)
    s.set('zmin', z1 * 1e-6); s.set('zmax', z2 * 1e-6)
    s.set('condition', 'inside')
    return s


box3('selAnchL', -90.01, -89.99, -70, -30, -0.5, 10.5, 2)
box3('selAnchR', 89.99, 90.01, -70, -30, -0.5, 10.5, 2)
u = comp.selection().create('selAnchors', 'Union')
u.set('entitydim', JI(2))
u.set('input', ['selAnchL', 'selAnchR'])
box3('selPlate', -50.01, 50.01, -100.01, -69.99, -0.5, 10.5, 3)

# ---- operators / variables ----
itg = comp.cpl().create('intAll', 'Integration')
itg.selection().all()
itg.set('opname', 'intAll')
avp = comp.cpl().create('avePlate', 'Average')
avp.selection().named('selPlate')
avp.set('opname', 'avePlate')
var = comp.variable().create('var1')
var.set('oopfrac', 'intAll(w^2)/intAll(u^2+v^2+w^2)')
var.set('uy_plate', 'avePlate(v)')
var.set('uz_plate', 'avePlate(w)')

# ---- physics / material ----
solid = comp.physics().create('solid', 'SolidMechanics', 'geom1')
fix = solid.create('fix1', 'Fixed', 2)
fix.selection().named('selAnchors')

mat = comp.material().create('matSi', 'Common')
mat.label('Si isotropic')
mat.propertyGroup('def').set('density', '2330[kg/m^3]')
enu = mat.propertyGroup().create('Enu', "Young's modulus and Poisson's ratio")
enu.set('E', '170[GPa]')
enu.set('nu', '0.28')

# ---- swept mesh (thin extruded layers) ----
mesh = comp.mesh().create('mesh1')
ftri = mesh.create('ftri1', 'FreeTri')
# source: bottom faces (z=0)
box3('selBot', -101, 101, -101, 1, -0.01, 0.01, 2)
ftri.selection().named('selBot')
sz = mesh.feature('size')
sz.set('custom', True)
sz.set('hmax', '4[um]')
sz.set('hmin', '0.05[um]')
szb = mesh.create('szBeam', 'Size')
szb.selection().geom('geom1', 2)
box3('selBotBeams', -91, 91, -70, -30, -0.01, 0.01, 2)
szb.selection().named('selBotBeams')
szb.set('custom', True)
szb.set('hmaxactive', True)
szb.set('hmax', '1.5[um]')
swe = mesh.create('swe1', 'Sweep')
dis = swe.create('dis1', 'Distribution')
dis.set('numelem', JI(2))
mesh.run()
print('mesh elements:', mesh.stat().getNumElem())

# ---- eigenfrequency study ----
std = jm.study().create('stdEig')
eig = std.create('eig', 'Eigenfrequency')
eig.set('neigsactive', True)
eig.set('neigs', JInt(8))
eig.set('shift', '1[kHz]')

gev = None
dstag = None
rows = []
for t_um in T_SWEEP:
    model.parameter('t_si', f'{t_um}[um]')
    jm.study('stdEig').run()
    if gev is None:
        for tg in [str(x) for x in jm.sol().tags()]:
            try:
                if str(jm.sol(tg).study()) == 'stdEig':
                    soltag = tg
            except Exception:
                pass
        ds = jm.result().dataset().create('dsEig', 'Solution')
        ds.set('solution', soltag)
        gev = jm.result().numerical().create('gevF', 'Global')
        gev.set('data', 'dsEig')
    gev.set('expr', ['freq'])
    freqs = np.array(gev.getReal()[0], dtype=float)
    gev.set('expr', ['oopfrac'])
    oop = np.array(gev.getReal()[0], dtype=float)
    rows.append([t_um] + [x for pair in zip(freqs[:6], oop[:6]) for x in pair])
    lab = '  '.join(f'{f/1e3:7.1f}k({o:.2f})' for f, o in zip(freqs[:6], oop[:6]))
    print(f't={t_um:4.1f} um: {lab}')

hdr = 't_um,' + ','.join(f'f{i}_Hz,oop{i}' for i in range(1, 7))
np.savetxt(f'{DATA}/p3b_modes_vs_t.csv', np.array(rows), delimiter=',',
           header=hdr, comments='')

# ---- final t: mode maps ----
model.parameter('t_si', f'{T_FINAL}[um]')
jm.study('stdEig').run()
gev.set('expr', ['freq'])
freqs = np.array(gev.getReal()[0], dtype=float)
gev.set('expr', ['oopfrac'])
oop = np.array(gev.getReal()[0], dtype=float)
print('final t modes:', ', '.join(f'{f/1e3:.1f} kHz (oop {o:.2f})'
                                  for f, o in zip(freqs[:5], oop[:5])))

pg = jm.result().create('pgM', 'PlotGroup3D')
pg.set('data', 'dsEig')
s1 = pg.create('s1', 'Surface')
s1.set('expr', 'solid.disp')
s1.set('unit', 'nm')
dfm = s1.create('def1', 'Deform')
dfm.set('expr', ['u', 'v', 'w'])
dfm.set('scaleactive', False)
img = jm.result().export().create('imgM', 'Image')
img.set('sourceobject', 'pgM')
img.set('size', 'manualweb')
img.set('unit', 'px')
img.set('width', '1400')
img.set('height', '900')
for i in range(5):
    pg.set('looplevel', [str(i + 1)])
    pg.run()
    fn = f'{FIGS}/p3b_mode{i+1}_{freqs[i]/1e3:.0f}kHz.png'
    img.set('pngfilename', fn)
    img.run()
    print('exported', fn)

# ---- static spring constants at final t ----
bl = solid.create('blK', 'BodyLoad', 3)
bl.selection().named('selPlate')
std2 = jm.study().create('stdK')
std2.create('stat', 'Stationary')
ds2 = None
gev2 = None


def solve_k(direction):
    global ds2, gev2
    fexpr = ['0', '0', '0']
    fexpr[direction] = 'F_test/(100[um]*30[um]*t_si)'
    bl.set('forceReferenceVolume', fexpr)
    jm.study('stdK').run()
    for tg in [str(x) for x in jm.sol().tags()]:
        try:
            if str(jm.sol(tg).study()) == 'stdK':
                soltag2 = tg
        except Exception:
            pass
    if ds2 is None:
        ds2 = jm.result().dataset().create('dsK', 'Solution')
        gev2 = jm.result().numerical().create('gevK', 'Global')
        gev2.set('data', 'dsK')
    ds2.set('solution', soltag2)
    gev2.set('expr', ['uy_plate' if direction == 1 else 'uz_plate'])
    val = float(gev2.getReal()[0][0])
    return 1e-6 / abs(val)


k_ip = solve_k(1)
k_oop = solve_k(2)
t_m = T_FINAL * 1e-6
k_ip_ana = 4 * E * t_m * WB**3 / LB**3
k_oop_ana = 4 * E * t_m**3 * WB / LB**3
print(f'k in-plane : FEM {k_ip:8.2f} N/m, analytic {k_ip_ana:8.2f} N/m')
print(f'k out-plane: FEM {k_oop:8.2f} N/m, analytic {k_oop_ana:8.2f} N/m')
with open(f'{DATA}/p3b_spring_constants.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['t_um', 'k_ip_FEM', 'k_ip_analytic', 'k_oop_FEM',
                'k_oop_analytic'])
    w.writerow([T_FINAL, k_ip, k_ip_ana, k_oop, k_oop_ana])

model.parameter('t_si', f'{T_FINAL}[um]')
model.save(f'{MODELS}/phase_shifter_p3b_3d.mph')
print('DONE Part 3b 3D')
