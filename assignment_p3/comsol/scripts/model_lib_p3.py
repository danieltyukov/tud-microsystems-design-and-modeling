"""Shared builder for the ET4260 Part-3a optimized comb-drive phase shifter.

Geometry (2D top view, um, y=0 at bar top, motion = y, waveguide above bar):
  bar    : x in [-50, 50],  y in [-10,   0]   (L=100, Wm=10)
  column : x in [ -5,  5],  y in [-40, -10]   (Wm=10 wide, Lm=30 long)
  spine  : x in [-50, 50],  y in [-50, -40]
  beams  : Lb=40, Wb=1.85, 1 um slots from the plates, anchors at x=+-45
           upper pair y in [-12.85, -11], lower pair y in [-39, -37.15]
  movers : 29 fingers, Wc=0.7, Lc=5.5, y in [-55.5, -50], pitch 3.2/1.6 um
  fixed  : 28 inner (0.7) + 2 ends (3Wc=2.1), y in [-60, -54.5] (d0=1 overlap)
  air    : x in [-49.5, 49.5], y in [-60, -50] (stator floor at y=-60)

Mirrors assignment_p2/comsol/scripts/model_lib.py (verified recipe).
"""
import mph
from jpype import JInt

UM = 1e-6
BASE = '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p3'
MODELS = f'{BASE}/comsol/models'
DATA = f'{BASE}/data'
FIGS = f'{BASE}/figures'

# ----- final Part-3a dimensions (um) -----
P = dict(
    L=100, Wm=10, Lm=30, Lb=40, Wb=1.85,
    Lc=5.5, Wc=0.7, h=0.9, d0=1.0, g0=0.2,
    t_si=0.2, N=58,
)
MAT = dict(E='170[GPa]', nu='0.28', rho='2330[kg/m^3]')

N_MOV = 29           # moving fingers (N=58 gaps)
N_FIX_IN = 28        # inner fixed fingers (Wc)
WC_END = 3 * P['Wc']  # end fixed fingers 3*Wc = 2.1 um

Y_SPINE_BOT = -50.0
Y_MOV_TIP = Y_SPINE_BOT - P['Lc']            # -55.5
Y_FIX_TIP = Y_MOV_TIP + P['d0']              # -54.5 (d0 overlap)
Y_STATOR = Y_FIX_TIP - P['Lc']               # -60.0 (fixed fingers same Lc)
X_AIR = 49.5

ANCH_X = 45.0        # beam anchor |x|
BEAM_UP_Y = -12.85   # upper beam pair lower edge offset handled below
BEAM_LO_Y = -39.0


def comb_layout():
    """Return (moving, fixed) lists of (x_left, width) in um."""
    wc, h = P['Wc'], P['h']
    width = 2 * WC_END + N_MOV * wc + N_FIX_IN * wc + P['N'] * h
    x = -width / 2
    fixed = [(x, WC_END)]
    x += WC_END + h
    moving = []
    for i in range(N_MOV):
        moving.append((x, wc))
        x += wc + h
        if i < N_MOV - 1:
            fixed.append((x, wc))
            x += wc + h
    fixed.append((x, WC_END))
    assert abs((x + WC_END) - width / 2) < 1e-9, f'comb end mismatch: {x + WC_END}'
    return moving, fixed


def rect(geom, tag, x, y, w, hgt, contribute=None):
    r = geom.create(tag, 'Rectangle')
    r.set('pos', [f'{x}[um]', f'{y}[um]'])
    r.set('size', [f'{w}[um]', f'{hgt}[um]'])
    if contribute:
        r.set('contributeto', contribute)
    return r


def build_geometry(model, with_air):
    jm = model.java
    comp = jm.component().create('comp1', True)
    geom = comp.geom().create('geom1', 2)
    geom.lengthUnit('m')

    for cs in ['cselShuttle', 'cselMass', 'cselBar', 'cselBeams',
               'cselFingers', 'cselColumn'] + (['cselAir'] if with_air else []):
        geom.selection().create(cs, 'CumulativeSelection')

    rect(geom, 'bar', -50, -10, 100, 10, 'cselBar')
    rect(geom, 'col', -5, -40, 10, 30, 'cselColumn')
    rect(geom, 'mas', -50, -50, 100, 10, 'cselMass')
    # beams attach to the column flanks; 1 um release slots from the plates
    rect(geom, 'bUL', -ANCH_X, -12.85, 40, 1.85, 'cselBeams')
    rect(geom, 'bUR', 5, -12.85, 40, 1.85, 'cselBeams')
    rect(geom, 'bLL', -ANCH_X, -39.0, 40, 1.85, 'cselBeams')
    rect(geom, 'bLR', 5, -39.0, 40, 1.85, 'cselBeams')

    moving, fixed = comb_layout()
    for i, (xl, w) in enumerate(moving):
        rect(geom, f'mf{i}', xl, Y_MOV_TIP, w, P['Lc'], 'cselFingers')

    if with_air:
        rect(geom, 'airR', -X_AIR, Y_STATOR, 2 * X_AIR, -Y_STATOR + Y_SPINE_BOT)
        fix_tags = []
        for i, (xl, w) in enumerate(fixed):
            tag = f'ff{i}'
            rect(geom, tag, xl, Y_STATOR - 1.0, w, P['Lc'] + 1.0)
            fix_tags.append(tag)
        dif = geom.create('difAir', 'Difference')
        dif.selection('input').set('airR')
        dif.selection('input2').set(*fix_tags)
        dif.set('contributeto', 'cselAir')
        dif.set('keep', False)
    geom.run('fin')
    return geom


def add_component_selections(model, with_air):
    """Box selections at component level (meters)."""
    jm = model.java
    comp = jm.component('comp1')

    def box(tag, x1, x2, y1, y2, dim):
        s = comp.selection().create(tag, 'Box')
        s.set('entitydim', JInt(dim))
        s.set('xmin', x1 * UM); s.set('xmax', x2 * UM)
        s.set('ymin', y1 * UM); s.set('ymax', y2 * UM)
        s.set('condition', 'inside')
        return s

    # anchors: beam outer end edges at x=+-45
    box('selAnchorL', -ANCH_X - 0.01, -ANCH_X + 0.01, -40, -10, 1)
    box('selAnchorR', ANCH_X - 0.01, ANCH_X + 0.01, -40, -10, 1)
    u = comp.selection().create('selAnchors', 'Union')
    u.set('entitydim', JInt(1))
    u.set('input', ['selAnchorL', 'selAnchorR'])
    # bar top edge (waveguide-facing)
    box('selBarTop', -50.01, 50.01, -0.01, 0.01, 1)

    if with_air:
        _, fixed = comb_layout()
        tags = []
        s = box('selStatLine', -X_AIR - 0.01, X_AIR + 0.01,
                Y_STATOR - 0.01, Y_STATOR + 0.01, 1)
        tags.append('selStatLine')
        for i, (xl, w) in enumerate(fixed):
            box(f'selFF{i}', xl - 0.2, xl + w + 0.2,
                Y_STATOR - 0.21, Y_FIX_TIP + 0.21, 1)
            tags.append(f'selFF{i}')
        u = comp.selection().create('selTerminal', 'Union')
        u.set('entitydim', JInt(1))
        u.set('input', tags)
        # air side openings (moving-mesh symmetry)
        box('selAirSideL', -X_AIR - 0.01, -X_AIR + 0.01,
            Y_STATOR - 0.01, Y_SPINE_BOT + 0.01, 1)
        box('selAirSideR', X_AIR - 0.01, X_AIR + 0.01,
            Y_STATOR - 0.01, Y_SPINE_BOT + 0.01, 1)
        u2 = comp.selection().create('selAirSidesU', 'Union')
        u2.set('entitydim', JInt(1))
        u2.set('input', ['selAirSideL', 'selAirSideR'])


def set_parameters(model, extra=None):
    pars = {
        'L': '100[um]', 'Wm': '10[um]', 'Lm': '30[um]', 'Lb': '40[um]',
        'Wb': '1.85[um]', 'Lc': '5.5[um]', 'Wc': '0.7[um]', 'h_gap': '0.9[um]',
        'd0': '1[um]', 'g0': '200[nm]', 't_si': '200[nm]',
        'F_test': '1[uN]',
        'Q_target': '80.7',           # analytic damping model at f_r
        'f1_fem': '517[kHz]',         # placeholder, updated after modes study
        'beta_dk': '1/(2*pi*f1_fem*Q_target)',
        'alpha_dm': '0[1/s]',
        'Vd': '0[V]',
        'a_acc': '0[m/s^2]',
        'msc': '1',
    }
    if extra:
        pars.update(extra)
    for k, v in pars.items():
        model.parameter(k, v)


def add_solids_union(model):
    jm = model.java
    comp = jm.component('comp1')
    u = comp.selection().create('selSolids', 'Union')
    u.set('entitydim', JInt(2))
    u.set('input', ['geom1_cselBar_dom', 'geom1_cselColumn_dom',
                    'geom1_cselMass_dom', 'geom1_cselBeams_dom',
                    'geom1_cselFingers_dom'])


def add_operators(model):
    jm = model.java
    comp = jm.component('comp1')
    ave = comp.cpl().create('aveMass', 'Average')
    ave.selection().named('geom1_cselMass_dom')
    ave.set('opname', 'aveMass')
    aveb = comp.cpl().create('aveBarTop', 'Average')
    aveb.selection().geom('geom1', 1)
    aveb.selection().named('selBarTop')
    aveb.set('opname', 'aveBarTop')
    itg = comp.cpl().create('intShuttle', 'Integration')
    itg.selection().named('selSolids')
    itg.set('opname', 'intShuttle')
    var = comp.variable().create('var1')
    var.set('uy_mass', 'aveMass(v)')
    var.set('uy_bar', 'aveBarTop(v)')
    var.set('m_shuttle', 'intShuttle(solid.rho*t_si)')


def add_solid(model):
    jm = model.java
    comp = jm.component('comp1')
    solid = comp.physics().create('solid', 'SolidMechanics', 'geom1')
    solid.selection().named('selSolids')
    solid.prop('Type2D').set('Type2D', 'PlaneStress')
    solid.prop('d').set('d', 't_si')
    fix = solid.create('fix1', 'Fixed', 1)
    fix.selection().named('selAnchors')
    return solid


def add_material(model, with_air):
    jm = model.java
    comp = jm.component('comp1')
    mat = comp.material().create('matSi', 'Common')
    mat.label('Si isotropic')
    mat.propertyGroup('def').set('density', MAT['rho'])
    enu = mat.propertyGroup().create('Enu', "Young's modulus and Poisson's ratio")
    enu.set('E', MAT['E'])
    enu.set('nu', MAT['nu'])
    mat.selection().named('selSolids')
    if with_air:
        mat.propertyGroup('def').set('relpermittivity',
                                     ['11.7', '0', '0', '0', '11.7', '0', '0', '0', '11.7'])
        air = comp.material().create('matAir', 'Common')
        air.label('Air')
        air.propertyGroup('def').set('relpermittivity',
                                     ['1', '0', '0', '0', '1', '0', '0', '0', '1'])
        air.propertyGroup('def').set('density', '1.2[kg/m^3]')
        air.selection().named('selAirOnly')
    return mat


def add_mesh(model, with_air):
    jm = model.java
    comp = jm.component('comp1')
    mesh = comp.mesh().create('mesh1')
    sz = mesh.feature('size')
    sz.set('custom', True)
    sz.set('hmax', '2.5[um]*msc')
    sz.set('hmin', '0.02[um]')
    szf = mesh.create('sizeFine', 'Size')
    if with_air:
        u = comp.selection().create('selFine', 'Union')
        u.set('entitydim', JInt(2))
        u.set('input', ['geom1_cselFingers_dom', 'geom1_cselAir_dom'])
        szf.selection().geom('geom1', 2)
        szf.selection().named('selFine')
    else:
        szf.selection().named('geom1_cselFingers_dom')
    szf.set('custom', True)
    szf.set('hmaxactive', True)
    szf.set('hmax', '0.25[um]*msc')
    szb = mesh.create('sizeBeam', 'Size')
    szb.selection().named('geom1_cselBeams_dom')
    szb.set('custom', True)
    szb.set('hmaxactive', True)
    szb.set('hmax', '0.5[um]*msc')
    tri = mesh.create('ftri1', 'FreeTri')
    mesh.run()
    return mesh


def build_mech_model(client, name='phase_shifter_p3_mech'):
    model = client.create(name)
    set_parameters(model)
    build_geometry(model, with_air=False)
    add_component_selections(model, with_air=False)
    add_solids_union(model)
    add_operators(model)
    add_solid(model)
    add_material(model, with_air=False)
    add_mesh(model, with_air=False)
    return model


def build_em_model(client, name='phase_shifter_p3_em'):
    """Solids + air, ES + electromechanical coupling + moving mesh."""
    model = client.create(name)
    set_parameters(model)
    build_geometry(model, with_air=True)
    add_component_selections(model, with_air=True)
    add_solids_union(model)
    jm = model.java
    comp = jm.component('comp1')

    # pure air = complement of solids (cselAir from the boolean Difference is
    # contaminated with mover domains -- Part-2 gotcha, never use directly)
    selair = comp.selection().create('selAirOnly', 'Complement')
    selair.set('entitydim', JInt(2))
    selair.set('input', ['selSolids'])

    add_operators(model)

    solid = comp.physics().create('solid', 'SolidMechanics', 'geom1')
    solid.selection().named('selSolids')
    solid.prop('Type2D').set('Type2D', 'PlaneStress')
    solid.prop('d').set('d', 't_si')
    fix = solid.create('fix1', 'Fixed', 1)
    fix.selection().named('selAnchors')

    es = comp.physics().create('es', 'Electrostatics', 'geom1')
    es.prop('d').set('d', 't_si')
    dter = es.create('dter1', 'DomainTerminal', 2)
    dter.selection().named('selSolids')
    dter.set('TerminalType', 'Voltage')
    dter.set('V0', '0[V]')
    pot = es.create('pot1', 'ElectricPotential', 1)
    pot.selection().named('selTerminal')
    pot.set('V0', 'Vd')

    # moving mesh first, then the coupling (selection auto-detection order)
    dd = comp.common().create('free1', 'DeformingDomain')
    dd.selection().named('selAirOnly')
    dd.set('smoothingType', 'hyperelastic')
    sym = comp.common().create('sym1', 'Symmetry')
    sym.selection().named('selAirSidesU')

    emf = comp.multiphysics().create('eme1', 'ElectromechanicalForces', 'geom1')

    add_material(model, with_air=True)
    add_mesh(model, with_air=True)
    return model


if __name__ == '__main__':
    mov, fix = comb_layout()
    print(f'movers: {len(mov)}, fixed: {len(fix)}')
    print(f'comb x-extent: {mov[0][0] - P["h"] - WC_END:.2f} .. '
          f'{fix[-1][0] + WC_END:.2f} um')
    print(f'mover y: [{Y_MOV_TIP}, {Y_SPINE_BOT}], fixed y: '
          f'[{Y_STATOR}, {Y_FIX_TIP}], overlap d0 = {Y_FIX_TIP - Y_MOV_TIP}')
