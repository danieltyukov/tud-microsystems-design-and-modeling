"""Shared builder for the ET4260 Part-2 comb-drive phase shifter FEM models.

Geometry per GEOMETRY.md (all coordinates in um, y=0 at bar top, motion = y).
Two model flavours:
  - mech: suspended solids only (Q1-Q3)
  - em:   + air domain, Electrostatics, Electromechanics coupling, moving mesh (Q4-Q6)
"""
import mph
from jpype import JInt

mph.option('classkit', True)   # TU Delft Class Kit License: start COMSOL server with -ckl

UM = 1e-6
BASE = '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p2'
MODELS = f'{BASE}/comsol/models'
DATA = f'{BASE}/data'
FIGS = f'{BASE}/figures'

# ----- parameters (um where geometric) -----
P = dict(
    L=100, Wm=30, Lm=40, Lb=75, Wb=5,
    Lc=20, Wc=1, h=1, d0=19, g0=0.2,
    t_si=0.2, N=40,
)
MAT = dict(E='170[GPa]', nu='0.28', rho='2330[kg/m^3]')

N_MOV = 20          # moving fingers (N=40 gaps -> 20 movers, 21 fixed)
N_FIX_IN = 19       # inner fixed fingers (1 um)
WC_END = 3          # end fixed finger width = 3*Wc  (Part-2 change)


def comb_layout():
    """Return (moving, fixed) lists of (x_left, width) in um."""
    x = -42.5
    fixed = [(x, WC_END)]          # left end finger 3 um
    x += WC_END + 1                # + gap
    moving = []
    for i in range(N_MOV):
        moving.append((x, 1.0))
        x += 1 + 1                 # finger + gap
        if i < N_MOV - 1:
            fixed.append((x, 1.0))
            x += 1 + 1
    fixed.append((x, WC_END))      # right end finger
    assert abs((x + WC_END) - 42.5) < 1e-9, f'comb layout end mismatch: {x + WC_END}'
    return moving, fixed


def rect(geom, tag, x, y, w, hgt, contribute=None):
    r = geom.create(tag, 'Rectangle')
    r.set('pos', [f'{x}[um]', f'{y}[um]'])
    r.set('size', [f'{w}[um]', f'{hgt}[um]'])
    if contribute:
        r.set('contributeto', contribute)
    return r


def build_geometry(model, with_air):
    """Create geometry sequence; returns java geom object."""
    jm = model.java
    comp = jm.component().create('comp1', True)
    geom = comp.geom().create('geom1', 2)
    geom.lengthUnit('m')

    # cumulative selections
    for cs in ['cselShuttle', 'cselMass', 'cselBar', 'cselBeams',
               'cselFingers', 'cselColumn'] + (['cselAir'] if with_air else []):
        geom.selection().create(cs, 'CumulativeSelection')

    rect(geom, 'bar', -50, -30, 100, 30, 'cselBar')
    rect(geom, 'col', -15, -70, 30, 40, 'cselColumn')
    rect(geom, 'mas', -50, -100, 100, 30, 'cselMass')
    # 1 um release slots between beams and bar/mass plates (figure shows
    # the beams attached to the COLUMN only, separated from the plates)
    rect(geom, 'bUL', -90, -36, 75, 5, 'cselBeams')
    rect(geom, 'bUR', 15, -36, 75, 5, 'cselBeams')
    rect(geom, 'bLL', -90, -69, 75, 5, 'cselBeams')
    rect(geom, 'bLR', 15, -69, 75, 5, 'cselBeams')

    moving, fixed = comb_layout()
    mov_tags = []
    for i, (xl, w) in enumerate(moving):
        tag = f'mf{i}'
        rect(geom, tag, xl, -120, w, 20, 'cselFingers')
        mov_tags.append(tag)

    if with_air:
        # air bounding rect minus fixed fingers (cavities) minus moving fingers
        rect(geom, 'airR', -46, -121, 92, 21)
        fix_tags = []
        for i, (xl, w) in enumerate(fixed):
            tag = f'ff{i}'
            rect(geom, tag, xl, -121, w, 20)  # fixed fingers y in [-121,-101]
            fix_tags.append(tag)
        # copies of moving fingers to subtract from air
        dif = geom.create('difAir', 'Difference')
        dif.selection('input').set('airR')
        dif.selection('input2').set(*(fix_tags))
        dif.set('contributeto', 'cselAir')
        dif.set('keep', False)
        # moving fingers overlap the diff result; partition happens at form union
    geom.run('fin')
    return geom


def add_shuttle_union_selection(model):
    """Selection helpers at component level (boxes use meters)."""
    jm = model.java
    comp = jm.component('comp1')

    def box(tag, x1, x2, y1, y2, dim):
        s = comp.selection().create(tag, 'Box')
        s.set('entitydim', JInt(dim))
        s.set('xmin', x1 * UM); s.set('xmax', x2 * UM)
        s.set('ymin', y1 * UM); s.set('ymax', y2 * UM)
        s.set('condition', 'inside')
        return s

    # anchors: beam outer end edges at x=+-90
    box('selAnchorL', -90.01, -89.99, -71, -29, 1)
    box('selAnchorR', 89.99, 90.01, -71, -29, 1)
    u = comp.selection().create('selAnchors', 'Union')
    u.set('entitydim', JInt(1))
    u.set('input', ['selAnchorL', 'selAnchorR'])
    # bar top edge (waveguide-facing)
    box('selBarTop', -50.01, 50.01, -0.01, 0.01, 1)
    return comp


def add_terminal_selection(model):
    """Stator (electrode) boundaries: cavity walls of 21 fixed fingers + stator line."""
    jm = model.java
    comp = jm.component('comp1')
    _, fixed = comb_layout()
    tags = []

    def box(tag, x1, x2, y1, y2):
        s = comp.selection().create(tag, 'Box')
        s.set('entitydim', JInt(1))
        s.set('xmin', x1 * UM); s.set('xmax', x2 * UM)
        s.set('ymin', y1 * UM); s.set('ymax', y2 * UM)
        s.set('condition', 'inside')
        tags.append(tag)

    box('selStatLine', -46.01, 46.01, -121.01, -120.99)
    for i, (xl, w) in enumerate(fixed):
        box(f'selFF{i}', xl - 0.2, xl + w + 0.2, -121.2, -100.8)
    u = comp.selection().create('selTerminal', 'Union')
    u.set('entitydim', JInt(1))
    u.set('input', tags)
    # air side openings (mesh roller)
    s = comp.selection().create('selAirSides', 'Box')
    s.set('entitydim', JInt(1))
    s.set('xmin', -46.01 * UM); s.set('xmax', -45.99 * UM)
    s.set('ymin', -121.01 * UM); s.set('ymax', -99.9 * UM)
    s.set('condition', 'inside')
    s2 = comp.selection().create('selAirSides2', 'Box')
    s2.set('entitydim', JInt(1))
    s2.set('xmin', 45.99 * UM); s2.set('xmax', 46.01 * UM)
    s2.set('ymin', -121.01 * UM); s2.set('ymax', -99.9 * UM)
    s2.set('condition', 'inside')
    u2 = comp.selection().create('selAirSidesU', 'Union')
    u2.set('entitydim', JInt(1))
    u2.set('input', ['selAirSides', 'selAirSides2'])


def set_parameters(model, extra=None):
    pars = {
        'L': '100[um]', 'Wm': '30[um]', 'Lm': '40[um]', 'Lb': '75[um]',
        'Wb': '5[um]', 'Lc': '20[um]', 'Wc': '1[um]', 'h_gap': '1[um]',
        'd0': '19[um]', 'g0': '200[nm]', 't_si': '200[nm]',
        'F_test': '1[uN]',           # Q1 test force
        'Q_target': '83',            # Q3
        'f1_fem': '500[kHz]',        # placeholder, set after Q2
        'beta_dk': '1/(2*pi*f1_fem*Q_target)',
        'alpha_dm': '0[1/s]',
        'Vd': '0[V]',                # Q4 drive voltage
        'a_acc': '0[m/s^2]',         # Q5 acceleration
        'msc': '1',                  # mesh scale factor (convergence study)
    }
    if extra:
        pars.update(extra)
    for k, v in pars.items():
        model.parameter(k, v)


def add_operators(model):
    """Integration/average operators for result extraction + variables."""
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
    itg.selection().named('geom1_cselShuttle_dom')
    itg.set('opname', 'intShuttle')
    var = comp.variable().create('var1')
    var.set('uy_mass', 'aveMass(v)')
    var.set('uy_bar', 'aveBarTop(v)')
    var.set('m_shuttle', 'intShuttle(solid.rho*t_si)')


def collect_shuttle_into_cumulative(model):
    """cselShuttle = union of all solid parts (domain level, via component Union)."""
    jm = model.java
    comp = jm.component('comp1')
    u = comp.selection().create('selShuttle', 'Union')
    u.set('entitydim', JInt(2))
    u.set('input', ['geom1_cselBar_dom_uni', ])  # replaced below


def add_solid(model, with_air):
    jm = model.java
    comp = jm.component('comp1')
    solid = comp.physics().create('solid', 'SolidMechanics', 'geom1')
    # restrict to suspended solids (everything except air)
    if with_air:
        u = comp.selection().create('selSolids', 'Union')
        u.set('entitydim', JInt(2))
        u.set('input', ['geom1_cselBar_dom', 'geom1_cselColumn_dom',
                        'geom1_cselMass_dom', 'geom1_cselBeams_dom',
                        'geom1_cselFingers_dom'])
        solid.selection().named('selSolids')
    # thickness + 2D approximation (thin SOI device layer -> plane stress)
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
    if with_air:
        # assign Si to solids only; air gets vacuum permittivity material
        comp.material('matSi').selection().named('selSolids')
        air = comp.material().create('matAir', 'Common')
        air.label('Air')
        air.propertyGroup('def').set('relpermittivity',
                                     ['1', '0', '0', '0', '1', '0', '0', '0', '1'])
        air.propertyGroup('def').set('density', '1.2[kg/m^3]')
        air.selection().named('geom1_cselAir_dom')
        # Si needs relative permittivity for ES domain terminal coverage
        mat.propertyGroup('def').set('relpermittivity',
                                     ['11.7', '0', '0', '0', '11.7', '0', '0', '0', '11.7'])
    return mat


def add_mesh(model, with_air):
    jm = model.java
    comp = jm.component('comp1')
    mesh = comp.mesh().create('mesh1')
    # global default size
    sz = mesh.feature('size')
    sz.set('custom', True)
    sz.set('hmax', '3[um]*msc')
    sz.set('hmin', '0.05[um]')
    # fine size in fingers (+ air if present)
    szf = mesh.create('sizeFine', 'Size')
    if with_air:
        szf.selection().geom('geom1', 2)
        # named union: fingers + air
        u = comp.selection().create('selFine', 'Union')
        u.set('entitydim', JInt(2))
        u.set('input', ['geom1_cselFingers_dom', 'geom1_cselAir_dom'])
        szf.selection().named('selFine')
    else:
        szf.selection().named('geom1_cselFingers_dom')
    szf.set('custom', True)
    szf.set('hmaxactive', True)
    szf.set('hmax', '0.45[um]*msc')
    # beams: at least 3-4 elements across width
    szb = mesh.create('sizeBeam', 'Size')
    szb.selection().named('geom1_cselBeams_dom')
    szb.set('custom', True)
    szb.set('hmaxactive', True)
    szb.set('hmax', '1.2[um]*msc')
    tri = mesh.create('ftri1', 'FreeTri')
    mesh.run()
    return mesh


def build_mech_model(client, name='phase_shifter_mech'):
    model = client.create(name)
    set_parameters(model)
    build_geometry(model, with_air=False)
    add_shuttle_union_selection(model)
    # cumulative selection union for shuttle = all domains (mech model has no air)
    jm = model.java
    comp = jm.component('comp1')
    u = comp.selection().create('selShuttle', 'Union')
    u.set('entitydim', JInt(2))
    u.set('input', ['geom1_cselBar_dom', 'geom1_cselColumn_dom',
                    'geom1_cselMass_dom', 'geom1_cselBeams_dom',
                    'geom1_cselFingers_dom'])
    add_operators_mech(model)
    add_solid(model, with_air=False)
    add_material(model, with_air=False)
    add_mesh(model, with_air=False)
    return model


def build_em_model(client, name='phase_shifter_em'):
    """Full electromechanical model: solids + air, ES + EMF coupling + moving mesh.

    Recipe mirrors MEMS-module tutorials (electrostatically_actuated_cantilever,
    biased_resonator_2d_pull_in): ES on air + shuttle solids, shuttle = Domain
    Terminal (V=0, gives charge readout), stator cavity boundaries = Electric
    Potential Vd, DeformingDomain on air with hyperelastic smoothing.
    """
    model = client.create(name)
    set_parameters(model)
    build_geometry(model, with_air=True)
    add_shuttle_union_selection(model)
    add_terminal_selection(model)
    jm = model.java
    comp = jm.component('comp1')

    # selSolids = all suspended silicon domains
    u = comp.selection().create('selSolids', 'Union')
    u.set('entitydim', JInt(2))
    u.set('input', ['geom1_cselBar_dom', 'geom1_cselColumn_dom',
                    'geom1_cselMass_dom', 'geom1_cselBeams_dom',
                    'geom1_cselFingers_dom'])
    # pure air = complement of solids. NOTE: cselAir from the geometry
    # Difference output is CONTAMINATED with the moving-finger domains
    # (cumulative selections inherit all final domains carved from the
    # boolean output object) - never use it for physics selections.
    selair = comp.selection().create('selAirOnly', 'Complement')
    selair.set('entitydim', JInt(2))
    selair.set('input', ['selSolids'])

    # operators / variables
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

    # solid mechanics on solids only
    solid = comp.physics().create('solid', 'SolidMechanics', 'geom1')
    solid.selection().named('selSolids')
    solid.prop('Type2D').set('Type2D', 'PlaneStress')
    solid.prop('d').set('d', 't_si')
    fix = solid.create('fix1', 'Fixed', 1)
    fix.selection().named('selAnchors')

    # electrostatics on air + solids (solids become domain terminal)
    es = comp.physics().create('es', 'Electrostatics', 'geom1')
    es.prop('d').set('d', 't_si')
    dter = es.create('dter1', 'DomainTerminal', 2)
    dter.selection().named('selSolids')
    dter.set('TerminalType', 'Voltage')
    dter.set('V0', '0[V]')
    pot = es.create('pot1', 'ElectricPotential', 1)
    pot.selection().named('selTerminal')
    pot.set('V0', 'Vd')

    # moving mesh FIRST (its selection influences the coupling's auto-selection):
    # deforming air domain + sliding side openings
    dd = comp.common().create('free1', 'DeformingDomain')
    dd.selection().named('selAirOnly')
    dd.set('smoothingType', 'hyperelastic')
    sym = comp.common().create('sym1', 'Symmetry')
    sym.selection().named('selAirSidesU')

    # electromechanical forces coupling (auto-selection = all solid domains,
    # i.e. solid n es minus the deforming-mesh region)
    emf = comp.multiphysics().create('eme1', 'ElectromechanicalForces', 'geom1')

    # materials
    mat = comp.material().create('matSi', 'Common')
    mat.label('Si isotropic')
    mat.propertyGroup('def').set('density', MAT['rho'])
    mat.propertyGroup('def').set('relpermittivity',
                                 ['11.7', '0', '0', '0', '11.7', '0', '0', '0', '11.7'])
    enu = mat.propertyGroup().create('Enu', "Young's modulus and Poisson's ratio")
    enu.set('E', MAT['E'])
    enu.set('nu', MAT['nu'])
    mat.selection().named('selSolids')
    air = comp.material().create('matAir', 'Common')
    air.label('Air')
    air.propertyGroup('def').set('relpermittivity',
                                 ['1', '0', '0', '0', '1', '0', '0', '0', '1'])
    air.propertyGroup('def').set('density', '1.2[kg/m^3]')
    air.selection().named('selAirOnly')

    add_mesh(model, with_air=True)
    return model


def add_operators_mech(model):
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
    itg.selection().named('selShuttle')
    itg.set('opname', 'intShuttle')
    var = comp.variable().create('var1')
    var.set('uy_mass', 'aveMass(v)')
    var.set('uy_bar', 'aveBarTop(v)')
    var.set('m_shuttle', 'intShuttle(solid.rho*t_si)')
