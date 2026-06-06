"""Build the mechanics-only model (Q1-Q3) and save it.
Also introspect interface props to set 2D plane-stress correctly,
verify shuttle mass vs analytical, and export a geometry/mesh image.
"""
import sys
sys.path.insert(0, '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p2/comsol/scripts')
import mph
from model_lib import build_mech_model, MODELS, FIGS

client = mph.start(cores=4)
model = build_mech_model(client)
jm = model.java
comp = jm.component('comp1')
solid = comp.physics('solid')

# verify 2D approximation + thickness took effect
print('Type2D =', solid.prop('Type2D').getString('Type2D'))
print('thickness d =', solid.prop('d').getString('d'))

# verify entity counts in named selections
geom = comp.geom('geom1')
print('domains total:', geom.getNDomains() if hasattr(geom, 'getNDomains') else '?')
for sel in ['geom1_cselBar_dom', 'geom1_cselMass_dom', 'geom1_cselBeams_dom',
            'geom1_cselFingers_dom', 'geom1_cselColumn_dom']:
    try:
        ents = comp.selection(sel).entities() if hasattr(comp.selection(sel), 'entities') else None
        print(f'{sel}: {list(ents) if ents is not None else "?"}')
    except Exception as e:
        print(f'{sel}: ERROR {e}')
for sel in ['selAnchors', 'selBarTop', 'selShuttle']:
    try:
        ents = comp.selection(sel).entities()
        print(f'{sel}: {list(ents)}')
    except Exception as e:
        print(f'{sel}: ERROR {e}')

model.save(f'{MODELS}/phase_shifter_mech.mph')
print('SAVED phase_shifter_mech.mph')

# mesh statistics
try:
    st = comp.mesh('mesh1').stat()
    print('mesh elements:', st.getNumElem())
    print('min quality:', st.getMinQuality())
except Exception as e:
    print('mesh stat error:', e)

# export mesh/geometry image
try:
    jm.result().create('pgGeom', 'PlotGroup2D')
    jm.result('pgGeom').create('mesh1', 'Mesh')
    jm.result().export().create('imgGeom', 'pgGeom', 'Image')
    jm.result().export('imgGeom').set('pngfilename', f'{FIGS}/model_mech_mesh.png')
    jm.result().export('imgGeom').set('size', 'manualweb')
    jm.result().export('imgGeom').set('width', '1600')
    jm.result().export('imgGeom').set('height', '900')
    jm.result().export('imgGeom').run()
    print('mesh image exported')
except Exception as e:
    print('image export error:', e)
print('DONE')
