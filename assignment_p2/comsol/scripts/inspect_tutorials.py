"""Deep-inspect the two COMSOL tutorial models that Part 2 references,
to extract the exact recipe (moving mesh, EM coupling, pull-in global equation).
Output: comsol/TUTORIAL_NOTES.txt
"""
import mph

APPS = '/usr/local/comsol64/multiphysics/applications/MEMS_Module/Actuators'
OUT = '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p2/comsol/TUTORIAL_NOTES.txt'

client = mph.start(cores=4)
lines = []


def w(s=''):
    print(s)
    lines.append(str(s))


def dump_node(node, depth=0, maxdepth=3, props=False):
    pad = '  ' * depth
    try:
        t = node.type()
    except Exception:
        t = '?'
    w(f"{pad}- {node.name()}  [{t}]")
    if props:
        try:
            for p in node.properties():
                try:
                    v = node.property(p)
                    vs = str(v)
                    if len(vs) > 120:
                        vs = vs[:120] + '...'
                    w(f"{pad}    {p} = {vs}")
                except Exception:
                    pass
        except Exception:
            pass
    if depth < maxdepth:
        for c in node.children():
            dump_node(c, depth + 1, maxdepth, props)


for fname, deep_props in [('electrostatically_actuated_cantilever.mph', True),
                          ('biased_resonator_2d_pull_in.mph', True)]:
    w('=' * 80)
    w(f'MODEL: {fname}')
    w('=' * 80)
    model = client.load(f'{APPS}/{fname}')

    for branch in ['parameters', 'functions', 'components', 'geometries',
                   'selections', 'physics', 'multiphysics', 'materials',
                   'meshes', 'studies', 'solutions']:
        w(f'\n### {branch} ###')
        try:
            node = model / branch
            for c in node.children():
                dump_node(c, 1, maxdepth=3,
                          props=deep_props and branch in
                          ('physics', 'multiphysics', 'studies'))
        except Exception as e:
            w(f'  (error: {e})')

    # definitions live under the component node
    w('\n### component definitions (moving mesh etc.) ###')
    try:
        comp = (model / 'components').children()[0]
        jcomp = comp.java
        # tags of definitions children
        for tag in jcomp.defs().tags() if hasattr(jcomp, 'defs') else []:
            w(f'  defs tag: {tag}')
    except Exception as e:
        w(f'  (error: {e})')
    # java-level: common.MovingMesh in COMSOL 6.x lives in model.common()
    try:
        jm = model.java
        w('  -- model.common() tags --')
        for tag in jm.common().tags():
            n = jm.common(tag)
            w(f'    {tag}: {n.label()}  type={n.getType() if hasattr(n,"getType") else "?"}')
            try:
                for p in n.properties():
                    vs = str(n.getString(str(p))) if True else ''
                    w(f'       {p} = {vs[:100]}')
            except Exception:
                pass
    except Exception as e:
        w(f'  (model.common error: {e})')

    client.remove(model)

with open(OUT, 'w') as f:
    f.write('\n'.join(lines))
w(f'\nSaved notes to {OUT}')
