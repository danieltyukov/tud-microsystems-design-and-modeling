"""Clear stored solution data from saved .mph models (keeps setup, shrinks file).

Run after the study scripts: ~/workspace/comsol-mcp/.venv/bin/python 05_clean_models.py
"""
import sys
sys.path.insert(0, '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p3/comsol/scripts')
import os
import mph
from model_lib_p3 import MODELS

client = mph.start(cores=2)
for name in ['phase_shifter_p3_mech', 'phase_shifter_p3_em',
             'phase_shifter_p3b_3d']:
    path = f'{MODELS}/{name}.mph'
    if not os.path.exists(path):
        continue
    size0 = os.path.getsize(path) / 1e6
    model = client.load(path)
    jm = model.java
    for t in [str(x) for x in jm.sol().tags()]:
        try:
            jm.sol(t).clearSolutionData()
        except Exception as e:
            print(f'{name}/{t}: {str(e)[:60]}')
    model.save(path)
    print(f'{name}: {size0:.1f} MB -> {os.path.getsize(path)/1e6:.1f} MB')
    client.remove(model)
print('DONE clean')
