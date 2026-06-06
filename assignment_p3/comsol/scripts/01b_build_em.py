"""Build the Part-3a EM model (solids + air + ES + EMF coupling + moving mesh).

Run: ~/workspace/comsol-mcp/.venv/bin/python 01b_build_em.py
"""
import sys
sys.path.insert(0, '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p3/comsol/scripts')
import mph
import model_lib_p3 as lib
from model_lib_p3 import MODELS

client = mph.start(cores=6)
model = lib.build_em_model(client)
nelem = model.java.component('comp1').mesh('mesh1').stat().getNumElem()
print(f'EM model built, mesh elements: {nelem}')
model.save(f'{MODELS}/phase_shifter_p3_em.mph')
print('saved', f'{MODELS}/phase_shifter_p3_em.mph')
