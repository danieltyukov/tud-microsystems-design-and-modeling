"""Q1: stationary spring constant + mesh convergence study.

Applies F_test (total) on the mass plate along -y (comb-pull direction),
solves stationary at 3 mesh refinements, k = F/|uy|.
"""
import sys
sys.path.insert(0, '/home/danieltyukov/workspace/tud/tud-microsystems-design-and-modeling/assignment_p2/comsol/scripts')
import mph
import csv
from model_lib import MODELS, DATA

client = mph.start(cores=4)
model = client.load(f'{MODELS}/phase_shifter_mech.mph')
jm = model.java
comp = jm.component('comp1')
solid = comp.physics('solid')

# ---- body load with specified TOTAL force on the mass plate ----
bl = solid.create('blQ1', 'BodyLoad', 2)
bl.selection().named('geom1_cselMass_dom')
# total F_test as force-per-volume over the mass plate (exact: V = L*Wm*t_si)
bl.set('forceReferenceVolume', ['0', '-F_test/(L*Wm*t_si)', '0'])

# ---- stationary study ----
std = jm.study().create('stdQ1')
stat = std.create('stat', 'Stationary')

rows = []
for msc in [2.0, 1.0, 0.5]:
    model.parameter('msc', str(msc))
    comp.mesh('mesh1').run()
    nelem = comp.mesh('mesh1').stat().getNumElem()
    jm.study('stdQ1').run()
    uy_mass = model.evaluate('uy_mass', 'm')
    uy_bar = model.evaluate('uy_bar', 'm')
    m_sh = model.evaluate('m_shuttle', 'kg')
    F = 1e-6
    k_mass = F / abs(float(uy_mass))
    k_bar = F / abs(float(uy_bar))
    rows.append(dict(msc=msc, nelem=int(nelem), uy_mass=float(uy_mass),
                     uy_bar=float(uy_bar), k_mass=k_mass, k_bar=k_bar,
                     m_shuttle=float(m_sh)))
    print(f'msc={msc}: nelem={nelem}, uy_mass={float(uy_mass):.6e} m, '
          f'k={k_mass:.4f} N/m, k_bar={k_bar:.4f}, m={float(m_sh):.4e} kg')

with open(f'{DATA}/q1_spring_constant.csv', 'w', newline='') as f:
    wcsv = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    wcsv.writeheader()
    wcsv.writerows(rows)

# analytical references
k_ana = 4 * 12 * 170e9 * (200e-9 * (5e-6)**3 / 12) / (75e-6)**3
print(f'\nAnalytical k (P1, 4 fixed-guided beams): {k_ana:.4f} N/m')
print(f'FEM k (finest): {rows[-1]["k_mass"]:.4f} N/m  '
      f'({(rows[-1]["k_mass"]/k_ana-1)*100:+.2f} %)')
print(f'Shuttle mass FEM: {rows[-1]["m_shuttle"]:.4e} kg '
      f'(figure-true 3.542e-12, P1 2.33e-12)')

model.save(f'{MODELS}/phase_shifter_mech.mph')
print('DONE Q1')
