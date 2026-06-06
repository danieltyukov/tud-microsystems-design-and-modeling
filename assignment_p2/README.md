# Assignment Part 2/3 — FEM modelling (COMSOL) of the optical phase shifter

Working directory for the ET4260 final assignment Part 2. Everything is
script-driven against COMSOL 6.4 via MPh (Python–Java bridge), so every result
is reproducible end-to-end.

## Layout

```
assignment_p2/
├── ANSWERS.md            <- collected answers + numbers for the report
├── README.md             <- this file
├── comsol/
│   ├── GEOMETRY.md       <- device geometry derivation (P1 figure-faithful)
│   ├── models/           <- saved .mph models (open these in the GUI)
│   │   ├── phase_shifter_mech.mph        Q1-Q3 (solid mechanics only)
│   │   ├── phase_shifter_em.mph          Q4a, Q4b, Q5 (ES + EM coupling)
│   │   └── phase_shifter_em_step.mph     Q6 (voltage-step transient)
│   └── scripts/
│       ├── model_lib.py          shared geometry/material/mesh builder
│       ├── 01_build_mech.py      mechanics model build
│       ├── 02_q1_spring.py       Q1 spring constant + mesh convergence
│       ├── 03_q2_modes.py        Q2 eigenfrequencies + mode maps
│       ├── 04_q3_step.py         Q3 Rayleigh damping + force-step transient
│       ├── 05_build_em.py        EM model build (ES + coupling + moving mesh)
│       ├── 06_q4a_sweep.py       Q4a displacement vs voltage sweep
│       ├── 07_q4b_pullin.py      Q4b GE method (fails: 6.4 solver-gen bug, kept for documentation)
│       ├── exp_q4b_vsweep.py     Q4b ACTUAL: voltage continuation to pull-in
│       ├── 11_q4b_eigcheck.py    Q4b verification: prestressed eigenfreqs vs bias
│       ├── 08_q5_accel.py        Q5 max acceleration before contact
│       ├── 09_q6_vstep.py        Q6 voltage-step transients
│       └── 10_make_plots.py      regenerate all matplotlib figures from CSVs
├── data/                 <- exported CSVs (consumable by MATLAB for report plots)
└── figures/              <- quick-look PNGs (COMSOL exports + matplotlib)
```

## Running

```bash
cd /home/danieltyukov/workspace/comsol-mcp     # venv with mph installed
.venv/bin/python <script>.py
```

Scripts are ordered; 01/05 build models, the rest load + extend + save them.

## GUI screenshot checklist (for the report appendix)

Open the models in the COMSOL GUI and capture:

**phase_shifter_mech.mph**
- [ ] Geometry node with the parametrized rectangle list + final geometry plot
- [ ] Materials: Si assignment (E=170 GPa, ν=0.28, ρ=2330 kg/m³), plane stress, t=200 nm
- [ ] Fixed Constraint selection (4 beam anchor ends)
- [ ] Mesh (physics-controlled, refinement parameter msc)
- [ ] Q1 study: Stationary + point load; Q2 study: Eigenfrequency (6 modes)
- [ ] Q3: Linear Elastic Material > Damping (Rayleigh, β_dK = 3.866e-9 s) +
      Time Dependent study settings (range(0,1e-7,4e-4), rtol 1e-5)

**phase_shifter_em.mph**
- [ ] Electrostatics: Domain Terminal (shuttle, 0 V), Electric Potential
      (stator cavity walls, Vd), air domain selection `selAirOnly`
- [ ] Multiphysics: Electromechanical Forces coupling (all 27 solid domains)
- [ ] Moving Mesh: Deforming Domain (air only) + Symmetry on side openings
- [ ] Q4a study: Stationary + Auxiliary sweep Vd = range(0,10,150)
- [ ] Q5: Body Load (rho*a_acc) + stationary study
- [ ] Q6: Step function (1 µs rise), pot1 V0 = Vd*stepQ6(t[1/s]),
      Time Dependent study with Vd = {1, 50} V sweep

**phase_shifter_em.mph (Q4b)**
- [ ] stdVs study: Stationary + Auxiliary sweep Vd = range(0,5,450) with
      "Run continuation for: Last parameter" → diverges past pull-in
- [ ] x(V) plot from the partially completed sweep (last converged ≈ 222 V)

**phase_shifter_em_step.mph (Q6)**
- [ ] Step function stepQ6 (location 1 µs, smoothing 2 µs)
- [ ] pot1 V0 = Vd*stepQ6(t[1/s])
- [ ] Rayleigh damping node (same β_dK as Q3)
- [ ] Time Dependent study, Vd = {1, 50} V auxiliary sweep
```
