# Assignment Part 3/3 — geometry optimization + 3D out-of-plane analysis

Working directory for the ET4260 final assignment Part 3. Analytics in Python,
FEM script-driven against COMSOL 6.4 via MPh (same reproducible workflow as
Part 2).

## Layout

```
assignment_p3/
├── ANSWERS.md                 <- collected answers + numbers for the report
├── README.md                  <- this file
├── analysis/
│   ├── p3a_model.py           lumped-element model (V_DC_MAX, axial + side
│   │                          pull-in channels, Q, sensitivity), calibrated
│   │                          against Part-1 analytics and Part-2 FEM
│   ├── p3a_optimize.py        constrained grid-search optimizer (a.3)
│   └── p3a_final_design.py    final dimensions + summary CSVs + figures
├── comsol/
│   ├── models/                saved .mph models (open in the GUI)
│   │   ├── phase_shifter_p3_mech.mph    a.4 (i),(ii),(v)
│   │   ├── phase_shifter_p3_em.mph      a.4 (iii),(iv) + pull-in checks
│   │   └── phase_shifter_p3b_3d.mph     part b 3D model
│   └── scripts/
│       ├── model_lib_p3.py    shared geometry/physics builder (new dims)
│       ├── 01_mech_all.py     spring constant, 5 modes, step response
│       ├── 01b_build_em.py    EM model build
│       ├── 02_em_all.py       x(V) sweep, prestressed eig vs bias, max accel
│       ├── 03_p3b_3d.py       part b: t-sweep, 3D modes, k_ip/k_oop
│       └── 04_make_plots.py   regenerate matplotlib figures from CSVs
├── data/                      <- exported CSVs
└── figures/                   <- PNGs (COMSOL exports + matplotlib) + PDFs
```

## Running

```bash
# analytics (no COMSOL needed)
cd analysis
~/workspace/comsol-mcp/.venv/bin/python p3a_model.py        # calibration checks
~/workspace/comsol-mcp/.venv/bin/python p3a_optimize.py     # design search
~/workspace/comsol-mcp/.venv/bin/python p3a_final_design.py # freeze + figures

# FEM (COMSOL 6.4 license required; run sequentially - single license seat)
cd ../comsol/scripts
~/workspace/comsol-mcp/.venv/bin/python 01_mech_all.py
~/workspace/comsol-mcp/.venv/bin/python 01b_build_em.py
~/workspace/comsol-mcp/.venv/bin/python 02_em_all.py
~/workspace/comsol-mcp/.venv/bin/python 03_p3b_3d.py
~/workspace/comsol-mcp/.venv/bin/python 04_make_plots.py
```
