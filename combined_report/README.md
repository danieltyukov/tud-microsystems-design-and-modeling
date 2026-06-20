# ET4260 Final Report — combined (Parts 1+2+3)

`5714699_Tyukov_ET4260_Final.pdf` = Part 1 (pp. 1–11) + Part 2 (pp. 12–24) +
Part 3 (pp. 25–39), merged from the per-part reports:

```bash
pdfunite assignment_p1/report/5714699_Tyukov_ET4260_Part1.pdf \
         assignment_p2/report/5714699_Tyukov_ET4260_Part2.pdf \
         assignment_p3/report/5714699_Tyukov_ET4260_Part3.pdf \
         combined_report/5714699_Tyukov_ET4260_Final.pdf
```

Final submission deadline: 24 June 2026. Corrections to Parts 1/2 made after
their first submissions are marked in **blue**.

## Corrections in this final merge (addressing reviewer feedback)

**Part 1 — "some parts of proof mass is missing"** (eqs. 6, 9). The lumped model
now counts both 100×30 µm plates (modulation bar **and** comb-spine mass plate)
and the 20 moving fingers (N=40 = comb *gaps*). Mass 2.33 → **3.54 pg**;
cascades to f₀ 662 → **537 kHz**, Q 96.8 → **81.5**, settling 139 → **145 µs**,
shock limits, Bode peak. Electrostatics (pull-in 347 V, sensitivity) correctly
unchanged.

**Part 2 — boundaries / 4b / 6:**
- **Boundary-condition table** added (fixed constraint on the four beam-end
  edges; body loads over whole domains; ES terminals; moving-mesh rollers).
- **Q4b pull-in** reframed to the proper prestressed-eigenvalue method
  (electrostatic thickness = 1). The parametric sweep is demoted (mesh/step
  dependent). Rocking side-instability **224 V** explained as a symmetry-
  *breaking* bifurcation; **confirmed mesh-converged** by a re-run on a finer
  mesh (51 k elem, thickness = 1) giving 223–225 V (script `16_q4b_lean.py`,
  data `q4b_eig_vs_bias_fine.csv`, figure overlay).
- **Q6 1 V settling**: root-caused the first submission's unphysical 75 µs decay
  (generalised-α numerical dissipation at picometre amplitude). A **direct 1 V
  re-run** (reduced ρ∞→1, tightened rtol) reproduces the correct final value
  (1.55 pm) and near-ideal **1.92× overshoot**, confirming the 75 µs was an
  artefact. Settling quoted as **152 µs** from the well-resolved 50 V trace
  (= 1 V by linearity ≡ Q3 158 µs ≈ analytic 159 µs), since the pm-scale tail
  stays numerically sensitive. Scripts `19_q6_light.py` (2-core), data
  `q6_vstep_1V_direct.csv`, figures regenerated.

## Re-simulation environment

COMSOL 6.4 via TU Delft Class-Kit network license (`-ckl`, single concurrent
session — serialize jobs). **Requires eduVPN** to reach
`comsol-ckl.lic.tudelft.nl:27000`. Re-run scripts:
`assignment_p2/comsol/scripts/{16_q4b_lean,18_q6_final}.py`,
plotting `{plot_q4b_softening,plot_q6_fine}.py`. The direct 1 V transient
(`18_q6_final.py`) reproduces the 152 µs settling once a COMSOL seat is
available; the report already states the value via linearity.
