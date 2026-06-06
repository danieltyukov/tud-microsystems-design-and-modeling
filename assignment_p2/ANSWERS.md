# ET4260 Final Assignment — Part 2/3: FEM modelling of the optical phase shifter

**Daniel Tyukov (5714699) — complete working answers, data and figure index**

All COMSOL models built and solved programmatically (COMSOL 6.4, MPh/Java API) for full
reproducibility. Scripts: `comsol/scripts/` · Models: `comsol/models/` · Data: `data/` ·
Figures: `figures/`. Assignment source: `2026_final_assignament_part2_FINAL.pdf`.

Geometry change vs Part 1 (per assignment page 2): **the fixed comb fingers on each end
are $3 W_c = 3\,\mu\text{m}$ wide** (all other dimensions identical to Part 1's table).

---

## Model setup (common, assignment p.3 "Mechanical model and dynamics")

> *"Build a 2D FEM model in COMSOL, including the relevant elements for the system
> dynamics: the suspended mass, the beams and comb-drive structure. Exclude the optical
> waveguide. … For Silicon, assume an isotropic model."*

**Geometry (2D, top view)** — full derivation in `comsol/GEOMETRY.md`:
modulation bar $L \times W_m = 100 \times 30\,\mu\text{m}$, central column
$W_m \times L_m = 30 \times 40\,\mu\text{m}$, comb-spine mass plate
$100 \times 30\,\mu\text{m}$, 4 suspension beams $L_b \times W_b = 75 \times 5\,\mu\text{m}$
(outer ends fixed = anchors), 20 moving + 21 fixed comb fingers ($h = 1\,\mu\text{m}$ gaps,
$N = 40$, $d_0 = 19\,\mu\text{m}$ overlap, $L_c = 20\,\mu\text{m}$, $W_c = 1\,\mu\text{m}$,
end fixed fingers $3\,\mu\text{m}$), tip gaps $L_c - d_0 = 1\,\mu\text{m}$.
The optical waveguide is excluded (only the $g_0 = 200$ nm clearance above the bar
matters, in Q5).

**Physics**: Solid Mechanics, **plane stress** with thickness $t = 200$ nm (thin SOI
device layer: free top/bottom surfaces $\Rightarrow \sigma_{zz} \approx 0$), isotropic Si:
$E = 170$ GPa, $\rho = 2330\,\text{kg/m}^3$, $\nu = 0.28$ (assumed; the assignment gives
only $E$, $\rho$). Mesh: physics-controlled triangular, refinement factor swept for
convergence (4 516 / 12 709 / 44 563 elements). Mesh view: `figures/model_mech_mesh.png`.

**Electrostatics (Q4–Q6)**, following the two linked tutorials
(electrostatically-actuated cantilever; biased_resonator_2d_pull_in): air cavity
around the comb (Charge Conservation – Free Space), shuttle = **Domain Terminal**
at 0 V (also yields capacitance), stator floor + fixed-finger walls = Electric
Potential $V_d$, **Electromechanical Forces** multiphysics coupling on all solids,
**Moving Mesh** (Deforming Domain, hyperelastic smoothing) on the air with
Symmetry/Roller at the two side openings. Field plot: `figures/diag_potential_50V.png`.

**Mass bookkeeping** (recurring source of P1 ↔ FEM differences):
- P1 lumped mass: $m = \rho t (L W_m + L_m W_m + N L_c W_c) = \mathbf{2.33}$ **pg**
- figure-faithful shuttle: $\rho t (2 L W_m + L_m W_m + 20 L_c W_c) = \mathbf{3.54}$ **pg**
  (P1 counted one $100 \times 30$ plate instead of two, and 40 instead of 20 moving fingers)
- + beams $0.70$ pg (partially effective); FEM-integrated total: $4.24$ pg
  (`data/q1_spring_constant.csv`, column m_shuttle).

---

## Question 1 — Spring constant (stationary)

> *"Perform a stationary (DC) analysis to determine the spring constant. Compare with
> Part 1 and discuss any differences."*

**Theory (P1)**: four fixed–guided beams in parallel, $I = t W_b^3 / 12$:

$$k = 4 \cdot \frac{12 E I}{L_b^3} = \frac{4 E t W_b^3}{L_b^3} = 40.30\ \text{N/m}$$

**FEM**: $1\,\mu\text{N}$ test force on the mass along $y$; $k = F/u_y$. Mesh convergence
(`data/q1_spring_constant.csv`):

| mesh | elements | $k$ (N/m) |
|---|---|---|
| coarse | 4 516 | 36.99 |
| normal | 12 709 | 36.87 |
| fine | 44 563 | 36.79 |

**$k_{FEM} = 36.8$ N/m, −8.7% vs P1** (monotone from above $\Rightarrow$ converged;
displacement FEM is too stiff on coarse meshes).

**Why** (verification in `data/q1_verification.csv`): re-solving with plates/column
made 1000× stiffer gives $k = 39.83$ N/m (−1.1% from analytical):
- ≈7.6% — compliance of the "rigid" bodies the analytical model assumes (the guided
  beam end rotates/translates slightly with the plates),
- ≈1.1% — shear deformation + 2D continuum effects absent in Euler–Bernoulli theory.

## Question 2 — Eigenfrequency analysis, first 5 modes

> *"Determine the fundamental resonance frequency and obtain 2D maps of the modal
> displacement of the first 5 modes."*

**Theory (P1)**: $f_0 = \frac{1}{2\pi}\sqrt{k/m} = 661.9$ kHz with P1's $k = 40.3$ N/m,
$m = 2.33$ pg.

**FEM** (`data/q2_eigenfrequencies.csv`, two meshes, 0.1% converged):

| mode | $f$ (FEM) | character | 2D map |
|---|---|---|---|
| 1 | **496.26 kHz** | $y$-translation of shuttle (actuation mode), S-shaped beams | `figures/q2_mode1_496kHz.png` |
| 2 | 2 542.7 kHz | in-plane rocking + in-phase finger bending | `figures/q2_mode2_2543kHz.png` |
| 3 | 3 308.1 kHz | comb-finger cantilever band, collective in-phase | `figures/q2_mode3_3308kHz.png` |
| 4 | 3 311.2 kHz | finger band, alternating pattern | `figures/q2_mode4_3311kHz.png` |
| 5 | 3 311.2 kHz | finger band, higher spatial pattern | `figures/q2_mode5_3311kHz.png` |

**$f_1 = 496.3$ kHz vs P1's 661.9 kHz (−25%) — fully reconciled**:

$$\frac{f_{FEM}}{f_{P1}} = \sqrt{\frac{k_{FEM}}{k_{P1}} \cdot \frac{m_{P1}}{m_{eff}}} = \sqrt{0.913 \times 0.613} = 0.748 \;\Rightarrow\; 495\ \text{kHz} \checkmark$$

with $m_{eff} = k/(2\pi f_1)^2 = 3.79$ pg = figure-true shuttle $3.54$ pg + ≈36% of beam
mass (consistent with distributed-beam kinetic energy).
Hand check of the finger band: fixed-free cantilever $20 \times 1\,\mu\text{m}$
$\rightarrow f \approx 3.4$ MHz ✓. Near-degeneracy of modes 3–5: 20 nearly identical
fingers $\Rightarrow$ a band of closely spaced collective modes.

*2D-model limitation worth stating in the report*: out-of-plane modes are excluded by
construction; with $t = 200$ nm the out-of-plane beam stiffness is $(t/W_b)^2 = 1/625$
of in-plane $\Rightarrow$ the physical device has a much lower out-of-plane mode
(~20 kHz). The assignment explicitly requests the 2D model, so this is a documented
limitation.

## Question 3 — Rayleigh damping ($Q = 83$) + force-step response

> *"Include damping such that Q = 83 (hint: Rayleigh). Time-domain step response with a
> force at the suspended mass; compare with the analytically predicted settling time."*

**Damping design**: Rayleigh model $\zeta(f) = \frac{\alpha}{4\pi f} + \beta \pi f$.
Choosing pure stiffness-proportional damping ($\alpha = 0$) targeted at $f_1$:

$$\beta_{dK} = \frac{1}{2\pi f_1 Q} = \frac{1}{2\pi \cdot 496.26\ \text{kHz} \cdot 83} = 3.866 \times 10^{-9}\ \text{s} \;\Rightarrow\; \zeta(f_1) = \frac{1}{2 \cdot 83}$$

$\beta$-damping also damps the 3.3 MHz finger band $\propto f$ $\Rightarrow$ clean
single-mode response.

**Transient**: $1\,\mu\text{N}$ step (body load on the mass) from rest, 0–400 µs,
rtol $10^{-5}$. Data `data/q3_step_response.csv`, figure `figures/q3_step_response.png`.

| quantity | value |
|---|---|
| final displacement | 27.11 nm (static $F/k$ = 27.12 nm, 0.05% ✓) |
| overshoot | 1.96× (ideal underdamped $\approx 2 - \pi/(2Q)$) |
| **settling time (5%), FEM** | **158.4 µs** |
| analytic $t_s = \ln(20)/(\zeta\omega_1)$ with FEM $f_1$ | 159.5 µs (0.7% ✓) |
| P1 predicted ($Q = 96.8$, $\omega_0 = 4.16 \times 10^6$ rad/s) | 139.4 µs |

Same formula, different inputs: $t_s \propto Q/\omega_1$ — P1 used its computed
$Q = 96.8$ (the assignment now imposes $Q = 83$) and its higher $\omega_0$
(mass/stiffness bookkeeping, Q1–Q2).

## Question 4 — Displacement vs voltage; pull-in

> *"Include Electrostatics + moving mesh (tutorials linked). (a) Sweep 0–150 V, compare
> with P1 analytical on the same graph. (b) Increase voltage until pull-in; compare with
> P1's approximation (hint: biased_resonator_2d_pull_in tutorial)."*

**Theory (P1)**: ideal comb force $F = \frac{1}{2}V^2 \frac{dC}{dx}$ with
$\frac{dC}{dx} = \frac{N \epsilon_0 t}{h} = 7.083 \times 10^{-11}$ F/m:

$$x(V) = \frac{V^2}{2k}\frac{N \epsilon_0 t}{h}\,; \quad x(150\ \text{V}) = 19.8\ \text{nm}$$

Pull-in (P1, mover-tip gap model): $V_{PI} = 347$ V at $x_{PI} = 404$ nm.

### Q4a — $x(V)$, 0…150 V

Stationary auxiliary sweep $V_d$ = 0:10:150 V (continuation). Data
`data/q4a_displacement_vs_voltage.csv`, **combined graph (assignment requirement)**
`figures/q4a_displacement_vs_voltage.png`.

| $V_d$ | FEM $\lvert u_y \rvert$ | P1 analytical | ratio |
|---|---|---|---|
| 50 V | 4.04 nm | 2.20 nm | 1.84 |
| 100 V | 16.51 nm | 8.79 nm | 1.88 |
| 150 V | 38.60 nm | 19.78 nm | 1.95 |

FEM ≈ 1.85× analytical, decomposing **exactly** as:
- ×1.09 — softer FEM spring (Q1);
- ×1.68 — larger force gradient: the two **tip-gap families** (mover-tip→stator and
  fixed-tip→mass-plate, both $1\,\mu\text{m}$ and both closing with engagement) add
  $\approx \epsilon_0 W_c t / g^2$ per tip to the flank term. Force-implied $dC/dx$ at
  50 V: $1.19 \times 10^{-10}$ vs ideal-comb $7.08 \times 10^{-11}$ F/m.
- Terminal capacitance: $C(0\ \text{V}) = 1.511$ fF vs parallel-plate
  $N \epsilon_0 d_0 t / h = 1.346$ fF (+12% fringing); $C$ grows with $V$ (deeper
  engagement) ✓.
- The growing ratio (1.84→1.95) is the onset of the $1/g^2$ tip nonlinearity — the same
  physics that causes pull-in.

### Q4b — pull-in voltage

**Method used**: stationary continuation sweep $V_d$ = 0:5:450 V — each step initialized
from the previous solution tracks the physical stable branch; past the fold no stable
equilibrium exists and the solve diverges. Last converged voltage brackets $V_{PI}$.
Data `data/q4b_vsweep_curve.csv`, figure `figures/q4b_vsweep_branch.png`.

**Result: $V_{PI}$(FEM) = 222–227 V** (last converged 222.2 V at $u_y = -94$ nm), vs
**P1's 347 V**: ≈36% lower. Scaling P1's estimate with the FEM corrections:

$$V_{PI} \approx 347 \cdot \sqrt{\frac{k_{FEM}}{k_{P1}} \cdot \frac{(dC/dx)_{P1}}{(dC/dx)_{FEM}}} \approx 240\ \text{V}$$

— within ~8% of the FEM bracket, but see the mechanism below.

**Mechanism (verified by prestressed eigenfrequency analysis,
`data/q4b_eig_vs_bias.csv`, script `comsol/scripts/11_q4b_eigcheck.py`)**: solving the
biased eigenproblem at increasing $V_d$ shows which mode loses stiffness first
(instability $\Leftrightarrow f^2 \to 0$):

| $V_d$ | $f_1$ (shuttle, axial) | $f_2$ (rocking) | $f_3$ (finger band) |
|---|---|---|---|
| 0 V | 496.8 kHz | 2 544.5 kHz | 3 314 kHz |
| 100 V | 489.7 (−1.4%) | 2 379.0 (−6.5%) | 3 007 (−9.3%) |
| 180 V | 470.3 (−5.3%) | 1 685.1 (−34%) | 2 168 (−35%) |
| 210 V | 457.3 (−8.0%) | 989.7 (−61%) | 1 560 (−53%) |
| 220 V | 449.2 (−9.6%) | **508.7 (−80%)** | 1 265 (−62%) |

Softening plot: `figures/q4b_eig_softening.png`; rocking-mode shape at 220 V:
`figures/q4b_critical_mode.png`.

**$f_2^2 \to 0$ extrapolates to ≈ 224 V — matching the 222.2 V divergence exactly.**
The instability is the **in-plane rocking mode** (shuttle rotation with in-phase
lateral finger bending, the Q2 mode-2 shape): a small rotation moves every finger
toward one gap wall, and the lateral electrostatic negative stiffness
($\propto V^2 d_0 / h^3$ per finger, acting on long moment arms) overwhelms the
rotational suspension stiffness. The axial shuttle mode would survive to ≈445 V
($f_1^2$ extrapolation), and the finger band to ≈260 V — the rocking mode goes first.

This explains why the fold occurs at only 94 nm $\ll g/3 = 333$ nm, and why P1's axial
tip-gap model (347 V) — and any axial-only correction (≈240 V) — mispredicts both the
voltage and the character of pull-in: **the failure is rotational side pull-in, a
failure mode the 1-DOF analytical model cannot represent at all.** Design remedy:
stiffer rotational suspension (wider beam spacing) or shorter/stiffer fingers.

*Method note (for the report appendix)*: the assignment-linked displacement-controlled
continuation (Global Equation solving $V(x)$) was replicated exactly — point probe
(material frame), GE-disabled baseline at warm bias, segregated solver groups
[$V$+ODE]/[$u$]/[spatial] — but COMSOL 6.4's freshly generated solver sequences fail
with NaN on assembly for any global equation referencing displacement. Verified on
COMSOL's own biased_resonator_2d_pull_in model: its shipped solver sequence solves; a
fresh study with identical settings fails. The voltage-continuation method above is
unaffected. (Script kept for documentation: `comsol/scripts/07_q4b_pullin.py`.)

## Question 5 — Maximum acceleration before mass–waveguide contact

> *"Determine the maximum acceleration the system can withstand without contact between
> the suspended mass and the waveguide, along the direction of motion, for
> $V_d$ = 0, 1, 10, 50 V."*

**Theory**: contact when the modulation bar closes the $g_0 = 200$ nm gap:

$$a_{max} \approx \frac{k \,(g_0 + \lvert x_{dc}(V) \rvert)}{m_{eff}}$$

(the bias pulls the shuttle *away* from the waveguide).

**FEM**: body load $\rho a$ on all suspended Si ($+y$, toward waveguide); $a_{max}$ from
a linear solve pair, then **verification solve at $a_{max}$** (residual < 0.2%). Data
`data/q5_max_acceleration.csv`, figure `figures/q5_max_acceleration.png`.

| $V_d$ | $u_y$ at rest | $a_{max}$ | in g |
|---|---|---|---|
| 0 V | 0 | $1.9026 \times 10^6$ m/s² | 193 950 |
| 1 V | −0.002 nm | $1.9027 \times 10^6$ m/s² | 193 950 |
| 10 V | −0.160 nm | $1.9037 \times 10^6$ m/s² | 194 060 |
| 50 V | −4.03 nm | $1.9295 \times 10^6$ m/s² | 196 690 |

- $a_{max}$ **increases** with $V_d$ (+1.4% at 50 V): bias adds clearance and the
  tip-gap attraction relaxes as the shuttle rises. (P1 assumed the displacement reduces
  the gap — sign worth correcting in the revised P1.)
- vs naive $k g_0 / m_{shuttle} = 2.08 \times 10^6$ m/s²: FEM 8.7% lower because the
  distributed inertial load also bends the beams — $k g_0 / a_{max} \Rightarrow
  m_{eff} = 3.88$ pg ≈ shuttle + ½ beam mass.
- P1 reported $3.46 \times 10^6$ m/s² (353 000 g): difference dominated by P1's smaller
  mass (2.33 vs 3.54 pg) and stiffer $k$.

## Question 6 — Voltage-step response (1 V, 50 V) vs Q3

> *"Time-domain step response using a voltage applied to the comb-drive (1 V, 50 V).
> Compare with Question 3 and discuss differences in the dynamic response."*

Coupled solid + ES + moving-mesh transient, Rayleigh damping as Q3, smoothed voltage
step (1 µs rise $\ll$ 2 µs oscillation period; sharp ns steps stall the coupled DAE
integrator). Data `data/q6_vstep_1V.csv`, `data/q6_vstep_50V.csv`; figures
`figures/q6_vstep_responses.png`, normalized overlay vs Q3
`figures/q6_vs_q3_normalized.png`.

| quantity | 1 V step | 50 V step |
|---|---|---|
| final displacement | −1.61 pm | −3.97 nm |
| Q4a static value | −1.6 pm ✓ | −4.04 nm (−1.7%) |
| ringing frequency | 496 kHz = $f_1$ ✓ | 496 kHz = $f_1$ ✓ |
| settling (5%) | tolerance-limited (pm scale) | ≈ 150 µs |

**Comparison with Q3 (158.4 µs) — the requested discussion**:
1. **Dynamics are the same**: both inputs excite mode 1; ringing at $f_1$ and settling
   ≈ Q3. At 1 V and 50 V the device sits far below pull-in (222 V), so the system is
   effectively linear and the voltage step is simply a force step of amplitude
   $F = \frac{1}{2}V^2\,dC/dx$.
2. **Quadratic input scaling** (the fundamental difference from Q3): response amplitude
   $\propto V^2$, shown by $-1.6$ pm $\to -4.0$ nm for 50× voltage (2500× force ✓). A
   force input scales linearly; a voltage input quadratically — relevant for modulation
   linearity in the application.
3. **Lower overshoot (1.62× vs 1.96×)**: consequence of the 1 µs input rise filtering
   the step, not of different physics. Electrostatic spring softening at 50 V is
   negligible ($\Delta f_1/f_1 \approx -0.3$% from `data/q4b_eig_vs_bias.csv`),
   consistent with the unchanged ringing frequency.
4. *Numerical note*: the 50 V trace is quantitative (envelope decay within ≈7% of
   $\zeta\omega_1$). The 1 V response (−1.6 pm) sits below the integrator's effective
   absolute error floor: its waveform and final value are correct, but its apparent
   decay is dominated by numerical dissipation even at rtol $10^{-5}$ (verified by
   re-running — decay unchanged). Since the system is linear at these amplitudes
   (point 2), the physical 1 V transient is exactly the 50 V trace scaled by 1/2500,
   with the same ≈150–158 µs settling. Resolving picometer transients directly would
   require manual per-variable absolute tolerances — a known FEM practicality worth one
   sentence in the report.

---

## File index

| artifact | contents |
|---|---|
| `comsol/models/phase_shifter_mech.mph` | Q1–Q3 setup (solutions cleared to keep file small) |
| `comsol/models/phase_shifter_em.mph` | + ES/EM coupling/moving mesh; studies Q4a, Q5, **Q4b sweep (stdVs, with solutions)** |
| `comsol/models/phase_shifter_em_step.mph` | Q6 voltage-step transient setup |
| `data/q1_*.csv, q2_*.csv` | spring constant, convergence, verification, eigenfrequencies |
| `data/q3_step_response.csv` | force-step transient $u_y(t)$ |
| `data/q4a_*.csv, q4b_*.csv` | $x(V)$ sweep, pull-in branch, eigenfrequency-vs-bias |
| `data/q5_max_acceleration.csv` | $a_{max}(V_d)$ + verification residuals |
| `data/q6_vstep_{1V,50V}.csv` | voltage-step transients |
| `figures/*.png` | all report figures (regenerate: `comsol/scripts/10_make_plots.py`) |

GUI screenshot checklist for the Word report: see `README.md`.
