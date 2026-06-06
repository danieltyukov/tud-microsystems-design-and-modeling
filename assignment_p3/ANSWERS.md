# ET4260 Final Assignment — Part 3/3: Geometry optimization + 3D out-of-plane analysis

**Daniel Tyukov (5714699) — complete working answers, data and figure index**

Analytical model + optimizer: `analysis/` (Python). All COMSOL models built and solved
programmatically (COMSOL 6.4, MPh/Java API): `comsol/scripts/`, models in
`comsol/models/`. Data: `data/` · Figures: `figures/`. Assignment source:
`2026_final_assignament_part3_FINAL.pdf`.

Technology update vs Parts 1/2: **minimum feature size 200 nm** (all widths and gaps),
features $\le$ 1 mm, $2 < N < 200$, suspended thickness fixed at $t=200$ nm in part a.

---

## Part a) 2D optimization of the phase shifter

### a.1 — Low-frequency modulation sensitivity vs $V_{DC}$; analytical $V_{DC\_MAX}$

From Part 1 (Q8), the small-signal transfer function is

$$H(s)=\frac{\Delta\phi(s)}{V_{AC}(s)}=\frac{S_\phi\,\eta}{ms^2+bs+k},\qquad
|H|_{f\to0}(V_{DC})=\frac{S_\phi(g_{op})\,\eta(V_{DC})}{k}$$

with the electromechanical coupling and the optical phase sensitivity

$$\eta=V_{DC}\frac{dC}{dx}=V_{DC}\frac{N\varepsilon_0 t}{h},\qquad
S_\phi=\frac{2\pi L}{\lambda}\,\Delta\eta_{eff}\,b\,e^{-b\,g_{op}},\qquad
\Delta\eta_{eff}=\eta_{eff,0}-\eta_{eff,\infty}=0.2 .$$

The comb drive pulls the shuttle **into** the comb, i.e. *away* from the waveguide
(sign verified by the Part-2 FEM), so the operating gap **grows** with bias:

$$g_{op}=g_0+x(V_{DC}),\qquad x(V_{DC})=\frac{V_{DC}^2}{2k}\frac{N\varepsilon_0 t}{h}
\equiv\alpha V_{DC}^2 .$$

The sensitivity is therefore a competition between the coupling ($\propto V_{DC}$)
and the exponentially decaying optical overlap:

$$\boxed{\;|H|_{f\to0}(V_{DC})=\frac{2\pi L}{\lambda}\Delta\eta_{eff}\,b\,
e^{-b(g_0+\alpha V_{DC}^2)}\,\frac{N\varepsilon_0 t}{kh}\,V_{DC}\;}$$

Setting $d|H|/dV_{DC}=0$ gives $1-2b\alpha V_{DC}^2=0$, i.e.

$$\boxed{\;V_{DC\_MAX}=\frac{1}{\sqrt{2b\alpha}}=\sqrt{\frac{k\,h}{b\,N\,\varepsilon_0\,t}}\;}$$

Two remarkable corollaries:

1. **The optimum displacement is a purely optical quantity**:
   $x(V_{DC\_MAX})=\alpha V_{DC\_MAX}^2=\dfrac{1}{2b}=1.67\,\mu\text{m}$,
   independent of every geometric parameter.
2. The achievable peak sensitivity has the closed form
   $$|H|_{max}=\frac{2\pi L}{\lambda}\,\Delta\eta_{eff}\,
   \frac{e^{-b g_0-1/2}}{V_{DC\_MAX}},$$
   so maximizing sensitivity $=$ maximizing $L$ and *minimizing* $V_{DC\_MAX}$
   (i.e. minimizing $kh/N$) — within the constraints below.

For the Part-1 dimensions: $V_{DC\_MAX}=\sqrt{40.3\cdot10^{-6}/(0.3\cdot10^{6}\cdot40
\cdot\varepsilon_0\cdot200\,\text{nm})}=1377$ V, far above the derived pull-in
(347 V analytical, 222 V in the Part-2 FEM) — exactly the mismatch the assignment
points out. Figure: `figures/p3a_sensitivity_vs_V.png` (left panel).

### a.2 — Conditions for $V_{DC\_MAX}<V_{Pull-In}$

With the Part-1 finger-tip pull-in model ($u$ = tip gap at pull-in,
$g_{t0}=L_c-d_0$ the tip clearance at rest):

$$V_{PI}=\sqrt{\frac{k\,u^3}{N\varepsilon_0 W_c t}},\qquad
g_{t0}=\frac{3u}{2}+\frac{u^3}{2hW_c}.$$

Taking the ratio with $V_{DC\_MAX}$, **the stiffness $k$, comb size $N$, thickness
$t$ and $\varepsilon_0$ all cancel**:

$$\left(\frac{V_{DC\_MAX}}{V_{PI}}\right)^2=\frac{h\,W_c}{b\,u^3}
\quad\Longrightarrow\quad
\boxed{\;V_{DC\_MAX}<V_{PI}\iff b\,u^3>h\,W_c\;}$$

The equivalent (and more intuitive) displacement form: since $x(V)$ is monotone,
$V_{DC\_MAX}<V_{PI}$ exactly when the stable comb travel exceeds the optical optimum:

$$\boxed{\;x_{PI}=g_{t0}-u>\frac{1}{2b}=1.67\,\mu\text{m}\;}$$

In the deep-clearance limit ($u^3/2hW_c\gg 3u/2$, $u^3\to 2hW_c\,g_{t0}$) the condition
collapses to simply $g_{t0}>1/(2b)$.

**Critical terms and the resulting design constraints:**

| term | role | constraint |
|---|---|---|
| $g_{t0}=L_c-d_0$ | sets the available travel | must be $\gg 1/(2b)$: Part 1 had $1\,\mu$m $\to x_{PI}=0.40\,\mu$m (fails); need $g_{t0}\gtrsim 4\,\mu$m |
| $h\,W_c$ | strength of the destabilizing tip force | smaller $hW_c$ moves pull-in deeper into the travel |
| $b$ | optics: sets both $x_{opt}=1/(2b)$ and the margin | material/optics constant, not adjustable |
| $k,N,t$ | **cancel** | cannot fix the problem by stiffening or resizing the comb |

The cancellation is the key insight: Part 1's design cannot be rescued by a stiffer
spring — only by **increasing the tip clearance** $L_c-d_0$ (and/or reducing $hW_c$).

**Side pull-in (lesson from the Part-2 FEM).** The axial tip model is necessary but
not sufficient: Part 2 showed the real instability of the old design is the lateral
*rocking* mode at 222 V (vs 347 V axial). For the redesign, four additional
side-stability channels are modelled analytically and imposed as constraints
(all per finger-gap negative stiffness $k_{neg}=2\varepsilon_0 t\,d_{op}V^2/h^3$,
engagement $d_{op}=d_0+1/(2b)$):

1. **individual finger bending** — Rayleigh quotient of a cantilever loaded by
   $k_{neg}$ over the engaged length;
2. **shuttle rocking** — beam axial stiffness on the lever arms of the two beam
   levels (column bending in series), vs $\sum k_{neg,i}x_i^2$ over the comb;
3. **spine bending** — the comb spine is a beam on a *negative elastic foundation*
   $\kappa=k_{neg}/p$ ($p$ = comb pitch); a cantilever half of length $\ell$ buckles
   at $\kappa_{crit}=(1.875/\ell)^4EI_{spine}$ — this is what forbids very slim
   plates;
4. **lateral shuttle translation** — beam axial stiffness vs $n_{mov}k_{neg}$.

Combined (Dunkerley) and calibrated once against the Part-2 FEM (model 247 V vs
FEM 222 V → factor 0.90), the model reproduces the old design's side pull-in to 1%.

### a.3 — Optimized dimensions

Maximize $|H|_{max}\propto L\,e^{-bg_0}/V_{DC\_MAX}$ subject to: the four assignment
constraints ($f_r>500$ kHz, $V_{PI}>200$ V, $V_{DC\_MAX}<V_{PI}$,
$a_{0,max}(\omega_r)>2\times10^4$ m/s²), the side-stability channels above with a
15% margin, internal-mode separation (bar/spine tips, fingers, beams $\ge2.5f_r$),
the comb fitting under the spine ($\approx N(W_c+h)\le L$), and the technology
window (features 200 nm–1 mm, $2<N<200$, $t=200$ nm). Full grid search
(`analysis/p3a_optimize.py`, $\sim10^6$ candidates) + rounding to clean values
(`analysis/p3a_final_design.py`).

**Final values (assignment table, page 2):**

| Feature | Symbol | Part 1/2 | **Part 3a final** |
|---|---|---|---|
| Length of modulation section | $L$ | 100 µm | **100 µm** |
| Modulation gap at rest | $g_0$ | 200 nm | **200 nm** |
| Length of beams | $L_b$ | 75 µm | **40 µm** |
| Width of beams | $W_b$ | 5 µm | **1.85 µm** |
| Length of mass central section | $L_m$ | 40 µm | **30 µm** |
| Width of mass sections | $W_m$ | 30 µm | **10 µm** |
| Comb-drive gap | $h$ | 1 µm | **0.9 µm** |
| Comb-drive overlap at rest | $d_0$ | 19 µm | **1.0 µm** |
| Comb-drive finger length | $L_c$ | 20 µm | **5.5 µm** |
| Comb-drive finger width | $W_c$ | 1 µm | **0.7 µm** |
| Bottom end fingers | $3W_c$ | 3 µm | **2.1 µm** |
| Thickness of suspended parts | $t$ | 200 nm | **200 nm** (fixed) |
| Number of comb gaps | $N$ | 40 | **58** |

Design logic in one paragraph: $g_0$ stays at the 200 nm minimum (sensitivity
$\propto e^{-bg_0}$); $L$ stays at 100 µm because the plate mass needed for spine
side-stability grows $\propto L^{7/3}$, making $|H|_{max}$ nearly flat in $L$ —
the extra optical length is eaten by the higher $V_{DC\_MAX}$; the comb is completely
re-proportioned: short stubby fingers ($5.5\times0.7$ µm, 32 MHz internal mode) at
small overlap $d_0=1$ µm with a large tip clearance $g_{t0}=4.5$ µm give
$x_{PI}=3.07$ µm $> 2\cdot$ the required travel, while $N=58$ gaps at $h=0.9$ µm
fill the spine; the suspension is softened ($k=13.5$ N/m vs 40.3) exactly to the
$f_r$ floor to minimize $V_{DC\_MAX}=\sqrt{kh/(bN\varepsilon_0t)}$; plates slim to
$W_m=10$ µm — the side-stability limit.

**Predicted performance and constraint check (analytical):**

| quantity | Part 1/2 | Part 3a | constraint |
|---|---|---|---|
| $k$ | 40.3 N/m | 13.5 N/m | — |
| $m_{eff}$ | 3.80 pg | 1.18 pg | — |
| $f_r$ | 518 kHz (FEM 496) | 539 kHz (FEM est. 517) | $>500$ kHz ✓ |
| $V_{DC\_MAX}$ | 1377 V | **627 V** | $<V_{PI}$ ✓ (margin 1.15) |
| $V_{PI}$ axial / side | 347 / 222 V | 744 / **723 V** | $>200$ V ✓ |
| $x_{PI}$ | 0.40 µm | 3.07 µm | $>1.67$ µm ✓ |
| $a_{0,max}(\omega_r)$, $Q$ | 2.7e4, 77 | **2.8e4 m/s²**, $Q=81$ | $>2\times10^4$ ✓ |
| $|H|_{f\to0}$ peak (achievable) | 8.7 mrad/V (at 222 V) | **73.9 mrad/V** (at 627 V) | maximized |

The achievable sensitivity improves **8.5× over the old design's best stable point**
and **183× over its quoted 10 V operating point**. Figures:
`figures/p3a_sensitivity_vs_V.png`, `figures/p3a_travel_condition.png`,
`figures/p3a_design_space.png` (the ~47 000 feasible designs; the chosen one sits on
the $S$–$V_{DC\_MAX}$ frontier, ~1% below the best raw grid point due to rounding).

*Practical note*: 627 V across a 0.9 µm gap is a $7\times10^8$ V/m field — fine for
the implied cryogenic/vacuum quantum application, but it would require vacuum
packaging; in air, Paschen breakdown would intervene. The voltage is high because
$t=200$ nm fixes a small actuator capacitance; part b's thicker device layer would
also reduce $V_{DC\_MAX}\propto1/\sqrt{t}$.

### a.4 — COMSOL validation of the new design

2D plane-stress model ($t=200$ nm), isotropic Si ($E=170$ GPa, $\nu=0.28$,
$\rho=2330$ kg/m³), same modelling recipe as Part 2 (Domain Terminal at 0 V,
stator electrode at $V_d$, Electromechanical Forces coupling, hyperelastic moving
mesh). Scripts: `comsol/scripts/01_mech_all.py`, `02_em_all.py`.

#### (i) Equivalent spring constant

$1\,\mu$N test force on the spine along $-y$; $k=F/|u_y|$. Mesh convergence
(`data/p3a_q4i_spring_constant.csv`):

| mesh | elements | $k$ (N/m) |
|---|---|---|
| coarse | 4 429 | 12.578 |
| normal | 12 731 | 12.525 |
| fine | 40 916 | 12.501 |

**$k_{FEM}=12.50$ N/m vs analytical $4EtW_b^3/L_b^3=13.45$ N/m: $-7.0\%$.**
Same physics as in Part 2 ($-8.7\%$ there): the analytical model assumes the
plates/column are rigid and the beam ends perfectly guided; in the FEM the
guided ends rotate slightly with the compliant plates, and shear deformation
adds compliance. The model-calibrated prediction used for the design
($0.92\,k_{ana}=12.38$ N/m) is within 1% of the FEM value.

#### (ii) First 5 eigenmodes (2D maps)

`data/p3a_q4ii_eigenfrequencies.csv`, two refinements (0.1% converged):

| mode | $f$ (FEM) | character | map |
|---|---|---|---|
| 1 | **521.7 kHz** | $y$-translation of shuttle (actuation mode), S-shaped beams | `figures/p3a_mode1_522kHz.png` |
| 2 | 2 634.8 kHz | spine in-plane rotation about the column, fingers swing laterally | `figures/p3a_mode2_2635kHz.png` |
| 3 | 2 786.1 kHz | bar in-plane rotation (symmetric partner of mode 2) | `figures/p3a_mode3_2786kHz.png` |
| 4 | 5 519.3 kHz | symmetric in-plane bending of both plate tips | `figures/p3a_mode4_5519kHz.png` |
| 5 | 7 616.5 kHz | higher-order plate bending + beam second modes | `figures/p3a_mode5_7616kHz.png` |

**$f_1=521.7$ kHz $>500$ kHz ✓** (analytic 538.6 kHz, $+3.2\%$ — exactly the
softer FEM spring: $\sqrt{12.50/13.45}=0.964$). The FEM modal mass
$k/\omega_1^2=1.163$ pg agrees with the lumped estimate (1.175 pg, $-1\%$).
Mode separation $f_2/f_1=5.1$: the redesigned short fingers (internal mode
$\sim$32 MHz) removed Part 2's 3.3 MHz finger band from the low spectrum
entirely, and modes 2/3 are the plate-rotation family that the side-stability
constraints deliberately keep stiff. Mode 4 corresponds to the analytic
bar/spine tip estimate (6.8 MHz, cantilever approximation $-19\%$ high — the
tips are loaded by the adjacent beam slots, not free).

#### (iii) Maximum DC displacement at the intended operating voltage

Stationary continuation sweep $V_d=0\to760$ V (steps 10 V, 627 V included).
Data: `data/p3a_q4iii_vsweep.csv`, figures: `figures/p3a_q4iii_vsweep.png`,
`figures/p3a_sensitivity_fem.png`.

**$x(627\,\text{V})_{FEM}=0.919\,\mu$m** toward the comb (away from the waveguide),
vs the analytical optimum $1/(2b)=1.667\,\mu$m: the FEM displacement is **45% below**
the lumped prediction. The sweep data decompose the deficit into two real physical
effects the ideal-comb model misses:

1. **Shallow-engagement (fringe-regime) comb constant.** The force-implied
   $(dC/dx)_{eff}=2k_{FEM}|u_y|/V^2$ is only $\approx71$ pF/m at low bias vs the
   ideal $N\varepsilon_0t/h=114$ pF/m (×0.62). With $d_0=1\,\mu$m the rest overlap
   is only $1.1h$ — the parallel-plate flank field is not yet developed (Part 2
   never saw this: its $d_0=19\,\mu$m comb was deeply engaged). The constant ideal
   $dC/dx$ assumption is simply not valid below engagements of $\sim(2\!-\!3)h$.
2. **Beam stress stiffening (Duffing hardening).** At micron-scale stroke the
   fixed–guided beams stretch ($x/W_b\approx0.5$–0.9): the prestressed
   eigenfrequency data show $f_1$ *rising* from 522 kHz (0 V) to 586 kHz (450 V),
   i.e. the tangent stiffness grows $\approx26\%$, further suppressing $dx/dV$.

Consequence for the operating point: the *true* sensitivity peak shifts to higher
voltage. The FEM-derived low-frequency sensitivity
$|H|(V)=S_\phi(g_0+x)\,dx/dV$ (figure `p3a_sensitivity_fem.png`) is 35.7 mrad/V at
the design's 627 V, and **54.9 mrad/V at 760 V and still rising** — the FEM-true
optimum sits near $x=1/(2b)$, reached at $V\approx870$ V by extrapolation of the
sweep. Because $S(x)\propto\sqrt{x}\,e^{-bx}$ is very flat around the optimum,
operating anywhere in the 750–900 V plateau recovers the design sensitivity to
within a few percent; the redesign goal (peak below pull-in) is preserved since
pull-in moved up even faster (see below).

#### (iv) Maximum DC acceleration without hitting the waveguide

Body load $\rho a$ on all suspended Si ($+y$, toward the waveguide) at $V_d=0$
(worst case: the bias only adds clearance); contact when the bar top closes
$g_0=200$ nm. Linear solve pair + verification solve at $a_{max}$.
Data: `data/p3a_q4iv_max_acceleration.csv`.

**$a_{max,DC}=2.12\times10^6$ m/s² $=2.16\times10^5$ g** (verification residual
$-0.6\%$, the small deficit being the onset of geometric stiffening).
Analytical $k_{FEM}\,g_0/m_{static}=2.10\times10^6$ m/s² (+1.3% ✓); the pure
lumped prediction $k_{ana}g_0/m=2.26\times10^6$ m/s² is 6% high (softer FEM
spring). Note this is the *quasi-static* limit; the binding design constraint was
the resonant-AC value $a_{0,max}(\omega_r)=g_0\omega_r^2/Q=2.8\times10^4$ m/s²
$>2\times10^4$ ✓ (worst-case frequency $f_r$, amplification $Q\approx81$).

#### (v) Step response and rise time

Time-domain analysis with Rayleigh damping (stiffness-proportional,
$\beta_{dK}=1/(2\pi f_1 Q)$, $Q=80.7$ from the Part-1 damping model at the new
geometry, $f_1=521.7$ kHz FEM), 1 µN force step on the spine from rest, 0–300 µs,
rtol $10^{-5}$. Data: `data/p3a_q4v_step_response.csv`, figure:
`figures/p3a_q4v_step_response.png`.

| quantity | FEM | analytical |
|---|---|---|
| final displacement | 79.71 nm | static $F/k_{FEM}$ = 80.0 nm (0.4%) ✓ |
| **rise time (10–90%)** | **0.30 µs** | underdamped $1-\cos\omega_1t$ rise: $1.02/\omega_1=0.31$ µs ✓ |
| **settling time (5%)** | **148.6 µs** | $\ln(20)/(\zeta\omega_1)=147.5$ µs (0.7%) ✓ |
| overshoot | $\approx2\times$ | underdamped $2-\pi/(2Q)$ |

The device is strongly underdamped ($\zeta=1/2Q=6.2\times10^{-3}$): the *rise*
is fast (a fraction of the 1.92 µs oscillation period — the comb force
accelerates the light 1.18 pg shuttle almost instantaneously on the modulation
timescale), but *settling* to 5% takes $\approx78$ oscillation periods. Compared
with Part 2 (158 µs at $f_1=496$ kHz, $Q=83$): $t_s\propto Q/f_1$, and the new
device's slightly higher $f_1$ and slightly lower $Q$ both shorten settling
marginally. For modulation use the ring-down is irrelevant (the device operates
as a resonant or sub-resonant modulator, not a step actuator), but it sets the
recovery time after a shock or bias step.

#### Pull-in verification (extra)

Two FEM checks (beyond the assignment list, since the design hinges on
$V_{DC\_MAX}<V_{PI}$):

1. **Static branch**: the continuation sweep converges through **760 V with no
   fold** — no pull-in anywhere in or near the operating range
   ($760/627=1.21\times$ the design voltage).
2. **Prestressed eigenfrequencies vs bias** (continuation-ramped stationary +
   eigenfrequency at each bias; `comsol/scripts/06_eigcheck_warm.py`,
   `data/p3a_eig_vs_bias.csv`, figure `figures/p3a_eig_softening.png`):

| $V_d$ | $f_1$ (actuation) | $f_2$ (spine rotation = side mode) |
|---|---|---|
| 0 V | 522.2 kHz | 2 636.7 kHz |
| 400 V | 557.4 (+6.7%) | 2 512.7 (−4.7%) |
| 627 V | **675.2 (+29%)** | 2 170.5 (−18%) |
| 760 V | 584.5 | 1 519.5 (−42%) |

$f_1$ is **non-monotonic**: beam stress-stiffening dominates up to ≈650 V
(+29%), then electrostatic softening takes over; $f_1^2\to0$ extrapolates to
≈1100 V. The critical channel is **mode 2 — the spine-rotation/side mode**, whose
$f_2^2\to0$ extrapolation gives **$V_{PI,FEM}\approx876$ V** (mode shape at 760 V:
`figures/p3a_softening_mode2.png`). So:

- the analytically designed side pull-in (723 V, Dunkerley-combined, P2-calibrated)
  is **17% conservative** vs the FEM — expected, since Dunkerley underestimates
  when several channels have comparable strength;
- the *character* of the instability (lateral/side, via the comb's negative
  lateral stiffness — not the axial tip fold) confirms the Part-2 lesson and the
  constraint set used in the optimization;
- ordering as designed: $V_{op}$ (627 V) < $V_{PI,side}$ (876 V) < $V_{PI,axial}$
  ($f_1^2\to0$: ~1100 V) — **all four assignment constraints hold in the FEM**:
  $f_r=521.7$ kHz > 500 kHz ✓, $V_{PI}\approx876$ V > 200 V ✓, $V_{DC\_MAX}<V_{PI}$ ✓
  (28% voltage margin at the design point), $a_{0,max}(\omega_r)=2.8\times10^4$ >
  $2\times10^4$ m/s² ✓.

**Summary of analytic ↔ FEM deviations (the a.4 discussion in one table):**

| quantity | analytic (design) | FEM | cause of deviation |
|---|---|---|---|
| $k$ | 13.45 N/m | 12.50 N/m (−7%) | plate compliance, shear |
| $f_r$ | 538.6 kHz | 521.7 kHz (−3%) | softer $k$ |
| $x(627\,\text{V})$ | 1.667 µm | 0.919 µm (−45%) | fringe-regime $dC/dx$ (×0.62) + stress stiffening (+26% $k_{tan}$) |
| peak $\|H\|$ | 73.9 mrad/V @ 627 V | ≈55 mrad/V @ 760 V (stable, still rising) | same two causes; flat $S(x)$ keeps the loss modest |
| $V_{PI}$ (side) | 723 V | ≈876 V (+17%) | Dunkerley conservatism; stiffening helps stability |
| $a_{max,DC}$ | 2.26×10⁶ m/s² | 2.12×10⁶ (−6%) | softer $k$ |
| $t_s$, rise | 147.5 µs, 0.31 µs | 148.6 µs, 0.30 µs | ✓ (≤1%) |

The two systematic FEM effects — shallow-engagement comb force and beam stress
stiffening — *shift* the sensitivity optimum from 627 V toward the 750–900 V
range without destabilizing it (pull-in moves up by the same stiffening physics).
A second design iteration could recover the low-voltage optimum by deepening the
rest overlap ($d_0\to2.5$–3 µm, restoring ideal $dC/dx$) and/or lengthening the
beams at constant $k$ (reducing $x/W_b$ and thus the Duffing hardening), at a
modest cost in side-stability margin — left as a documented refinement since all
assignment constraints are already met.

---

## Part b) 3D model to mitigate out-of-plane behavior

Device dimensions from **Part 2** (not the part-a redesign), comb fingers removed,
as instructed.

### b.1 — Why $t=200$ nm is problematic

A fixed–guided beam bending **in plane** has $I_{ip}=tW_b^3/12$; bending **out of
plane** only swaps the section dimensions, $I_{oop}=W_b t^3/12$. With four beams:

$$k_{ip}=\frac{4E\,t\,W_b^3}{L_b^3}=40.3\ \text{N/m},\qquad
k_{oop}=\frac{4E\,t^3\,W_b}{L_b^3}=k_{ip}\left(\frac{t}{W_b}\right)^2 .$$

$$\frac{k_{oop}}{k_{ip}}=\left(\frac{0.2}{5}\right)^2=\frac{1}{625}
\;\Rightarrow\; k_{oop}=0.064\ \text{N/m}.$$

Since the modal mass is the same, $f_{oop}=f_{ip}\,(t/W_b)\approx496/25\approx
20$ kHz. Consequences:

- the **true fundamental mode is out-of-plane flapping at ~20 kHz**, not the
  500 kHz actuation mode — ambient acoustic/structural vibration sits exactly in
  that band and is amplified by $Q$;
- the out-of-plane shock tolerance drops by the same factor 625: accelerations of
  only $\sim10^3$ m/s² move the bar by hundreds of nm vertically, modulating the
  evanescent overlap with the waveguide (parasitic amplitude/phase noise) and
  risking stiction to the substrate;
- electrostatic comb levitation (the asymmetric fringe field above the open comb)
  directly drives this soft DOF at $2f_m$;
- a 200 nm plate spanning 100 µm is also vulnerable to residual-stress warping.

### b.2 — Thickness for an in-plane fundamental mode; 3D eigenfrequency analysis

Both $k_{ip}$ and $m$ scale $\propto t$, so $f_{ip}$ is **independent of $t$**, while
$f_{oop}\propto t$. The fundamental mode becomes the in-plane actuation mode when

$$k_{oop}\ge k_{ip}\iff t\ge W_b=5\ \mu\text{m}\quad\text{(naive, beam-only)}$$

(at $t=W_b$ the two modes would be degenerate). This is the *hand* answer; the 3D
FEM below shows it is insufficient, because the lowest out-of-plane mode of the real
structure is a pitch (see-saw) mode that the single-DOF beam argument cannot see —
the actual crossover is at $t\approx7.5$ µm and the final choice is **$t=10$ µm**.

**3D COMSOL eigenfrequency analysis** (swept-mesh extruded model, Part-2 in-plane
dimensions, no comb fingers; `comsol/scripts/03_p3b_3d.py` + `03b_p3b_extend.py`;
data `data/p3b_modes_vs_t.csv`, figure `figures/p3b_modes_vs_t.png`):

| $t$ (µm) | lowest OOP mode | in-plane actuation mode | fundamental |
|---|---|---|---|
| 0.2 | **20.0 kHz** | ~512 kHz | OOP (×26 below!) |
| 1.0 | 97.3 kHz | 512.2 kHz | OOP |
| 5.0 | 379.5 kHz | 513.0 kHz | OOP — naive criterion fails! |
| 6.0 | 430.6 kHz | 513.3 kHz | OOP |
| 8.0 | 522.8 kHz | 514.0 kHz | **in-plane** (1.7% margin) |
| **10.0** | **607.9 kHz** | **514.8 kHz** | **in-plane (18% margin)** ✓ |

Key findings:

1. The in-plane actuation frequency is **independent of $t$** (512–515 kHz for all
   $t$) as predicted — both $k_{ip}$ and $m$ scale with $t$. (FEM 514.8 vs analytic
   531 kHz: −3%, the softer FEM spring.)
2. At $t=200$ nm the lowest mode is out-of-plane at **20.0 kHz — matching the
   analytical $f_{ip}(t/W_b)\approx20$ kHz exactly** and confirming the b.1 hazard.
3. **The naive crossover criterion $t\ge W_b=5$ µm fails in 3D.** The lowest OOP
   mode is *not* the pure flap the beam formula describes, but a **see-saw pitch
   mode** (bar and spine antiphase in $z$, rotation about the beam axis — figure
   `p3b_mode2_608kHz.png`). A lumped pitch estimate — four beams' out-of-plane
   bending stiffness $12EI_{oop}/L_b^3$ acting on the beam-level lever arms
   $\pm(L_m+W_b)/2$ against the plate inertia at $\approx35$ µm radius — gives
   $f_{pitch}\approx0.55$ MHz at $t=10$ µm (FEM 608 kHz; torsion of the beams,
   neglected, makes up the difference). Because the pitch mode rises more slowly
   with $t$ than the flap, it crosses the actuation mode only at
   $t\approx7.5$ µm.
4. **Final thickness: $t=10$ µm** (in-plane fundamental with 18% separation).
   Realistic SOI device-layer thickness; all other dimensions unchanged.

**First 5 3D eigenmodes at $t=10$ µm** (figures `p3b_mode1..5_*.png`):

| mode | $f$ (FEM) | character |
|---|---|---|
| 1 | **514.8 kHz** | in-plane $y$-translation — the actuation mode is now the fundamental ✓ |
| 2 | 607.9 kHz | out-of-plane see-saw pitch (bar/spine antiphase) |
| 3 | 927.9 kHz | out-of-plane flap (whole shuttle translates in $z$) |
| 4 | 1 894.1 kHz | higher out-of-plane (antisymmetric beam-pair/plate bending) |
| 5 | 2 858.9 kHz | in-plane shuttle rotation (t-independent, cf. 2D mode family) |

**Equivalent spring constants at $t=10$ µm** (`data/p3b_spring_constants.csv`):

| direction | FEM | analytic $4EI\cdot12/L_b^3$ | deviation |
|---|---|---|---|
| in-plane $y$ | 1 883 N/m | 2 015 N/m | −6.5% (plate compliance + shear, as in 2D) |
| out-of-plane $z$ (static, at the spine) | 1 641 N/m | 8 059 N/m | **−80%** |

The large out-of-plane discrepancy is the same physics as the mode ordering: a
static $z$-force on the spine plate is offset from the pitch axis, so the
measured compliance is dominated by the *pitch* DOF, not by the pure flap that
the $4Et^3W_b/L_b^3$ formula describes. A modal cross-check on the flap mode
itself works fine: $f_{flap,FEM}=927.9$ kHz vs $f_{ip}(t/W_b)=1\,030$ kHz (−10%,
beam-end rotation + shear at the now-stubby $t/L_b=0.13$ aspect). The lesson for the report:
in 3D, "the" out-of-plane spring constant is mode-specific, and the lowest
out-of-plane mode of a two-plate shuttle is rotational, not translational.

---

## File index

| artifact | contents |
|---|---|
| `analysis/p3a_model.py` | lumped model: $V_{DC\_MAX}$, pull-in (axial + 4 side channels), Q, sensitivities |
| `analysis/p3a_optimize.py` | constrained grid-search optimizer (a.3) |
| `analysis/p3a_final_design.py` | final dims, summary CSVs, analytical figures |
| `comsol/scripts/model_lib_p3.py` | COMSOL geometry/physics builder (2D, new dims) |
| `comsol/scripts/01_mech_all.py` | a.4 (i) spring, (ii) modes, (v) step response |
| `comsol/scripts/01b_build_em.py`, `02_em_all.py` | a.4 (iii) x(V), (iv) accel, pull-in check |
| `comsol/scripts/06_eigcheck_warm.py` | prestressed eig vs bias (continuation ramps) |
| `comsol/scripts/03_p3b_3d.py`, `03b_p3b_extend.py` | part b: 3D model, t-sweep, modes, springs |
| `comsol/scripts/04_make_plots.py`, `05_clean_models.py` | figures from CSVs; shrink .mph files |
| `data/*.csv` | all exported numbers |
| `figures/*` | all figures (analytical + COMSOL exports) |
