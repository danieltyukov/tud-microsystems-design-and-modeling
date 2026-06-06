# Device geometry spec — Part 2 FEM model (2D, top view)

Derived from the P2 assignment figure (full-res `word/media/image2.png` from the docx),
the parameter table, and cross-checks against Part 1's analytical model.

## Coordinate system
- x: horizontal (waveguide propagation direction), y: vertical (motion direction x_m).
- y = 0 at the TOP surface of the modulation bar. Waveguide sits at y = +g₀ (EXCLUDED from model).
- Comb pulls the shuttle DOWN (−y): overlap d increases, modulation gap g = g₀ + |u_y| increases.

## Suspended shuttle (Solid Mechanics domains, doped Si, grounded)
| Part | x-extent (µm) | y-extent (µm) | Size |
|---|---|---|---|
| Modulation bar | −50 … +50 | −30 … 0 | L×W_m = 100×30 |
| Central column | −15 … +15 | −70 … −30 | W_m×L_m = 30 wide × 40 tall |
| Mass plate (comb spine) | −50 … +50 | −100 … −70 | 100×30 (width not tabulated; figure shows = L) |
| Upper beams (×2) | ±15 … ±90 | −35 … −30 | L_b×W_b = 75×5, attach at column top |
| Lower beams (×2) | ±15 … ±90 | −70 … −65 | 75×5, attach at column bottom |
| Moving fingers (×20) | see comb layout | −120 … −100 | L_c×W_c = 20×1 |

Beam outer ends (x = ±90): Fixed Constraint (anchors assumed rigid).

## Comb layout (40 gaps, h = 1 µm)
Pattern left→right: F₃ g M g F₁ g M g F₁ … F₁ g M g F₃
- 20 moving fingers (1 µm), 19 inner fixed (1 µm), 2 end fixed (3 µm — **the P2 change**), 40 gaps (1 µm)
- Total width 85 µm → comb spans x ∈ [−42.5, +42.5]
- Moving fingers: y ∈ [−120, −100] (root at mass plate bottom)
- Fixed fingers: y ∈ [−121, −101] (root at stator surface y = −121); overlap d₀ = 19 ✓
- Tip gaps: moving tip ↔ stator surface = 1 µm; fixed tip ↔ mass bottom = 1 µm (= L_c − d₀) ✓

## Electrostatics (air domain only)
- Air rectangle: x ∈ [−46, +46], y ∈ [−121, −100], minus finger solids.
- Fixed fingers + stator = CAVITIES in air (not modeled as solid): their boundaries = Terminal V_d.
- Shuttle-facing boundaries (mass bottom, moving finger surfaces) = Ground.
- Open sides (x = ±46) = Zero Charge (default).
- Moving Mesh (Deforming Domain) on air; Electromechanics, Solid coupling for Maxwell stress.

## Material (isotropic Si per assignment)
E = 170 GPa, ρ = 2330 kg/m³, ν = 0.28 (assumed — not given; plane-stress k is ν-insensitive).
2D approximation: plane stress, thickness t = 200 nm (correct for thin SOI device layer).

## Mass bookkeeping cross-check
- P1 lumped: ρt(L·W_m + L_m·W_m + N·L_c·W_c) = ρt·5000 µm² = 2.33 pg → f₀ = 662 kHz
- Figure-faithful: ρt(2·L·W_m + L_m·W_m + 20·L_c·W_c) = ρt·7600 µm² = 3.54 pg → f₀(lumped) = 537 kHz
- Discrepancy: P1 omitted the comb-spine plate and counted 40 (not 20) moving fingers.
  This is the main expected source of FEM-vs-P1 difference in Q1–Q3. Discuss in report.

## Ambiguities resolved by judgement (flag in report)
1. Mass-plate width not tabulated → taken = L (figure proportions match the bar exactly).
2. Beam attachment heights → flush with column top/bottom (figure shows slots above/below).
3. ν not given → 0.28; plane stress makes k insensitive to this.
4. Stator bar height irrelevant (only its surface enters ES as boundary).
