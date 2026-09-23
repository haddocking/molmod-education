# MolMod-scripts — notebook conversion project

Teaching scripts (Utrecht MolMod course) that minimize/simulate 2D Lennard-Jones + Coulomb
charged particles. The originals are **Tkinter GUI** programs; this project ports selected ones
to **Jupyter notebooks** that keep the physics identical but replace the GUI with a headless run
loop plus matplotlib visualisations.

## Directory layout (reorganised 2026-07-20; repo renamed `MolMod-scripts` → `MolMod-education`)

- `README.md` (repo root) — user-facing overview, organised **by topic** (energy minimisation / MD /
  Monte Carlo / thermodynamics / algorithms & randomness / PyMOL Tetris), with a getting-started
  section and a Notebook↔script table per topic. All 14 notebooks + 13 scripts are listed.
  Topic §1 is **the interactions themselves** (`LJ-ELEC_Potentials`), added 2026-09-14 ahead of the
  method topics, which were renumbered 2–7.
- `src/` — the standalone `.py` programs (the 13 GUI scripts, incl. `pymoltris3.py`).
- `Notebooks/` — the `.ipynb` notebooks. Each notebook's paired source `.py` now lives in `src/`
  (not next to the notebook), e.g. `Notebooks/Potential-well.ipynb` ↔ `src/Potential-well.py`.
  One exception to the "notebooks only" rule lives here since 2026-09-23: the HADDOCK3
  antibody-antigen tutorial `HADDOCK3-antibody-antigen-lite.ipynb` (JupyterLite) and its 1.6 MB
  data directory `Notebooks/HADDOCK3-antibody-antigen-data/` — see the HADDOCK3 section at the end.
  The **full Colab tutorial it was derived from (`HADDOCK3-antibody-antigen.ipynb`) is no longer in
  this repo** (removed 2026-09-23, same day): only the lite version ships here.
- `requirements.txt` (repo root) — pip deps to run the **notebooks** (`matplotlib`, `numpy`,
  `jupyterlab`/`notebook`/`ipykernel`). Documents in comments that the **GUI scripts need no pip
  packages** (stdlib + Tkinter; Tk may need an OS package like `python3-tk`) and that
  `pymoltris3.py` needs a separate **PyMOL** install (conda-forge `pymol-open-source` / `brew`).
  Verified: a clean venv from `requirements.txt` executes a notebook with 0 errors and all 13
  `src/*.py` parse.

## Tooling

- Python/Jupyter used for building & verifying: `/Users/abonvin/haddock_git/haddock3/.venv-3.14/bin/python3`
  (has `nbformat`, `matplotlib`, `ipywidgets`, `jupyter`). The system default `tcsh` + bare python
  does **not** have these.
- Build notebooks programmatically with `nbformat` (see the approach below), not by hand.
- Verify every change by executing end-to-end and checking for zero error outputs:
  ```
  <venv>/bin/jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=900 \
      --output /tmp/out.ipynb <notebook>.ipynb
  ```
  then scan the executed copy for `output_type == 'error'`.

## Notebooks produced (all execute with 0 errors)

| Notebook | Source script | Method |
|---|---|---|
| `LJ-ELEC_Potentials.ipynb` | (no `.py`; new 2026-09-14) | **the two pair potentials on their own** — LJ split into its r⁻¹² / r⁻⁶ components, Coulomb, forces, cutoff, + **student exercises** |
| `LJ-ELEC_EM-steepest.ipynb` | `LJ-ELEC_EM-steepest.py` (Py3) | steepest descent |
| `LJ-ELEC_EM-conjugate.ipynb` | `LJ-ELEC_EM-conjugate.py` (Py3) | conjugate gradient (Fletcher–Reeves) + optional SD warm-up (`numsteep`) |
| `LJ-ELEC_EM-simplex.ipynb` | `LJ-ELEC_EM-simplex.py` (**Python 2** source, ported to 3) | downhill simplex (Nelder–Mead) |
| `LJ-ELEC_MD-Verlet.ipynb` | `LJ-ELEC_MD-Verlet.py` (Py3) | molecular dynamics (Verlet integrator) |
| `LJ-ELEC_MD-Verlet-noPBC.ipynb` | `LJ-ELEC_MD-Verlet-noPBC.py` (Py3) | effect of **omitting PBC** (minimum-image convention) in MD |
| `LJ-ELEC_MD-VelocityVerlet.ipynb` | `LJ-ELEC_MD-VelocityVerlet.py` (Py3) | **velocity-Verlet** integrator + **Berendsen weak-coupling thermostat** |
| `LJ-ELEC_MD-SoftCore.ipynb` | (derived from the velocity-Verlet notebook; no `.py`) | explains/illustrates the **soft core** — plots the vdW function ±clamp, MD run ±soft-core |
| `LJ-ELEC_MMC.ipynb` | `LJ-ELEC_MMC.py` (Py3) | Metropolis Monte Carlo (displacement + charge-swap moves) |
| `LJ-ELEC_MMC-dipole.ipynb` | `LJ-ELEC_MMC-dipole.py` (Py3) | Metropolis MC with oriented dipoles (displacement + swap + **rotation**, LJ+Coulomb+dipole–dipole energy) |

> ⚠️ **Reading the per-notebook histories below.** They were written before the 2026-09-14 LJ
> notation refactor, so where they quote notebook parameters as `Rmin=…` / `Epsilon=25` (or `=20`)
> they now read `Sigma=…` / `Epsilon=6.25` (or `5.0`) in the notebooks. The **`.py` scripts still
> use the old names and values**, so statements about `src/*.py` are still accurate as written.
> Same curve, same numbers either way — see the refactor section for the mapping.

The `.py` GUI scripts were also brought to Py3 and made consistent with their notebooks
(shared params, tuned stop thresholds `deltaE=1e-5`/`normFmin=1e-4` for the EM ones; the simplex
`.py` got the full Py2→Py3 port + tuned simplex params). Other `.py` scripts (MC variants,
`LJ-ELEC-MCf*.py`) are **not yet converted**.

### Standalone algorithm demos (not LJ particle systems)

Small self-contained teaching demos, same convert-and-clean workflow (fix `.py`, build notebook
with `nbformat`, verify 0 errors). These do **not** share the LJ system params/structure.

| Notebook | Source script | Illustrates |
|---|---|---|
| `MC-PI-random.ipynb` | `MC-PI-random.py` | Monte Carlo estimation of π (darts in a circle, 1/√N convergence) + Bailey/BBP series |
| `GA-optimisation.ipynb` | `GA-optimisation.py` | genetic algorithm evolving a target string (elitism, crossover, point + optional swap mutation) |
| `Random-number.ipynb` | `Random-number.py` | *how random are random numbers* — good RNG (Mersenne Twister) vs a bad LCG |
| `Potential-well.ipynb` | `Potential-well_Py3.py` | **enthalpy (well depth) vs entropy (well width)** — 1-D double-well Metropolis MC |

Random-number specifics: original is a Tkinter demo dropping random points in a box (each point =
two successive `random()` draws) to eyeball RNG uniformity. `.py` cleaned to Py3 (removed dead
`sqrt`/`randint` imports + unused `dist`). The notebook makes it *discriminating* by contrasting the
**Mersenne Twister** (`random`) against a deliberately **bad LCG** (`a=1229, c=1, m=2048`): plotted
the original's own way (consecutive draws → x,y) the good RNG fills the square while the bad one
collapses onto **diagonal lattice lines** — verified `empty50=0.88` grid occupancy vs `0.04`. The key
teaching sting: the bad LCG still **passes a 1-D uniformity/χ² test** (χ²≈0.64) yet **fails lag-1
autocorrelation** (0.199 vs ~0) — *uniform ≠ random*. Chosen `a=1229,c=1,m=2048` empirically (small
multiplier → few visible lines; large multipliers like `a mod m ≈ 20077` show no 2-D banding; RANDU
`a=65539,c=0,m=2**31` only bands in 3-D). `embed_limit=64` on the animation rc to keep all frames.

GA specifics: `MUTATION_RATE=0.05` (original hard-coded 0.10, very slow tail; 0.05 → converges ~gen
250 at Seed=100); `SEED=100`; optional `SWAP_RATE` (swap two gene positions) — swap alone can never
converge (can't introduce a missing char), only helps *with* point mutation; single-pair swap is best
(swapped-block length is a step-size knob, longer blocks slow a positional target).

MC-PI specifics: `.py` fixed so the π estimate/frac_change/error update **every** step (were inside the
`if inside` block). Notebook: `npoints=100000`, saves one animation frame every 1000 steps (`embed_limit=64`).

Potential-well specifics: `Potential-well_Py3.py` is a Tkinter GUI (already Py3) with a **single
particle** doing **1-D Metropolis MC** in a **static double-well** potential — meant to teach
**enthalpy = well depth vs entropy = well width** (free energy of a flat well `F = -b - kT ln w`;
populations `P1/P2 = (w1/w2)·exp((b1-b2)/kT)` = entropy × enthalpy). Defaults: well 1 **wide+shallow**
(`2a1=200`, `b1=100`), well 2 **narrow+deep** (`2a2=100`, `b2=120`); `cstboltz=8.34e-3` (really R in
kJ/mol/K, ~4× real kB — units are loose/mixed, a teaching toy). **Two issues found & fixed in the
`.py`:** (1) **the MC loop was not proper Metropolis** — an inner `while (not ACCEPTED)` kept proposing
until a move was accepted and only counted *accepted* configs, discarding rejects. This throws away
residence-time weighting and **over-populates the wider well** (measured ~2× off: e.g. b1=100/b2=120/
T=1000 gave frac1≈0.315 vs analytic 0.154). Fixed to **one trial per step, count the config every step**
(rejected move stays put and is counted) → now matches the analytic Boltzmann populations to 3 digits
(verified via the real `Go()` with Tk mocked: 0.154/0.422/0.433/0.667 vs analytic 0.154/0.424/0.434/
0.667). Energy is recomputed from position each step (it's a pure fn of x → no stale-`Ene` bugs);
both x and drawing-height restored on reject; dead `xx=random()` removed; unused globals trimmed.
(2) **The defaults hide the lesson:** depth gap (20) ≫ kT (2.5 at 300 K) so at 300 K the particle sits
in the deep well (frac1≈0.001) and entropy is invisible — you need ~2400 K before the width ratio
competes (that's itself a teaching point: a big ΔH needs lots of thermal energy for entropy to matter).
The **notebook** (`Potential-well.ipynb`, 25 cells, 0 errors, ~4.8 MB w/ animation) keeps the `.py`
defaults but reveals the physics: §8 **temperature scan** (frac1 vs log-T, MC dots vs analytic curve,
approaching the 2:1 width limit), §9 two controlled demos — **equal depths** `b1=b2` → pure entropy
(wide well holds 2/3 at *all* T) and **small gap** `b2=b1+5` → enthalpy↔entropy **crossover**, §10
**depth & width are also knobs** (at *fixed* T=1000: sweep well-2 depth `b2` → population drains
*exponentially*; sweep well-1 width `2a1` → gains *linearly*; MC dots vs analytic, plus a 2-D
analytic heatmap of frac1 over `b2`×`2a1` with the default marked). Also §6 energy/running-fraction
convergence, §7 occupancy histogram over the potential, §11 hopping animation (its own **warm
run at T=1200 K** so the disc actually visits both wells, frac1≈0.22 — the 300 K default is frozen
in the deep well and boring to watch), §12 take-homes.
`run_mc` is a headless single-particle Metropolis loop; `analytic_frac1` gives the exact Boltzmann
fraction for the overlays. Not an LJ system — no shared EM params.

MD notebook specifics: it swaps the minimizer for a **Verlet integrator** (`Verlet`/`Step1`/`CalcVel`)
and adds velocities + kinetic energy/temperature (`calc_temp`) and Maxwell–Boltzmann initial
velocities (`InitVel`). It now uses the **shared EM system parameters** (`Epsilon=25`, `qat=Radius`,
`Seed=100`, plus the common `nAtoms=20`, `Radius=25`, …) so MD starts from the same configuration the
EM notebooks minimize, and adds only the MD-specific `Mass=10`, `Temperature=300`, `timestep`,
`nsteps=2000`. §2 helpers / §3 energy / §4 forces cells are the shared EM cells (with `calc_temp`
appended to §3). No convergence test — records `traj`, per-step `Epot/Ekin/Etot/Temp`; §9 plots
energy conservation + temperature. PBC wraps each moved position **together with its predecessor**
(preserves the Verlet position difference); velocities read **before** wrapping. Result (Seed=100):
Etot ≈ 2800, well conserved (<1% drift over 2000 steps), ⟨T⟩ ≈ 290 K.
Note a remaining `.py`/notebook difference: the `.py` currently has `timestep=1e-3`, the notebook
`5e-3` (timestep is MD-specific, not part of the shared EM set).

MD-noPBC specifics: `LJ-ELEC_MD-Verlet-noPBC.py` is the MD script with the minimum-image line
(`tmp = tmp - SignR(halfbox,tmp-halfbox) - SignR(halfbox,tmp+halfbox)`) **removed** from the energy
(`Calc_Ene2`) and force (`Calc_Force`) loops — meant to show the effect of dropping PBC. **Key finding
(verified headlessly): with the standard params it does NOT actually show any effect** — the 20
particles settle into a caged arrangement (spread across the box, NOT a central droplet) and vibrate
in place (net displacement ~tens of units), and **no particle crosses a wall** (0 edge-crossings), so
the wrap-discontinuity never fires. Minimum-image IS technically active (~37 distant cross-boundary
pairs) but they sit near the cutoff where LJ+Coulomb are weak and, with mixed charges, cancel → PBC
on/off differ by only ~0.1% in energy (both <1% drift, visually identical configs, final positions
differ by <0.4).
The `.py` also had **two problems unrelated to PBC**, now fixed: (1) a `calc_temp` **bug** — the
`v2 = v2 + vx**2 + vy**2` accumulation was **dedented outside the loop**, so only the *last* atom's KE
counted; this (via `InitVel`'s scaling) was the real cause of its ≈25% apparent energy "drift", NOT the
missing PBC. (2) params drifted from the shared set (`Epsilon=2*Radius`, `qat=2*Radius`, `timestep=5e-3`,
`Seed=107`) → **aligned to the PBC `.py`** (`Epsilon=25`, `qat=Radius`, `timestep=1e-3`, `Seed=100`) so
the two scripts differ *only* in the minimum-image line. After the fix noPBC conserves energy (0.2%
drift, ⟨T⟩≈289) like PBC. The **notebook** (`LJ-ELEC_MD-Verlet-noPBC.ipynb`, 17 cells) demonstrates the
effect *properly* with a single `use_pbc` switch on shared energy/force code. It has a **§3 Parameters
section with the two EM-style tables** ("molecular system and its properties" + "molecular dynamics",
the latter incl. the `use_pbc` switch row) added for consistency with the other MD notebooks — since the
demos set params inline, the tables are a reference (Demo C uses the full shared set; A/B fix 1–2
particles). This pushed the demos to **§4 Demo A, §5 Demo B, §6 Demo C, §7 take-homes** (demos are
referenced by letter, so no cross-refs broke). The 3 demos: **(A)** static
pair-energy-vs-position (partner near a wall, slide a partner across → PBC potential is periodic/
continuous and wraps the interaction across the boundary; raw distance drops it discontinuously),
**(B)** a controlled **2-body collision *through* the boundary** (like charges `q=200`, dy=60, v=3,
`nsteps=7000`, straddling the wall → PBC: they repel across the wall, energy swings ≈3700, p0 deflects
≈−36 in y; noPBC: flat energy, cruise past undeflected). NB `nsteps` is capped at 7000 on purpose:
without minimum-image the two particles both drift to the same wall ~step 7500 and would spuriously
collide in *raw* coords, so the run stops before that to keep the noPBC contrast clean. Charges are large
(200) so Coulomb deflects strongly while min-separation stays > `Rmin` (no LJ blow-up); plot shows
y-position vs time, **(C)** the full 20-particle MD both ways → energy curves coincide (0 wall-crossings), **plus a
3-panel config visualisation** (initial / final-PBC / final-noPBC, discs to scale, white +, dark −)
showing PBC and noPBC end in the identical configuration and nothing touches a wall; explains PBC only
matters when the system crosses/straddles the boundary + notes the calc_temp red herring.
To *see* a PBC effect you need the system at the walls (dense fluid / gas filling the box / smaller box
/ net drift); hot or dense many-body regimes explode (LJ singularity, no `distsquare` floor in the MD
`.py`), which is why demos A/B use controlled setups.

Velocity-Verlet + thermostat: `LJ-ELEC_MD-VelocityVerlet.py` (Py3 Tkinter GUI) + `.ipynb` notebook.
The `.py` is `LJ-ELEC_MD-Verlet.py` re-implemented with the **velocity-Verlet** integrator and a
**Berendsen weak-coupling thermostat**. Integrator (self-starting, no `Old_Atom_Coord`/`Step1`/`CalcVel`):
`VelVerlet_Pos` does `r += v*h + 0.5*(F/m)*h²`, then forces are recomputed at the new positions, then
`VelVerlet_Vel` does `v += 0.5*(F_old+F_new)/m*h`; the current forces are kept in a global `Force`
(computed in `Go` whenever `Iterations==0`, i.e. on start/after a reset). Thermostat `Berendsen(vel,
T, T0, h, tau)` scales velocities by `lambda=sqrt(1+(h/tau)(T0/T-1))` (guards T>0 and clamps the sqrt
arg ≥0); applied every step after the velocity update. New param `tauT=0.1` (coupling time, same units
as `timestep`; GUI entry added). PBC wrapping is simplified — only positions are wrapped (velocity-
Verlet stores velocities explicitly, so there is no predecessor to shift, unlike the position-Verlet
script). The **soft core** was later ported from the notebook to the `.py` too: params `SoftCore=0.8*Rmin`
+ `SoftCoreSquare`, and `distsquare=max(distsquare,SoftCoreSquare)` inside the cutoff check of both
`Calc_Ene2` and `Calc_Force` (note the `.py` uses `SoftCoreSquare`, the notebook `SoftCore2`). Verified
headlessly: stable 50000 steps, ⟨T⟩≈307 (a collision spike to ~814 K absorbed+recovered), matching the
notebook. The `.py` is the interactive GUI, so it has no `nsteps`/`report_every`/`nsteps_demo` (those
are notebook run-loop settings). Shares all other params with the plain-Verlet `.py` (`Epsilon=25,
qat=Radius, Seed=100, timestep=1e-3`). **Verified headlessly:** thermostat off (huge tau) → velocity-Verlet conserves energy
(0.78% drift/4000 steps, ⟨T⟩≈289); thermostat on (`tauT=0.1`, target 300) → relaxes to 300 K from hot
(1191 K) and cold (77 K) starts; weaker coupling = slower relaxation (tau 0.02/0.1/0.5 → ~66/339/1973
steps to reach T≈300). The **notebook** (`LJ-ELEC_MD-VelocityVerlet.ipynb`, 29 cells, 0 errors) uses a
single velocity-Verlet integrator + `Berendsen` in a `run_md(..., use_thermostat, report_every)` headless
loop; sections: theory → shared helper/energy/force cells → integrator+thermostat → params (two EM-style
split tables — "molecular system and its properties" + "molecular dynamics and thermostat" — before the
param cell; adds `tauT=0.1`, `use_thermostat`, `SoftCore`, `nsteps`, `report_every`, `nsteps_demo`) →
init → run → energy+temperature plots → config viz + animation → §11 temperature-is-an-average →
§12 thermostat demo (hot start T≈1200→300 for tauT 0.02/0.1/0.5 + NVE reference) → §13 take-homes.
(A "velocity-Verlet vs position-Verlet" comparison section was **removed** — at these params both
Verlet variants gave visually identical energy so it only confused; the position-Verlet helpers and the
`integrator` arg were dropped with it.) **§9 deliberately overlays an NVE reference run** (thermostat off) against the
NVT run: NVE Etot is flat (~1% drift) while its T wanders (233-307 K); NVT holds ⟨T⟩≈298 but its Etot
swings ~24% — because the thermostat pins Ekin while Epot fluctuates, so `Etot=Ekin+Epot` moves with
Epot (here Ekin≈2850 ≫ |Epot|≈200). This is the answer to "why isn't Etot conserved in the run plot?":
a thermostat is *supposed* to exchange energy with a bath; it holds T, not E. velocity-Verlet's own
energy conservation is the flat NVE curve (the §9 grey line). **§11 "Temperature is an average"** addresses a
second common confusion (raised while reviewing the animation): individual particles visibly dart about
while T barely moves. This is correct — `run_md` also records per-particle `Speed`; the section shows
the Maxwell–Boltzmann speed histogram (fastest atom ~2.4× the mean, slowest ~0) and a fastest/mean/
slowest-vs-time plot, making the point that `T ∝ ⟨v²⟩` is a whole-system average, not the fastest
particle (and the thermostat pins that average). The §10 animation title also shows the fastest-atom
speed next to T. (Sections: temperature-is-an-average §11, thermostat demo §12, takeaways §13.)

**Main run is 50000 steps** with `run_md(..., report_every=1000)` printing block-averaged
Epot/Ekin/Etot(±sd)/T(±sd) every 1000 steps (§8), followed by a **3-panel plot** of the run over all
50000 steps: energies (Ekin/Epot/Etot) | **potential-energy split into van der Waals (LJ) `Evdw` and
electrostatic (Coulomb) `Ecoul`** (run_md records both) | temperature (soft-core collision spikes
visible). §9/§11 use that long run, but §9's NVE-vs-NVT
comparison and the §12 demo run at a short `nsteps_demo=4000` (fast + collision-free). **Soft core
added** (`SoftCore=0.8*Rmin`, floored into Calc_Ene2/Calc_Force via `distsquare=max(distsquare,
SoftCore2)`): over a long run an oppositely-charged pair eventually drifts into hard contact and the
r^-12 core blows up. Diagnosed: NOT the integrator (dt 1e-3/5e-4/2.5e-4 all explode at the SAME
physical time t≈34.4 — the collision is under-resolved at any practical dt) and NOT primarily Coulomb
strength (halving qat still explodes ~same time); all-repulsive charges are stable (never touch). The
clamp turns the fatal collision into a recoverable T spike (a block hit ~461 K, thermostat recovers).
Floor tuned: 0.5→still explodes; 0.6/0.7→stable but run heats to ⟨T⟩ 340-540 (a collapsed pair sits
deep and the non-conservative clamp pumps energy); 0.9→⟨T⟩≈300 but fires during normal thermal contact
(perturbs the short NVE demo). **0.8*Rmin** balances both: stable 50k, ⟨T⟩≈305, and §9 NVE drift 0.25%
(clamp silent in the 4000-step window). Clamp is mildly non-conservative (injects energy on a
collision), documented in §6/§9 — why the conservation demos use the short window. Short EM/MD
notebooks omit the floor (2000-4000 steps essentially never collide).

Soft-core explainer notebook: `LJ-ELEC_MD-SoftCore.ipynb` (14 cells, 0 errors, ~200 KB, no `.py`, no
animation) — a standalone explainer for the soft core introduced in the velocity-Verlet notebook, built
on the same engine. `Calc_Ene2`/`Calc_Force` take a `softcore2` arg (`SoftCore**2` = on, `0.0` = off).
§4 **plots the van der Waals function** (LJ energy + radial force vs r, true-vs-soft-cored) — shows the
true LJ diverging (runs off the **clipped linear** y-axis) while the clamp flattens both to a finite
plateau below `SoftCore` (at SoftCore=0.8*Rmin=44.8: U=268, F=85; at 0.5*Rmin=28: U=100800, F=43543, →∞).
It also overlays the pure-repulsion `ε(rmin/r)¹²` reference and marks Rmin (U=0): the "bend" near r=55 a
user may notice is the **attraction↔repulsion crossover** (LJ = repulsion − attraction, not pure r⁻¹²;
at r=55 the two are +31 and −28, net +3), NOT the soft core — the curve is smooth. (Do **not** use a
symlog axis here: its linear↔log switch adds a fake kink that reads as an inflection — use clipped
linear.) §5 **compares MD with
and without the soft core** (`run_md(softcore2, ...)` on the standard 20-particle system, velocity-
Verlet+Berendsen), plotting Etot/T/closest-pair-distance: **OFF explodes at step 34442** (t≈34.4), **ON
stable** ⟨T⟩≈308. **§6 is a side-by-side animation** of the two boxes (`run_md` records `traj`): identical
until t≈34, then the no-soft-core box's colliding pair is flung off-screen (explosion) while the soft-core
box bounces — takeaways renumbered to §7. Note the file is ~6.8 MB with the embedded animation. **Key gotcha:** the explosion only reproduces if `InitVel` *continues* InitConf's RNG
stream (like the GUI script) — reseeding it (as the velocity-Verlet notebook does) changes the initial
velocities and the collision no longer happens within the run; so this notebook's `InitVel` does NOT
reseed. `nsteps=40000` (captures the step-34442 collision).

MMC notebook specifics: **Metropolis Monte Carlo**, derivative-free (simplex-style layout, no force
section). §4 defines `run_mc` (called from §7). Two trial moves: a single-coordinate **displacement**
in `[-deltaRmax, deltaRmax]` (prob `1-frac_swap`) and a two-atom **position/charge swap** (prob
`frac_swap`); Metropolis accept `dE<0 or random()<exp(-dE/kT)` (random drawn only uphill, matching
the source's short-circuit). Rejected moves are undone by **restoring the saved old coordinate**
(cleaner than the source's `-= factor`, which drifts an atom out of the box on a rejected
boundary-crossing move — physically equivalent via periodic image, but the notebook keeps positions
in-box). Running `Ene/ELJ/ECoul` updated **only on acceptance** so histories match the visited
config. Uses the shared EM system params (`nAtoms=20`, `Epsilon=25`, `qat=Radius`, `Seed=100`) +
MC controls `deltaRmax=50`, `frac_swap=0.2`, `Temperature=300`, `cstboltz=8.3502e-3` (J/mol/K, the
MC's own value), `nsteps=20000`. No `distsquare` floor (random moves never *exactly* overlap atoms;
near-overlaps just cost energy and get rejected). Result (Seed=100): E anneals ≈ −112 → fluctuating
~ −340, acceptance ≈ 0.25.

MMC-dipole: variant of MMC where each particle carries an **oriented dipole** (`Atom_Orient[i]`,
moment `mm*(cos,sin)`), adding a **rotation** move and a 2D dipole–dipole term
`U = [m1·m2 − 3(m1·r̂)(m2·r̂)]/r³` (`DipoleEne`, folded into `Calc_EneDip`). The original `.py` was
**buggy**: translation ΔE used `Calc_EneSingle` (no dipole → dipole ignored on moves); rotation ΔE
compared a with-dipole "new" against a no-dipole "old" (wrong sign/magnitude); a debug `print` fired
per pair; single-particle energy fns had a `range(0,len-1)` off-by-one. **Fixed** by computing ΔE
from a **full `Calc_EneDip` recompute for every move type** (translation/swap/rotation all include the
dipole), `Ene=Ene_new` on accept, clean restore-on-reject, debug print removed, buggy single-particle
fns deleted. Made **consistent with `LJ-ELEC_MMC.py`** params (nAtoms=40, Epsilon=20, Rmin=2*Radius,
qat=Radius, deltaRmax=30, frac_swap=0.2, cstboltz=8.3502e-3, Seed=100) + dipole extras `frac_rot=0.25`
and **`mm` tuned 50000→1000** (50000 froze it at this scale; 1000 → acceptance ≈0.22, E anneals
≈ −191 → −979 as dipoles order). The `mm²/r³` scale means dipole strength must track the param scale.
The **notebook** (`LJ-ELEC_MMC-dipole.ipynb`) goes further for teaching clarity: particles are
**neutral by default** (`qat=0`, Coulomb off, `frac_swap=0`, `frac_rot=0.35`) so the ordering is
purely dipolar — set `qat=Radius` to add charges back. It uses a **local** order parameter
`local_order` = ⟨cos(θ_i−θ_j)⟩ over neighbour pairs (rises 0.02→0.36 as chains form), NOT the global
net polarization — which is *misleading* because 2D dipoles order into head-to-tail chains/rings with
~zero net polarization (it actually falls as they order). Dipole arrows drawn via matplotlib `quiver`.
Neutral run (Seed=100): E −109 → −738, local order 0.02→0.36. The **`.py` was also switched to neutral
defaults** (qat=0, frac_swap=0, frac_rot=0.35, grey particles) to match — this required fixing the
`Go` `if qat==0` branch, which used to force `frac_simple_move=1.0` and thereby **disable rotations**;
now it only sets `calc_elec=0`, so translations AND rotations stay active. The GUI `qat` slider still
lets a user add charges back.

## `LJ-ELEC_Potentials.ipynb` — the pair potentials on their own (new 2026-09-14)

Notebook-only (no `.py`). Unlike every other `LJ-ELEC_*` notebook it does **not** simulate: it plots
the two pair potentials that drive all the others. 35 cells, 0 errors, outputs stripped.

**Per-cell code visibility** (2026-09-15). Cells whose Python is incidental — the install/imports
guard and the pure-matplotlib figure cells — ship **collapsed**, via the native
`metadata.jupyter.source_hidden = true` that JupyterLab 4 / Notebook 7 honour. The reader gets a
clickable collapser per cell, *View → Collapse/Expand All Code* for the lot, and — the point — a
collapsed cell still **executes normally**. No dependency, no CSS, no JS.

The build marks them at assemble time from a `COLLAPSED_BY_DEFAULT` tuple of source prefixes in
`cell()`, with an assert that exactly as many cells matched as there are prefixes (so a reworded
first line fails the build instead of silently un-collapsing a cell). Currently 10 of the 15 code
cells; **left open on purpose**: parameters (§4), energy functions (§5), the approaching pair (§14)
and both exercise stubs (§16) — i.e. everything a student reads or edits.

Verified: nbformat valid with the metadata, all 10 collapsed cells report an `execution_count` and
emit their 9 figures after `nbconvert --execute`, and the flag survives execution.

*(Superseded: a global `SHOW_CODE` flag that injected a `<style>` block using CSS `:has()` to hide
every input but its own. It worked, but hiding an input with `display:none` makes the cell
unclickable — you could no longer run it individually. Don't reintroduce it.)*

Section map: 1–2 theory · 3 imports/palette **+ how to show/hide code** · 4 parameters · 5 functions · 6 **the decomposition
figure** · 7 r⁻¹² vs r⁻⁶ ranges · 8 ε/σ + four real atom types · 9 force · 10 kT & Boltzmann ·
11 Coulomb + screening · 12 range & the cutoff · 13 LJ+Coulomb combined · 14 an approaching pair,
animated · 15 take-home · 16 **exercises**. (A 2-D energy-landscape section sat at 14 until
2026-09-14, when it was dropped for simplicity and everything after it shifted down one.)

### Key teaching numbers (all verified against the executed notebook)

| | value |
|---|---|
| C···C well vs thermal energy | ε = 0.1094 = **0.18 kT** at 300 K; Boltzmann peak only **1.20** |
| squeeze to 3.0 Å | +1.04 kcal/mol (≈10× the well depth, for <1 Å) |
| strongest attractive force | 0.077 kcal/(mol·Å) ≈ **5.4 pN**, at 1.24 σ |
| ±1 e pair at contact, vacuum | −87 kcal/mol = 146 kT = **800× the vdW well** |
| same pair in water (ε_r=80) | −1.09 kcal/mol ≈ **1.8 kT** — the textbook salt-bridge strength |
| at a 10 Å cutoff | vdW has **0.31 %** left (0.001 kT); Coulomb **38 %** (55.7 kT vacuum, 0.70 kT water) |
| ±0.1 e partial charges | vacuum: well 9× deeper, like-pair minimum destroyed. Water: ±10 % correction |

### §14 — deliberate vocabulary constraint

The animated approach must **not** use *kinetic energy*, *velocity-Verlet*, *thermostat* or
*energy conservation* — students have not met them at this point in the course. It is told entirely
in **separation and speed**: the pair is "set moving at about 16 Å/ps, a typical speed for a carbon
atom at 300 K", speeds up slightly in the well (16.3 Å/ps), and stops at the **closest approach**
(2.90 Å, well inside σ). The readout shows `t / r / speed`; there is no total-energy line. The
integrator underneath is still velocity-Verlet — just never named. A parenthetical points at the MD
notebooks via a link whose *visible text* avoids the word. **Keep it that way when editing §14.**

### §16 — the student exercises

Two exercises plus a "going further" list. Mechanism: stubs return a `TODO(r)` placeholder (NaN), a
`solved()` guard checks for it, and every downstream cell degrades gracefully — so the notebook
**executes cleanly in its un-filled state** (prints "not yet: replace the TODO(r) placeholders").
The 12-6 curve always plots; the student's 8-4 overlays once it works. Worked solutions live in
folded `<details>` blocks. Verified both paths: un-filled (0 errors) and filled-in (all checks `OK`,
every number matching the solution text).

1. **8-4 vs 12-6.** General *n-m* prefactor `C(n,m) = n/(n−m)·(n/m)^(m/(n−m))`, chosen so the zero
   crossing stays at σ and the depth stays −ε whatever the exponents. Nice result: **C(8,4) = 4
   exactly**, same as C(12,6). The 8-4 minimum moves to 2^(1/4)σ = 4.043 Å, same depth. Softer three
   ways: wall 2.1×/3.1×/4.8× cheaper at 0.9/0.8/0.7 σ; well **60 % wider** at half depth (1.96 vs
   1.19 Å); tail ratio **exactly (r/σ)² → 4, 9, 25, 100** at 2/3/5/10 σ. *(The table must compare the
   `*_attraction` terms, not the totals, or the ratios come out as 3.8/8.9 instead of 4/9.)*
   Part (e) is discussion-only: soft potentials for docking / coarse-graining / early refinement.
2. **Balancing dielectric.** ε_r such that |U_Coul| at Rmin equals the well depth:
   `eps_r = K·q²/(Rmin·ε)`. Scales as **q²**. Results: ±0.1 e → **8.0** (THF 7.6 / DCM 8.9);
   ±0.2 e → 31.8 (methanol 32.7); ±0.5 e → 198.8 and ±1 e → **795** — *no solvent exists*, which is
   the point. Self-check (d): at the balancing ε_r the well lands on −2ε (2.01× deeper, minimum
   pulled in to 3.77 Å). A 16-solvent `SOLVENTS` dict does the nearest-match lookup.

## LJ notation refactor — `ε·Z(Z−1)` → `4ε[(σ/r)¹² − (σ/r)⁶]` (2026-09-14)

Applied to **all 9 `LJ-ELEC_*` simulation notebooks**. Motivation: the `Z` substitution is opaque to
students, *and the old theory text was wrong* — it claimed the minimum was −ε at `Rmin`, but for
`ε·Z(Z−1)` it is **−ε/4 at 2^(1/6)·Rmin**; `Rmin` was really σ, the zero crossing.

What changed, per notebook:

- markdown theory block → the standard 4ε/σ form, with σ and the true minimum stated correctly
  (7 notebooks had a `### Lennard-Jones` block; `MD-VelocityVerlet` and `MD-Verlet-noPBC` had **no**
  statement of the formula at all and were given a one-line reminder + link to `LJ-ELEC_Potentials`)
- `Rmin` → `Sigma` (same value), `rmin_exp6` → `sigma_exp6`, `rmin6` → `sigma6`, `rmin` → `sigma`
- `Epsilon` **rescaled 25.0 → 6.25** (20.0 → 5.0 in `MMC-dipole`, `eps` 25.0 → 6.25 in
  `MD-Verlet-noPBC`) so that `Epsilon` genuinely *is* the well depth
- `LJ2` / `ForceLJ2` bodies rewritten with `u = (sigma/r)^6` and the plain form in the comment;
  the force chain-rule markdown now uses `dedu`/`dudr` with `E = 4ε(u²−u)`
- parameter tables: "position of the LJ energy minimum" → "separation where E_LJ = 0 (minimum at
  2^(1/6)σ)"; suggested `Epsilon` range 1–100 → 0.25–25

**Critical constraint: the arithmetic had to stay bit-identical.** `4*Epsilon_new == Epsilon_old`
exactly in binary (6.25 and ×4 are powers of two), so keeping the **multiplication order** preserved
gives identical floats: `4*epsilon * u * (u-1)` ≡ `epsilon_old * Z * (Z-1)`. The fully expanded form
`4*epsilon*(u*u - u)` is **NOT** bit-identical (≈40 % of evaluations differ in the last ulp) — and
these runs are chaotic, so that would have moved every MD/MC trajectory. `LJ-ELEC_MD-SoftCore`
depends on the no-soft-core run exploding at a *specific* step, so this was not academic.

Verification method (reuse it for any future numerical refactor): execute all 9 notebooks **before**
the change into `base/`, again **after** into `after/`, then diff the concatenated stdout (masking
`0x…` repr addresses) and md5 every output PNG. Result: all 9 numerically identical — same
`EXPLODED at step 34442`, `Final Epot = -409.45`, `Final E = -348.96` — and 21/22 figures
byte-identical (22 of 23). The one differing figure is SoftCore's LJ energy/force panel, where only the legend
text changed (`ε(r_min/r)¹²` → `4ε(σ/r)¹²`); curves confirmed pixel-identical by eye.

Gotcha found this way: `SoftCore`'s standalone plotting helpers (`F_lj`, and the "repulsion only"
curve) take `Epsilon` *directly* rather than via `LJ2`, so they needed their own coefficient fix
(`6*Epsilon` → `24*Epsilon`, `Epsilon*(Sigma/r)**12` → `4*Epsilon*…`). The diff caught it as an
exact factor-of-4 error in the printed force. **Grep for arithmetic uses of `Epsilon`/`eps` outside
`LJ2`/`ForceLJ2`/`Calc_*` when touching this again.**

The `src/*.py` GUI scripts were **not** converted — they still use `ε·Z(Z−1)`, `Rmin`, `Epsilon=25`.
Physics is unchanged so the README's "physics identical" pairing claim still holds, but the notation
now differs between a notebook and its script; the README Notes section documents this explicitly.

## Notebook structure (consistent across all three)

Section order — note **Parameters sits just before Initialisation & Run** (an explore-and-rerun
block), with reusable functions defined above it:

1. Imports — merged cell: Colab-friendly install-guard (`importlib.util.find_spec` → pip) + imports
   + `rc('animation', html='jshtml')` + `%matplotlib inline`. Only 3rd-party dep is **matplotlib**
   (stdlib `math`/`random` otherwise).
2. Helper functions (`dist`, `SignR` minimum-image, `charge_color`)
3. Energy functions (LJ + Coulomb from squared distance, PBC nearest-image)
4. Force functions (EM only) / Simplex minimizer (simplex)
5. The minimizer(s)
6. **Parameters** (split tables: "System and its properties" + "Minimizer")
7. Initialisation (`InitConf`, then builds `Atom_Coord`)
8. Run the minimization (headless loop replacing the Tk `Go` callback; records `traj`, `E_hist`)
9. Energy convergence plot
10. Visualise (initial vs minimized, white=+, dark=−)
11. Animation (`FuncAnimation`, `jshtml`; `stride` caps ~120 frames)
12. Comparison of the three minimizers (shared markdown note)

(Simplex has no Force section, so its numbers are one lower.)

## Shared conventions / decisions

- **System parameters are identical** across all three (for comparison): `nAtoms=20`, `Radius=25`,
  `Sigma=2.24*Radius`, `BoxDim=[500,500]`, `Epsilon=6.25`, `Dielec=1`, `qat=Radius`, `frac_neg=0.5`,
  `CutOff=250`, `Seed=100`. `InitConf` seeds with `Seed` for reproducibility. (`Sigma`/`Epsilon`
  were `Rmin`/`25` before the 2026-09-14 notation refactor — see below; same curve, same numbers.)
- Convergence threshold is named **`deltaE`** in all three (simplex's original `cc1` was renamed for
  consistency).
- **LJ is written as `4ε[(σ/r)¹² − (σ/r)⁶]` everywhere** (notation refactor 2026-09-14, see the
  dedicated section below). The old `E = ε·Z(Z−1)` with `Z = (r_min/r)⁶` is gone from the notebooks;
  the `src/*.py` scripts still use it.
- **Divide-by-zero floor** `distsquare = max(distsquare, 1e-6)` exists **only in `simplex.ipynb`**
  (its moves can overlap atoms); the EM notebooks deliberately omit it (forces keep atoms apart) —
  keep it that way for fidelity.
- Simplex keeps the original's quirks faithfully: compact simplex of `2*nAtoms` vertices (one per
  displaced coordinate) and a **cumulative** (not consecutive) `n2conv` convergence counter.

## Tuning results (with each notebook's current default parameters, Seed=100)

| Method | Final energy | Steps | Key tuned params |
|---|---|---|---|
| Steepest descent | ≈ −409 | ~1100 | `deltaE=1e-5`, `normFmin=1e-4` (tightened from 1e-3; descent params left default) |
| Conjugate gradient | ≈ −383 | ~970 | same stop thresholds; `numsteep=0` (warm-up hurt) |
| Simplex | ≈ −187 | ~620 | `Simplex_step=300`, `FracShrimp1=0.5`, `FracShrimp2=0.4`, `FracExpend=2.0`, `n2conv=600` |

Findings: for the EM methods, aggressive step-size growth finds *shallower* basins — the real win
was tighter stop thresholds (defaults halted early on a plateau). The simplex scales poorly to the
40-D search space and stays far shallower even when tuned (illustrates gradient methods winning in
high dimensions). These are **local** minima; ordering can shift with the seed.

Caveat learned the hard way: a standalone re-implementation of the CG loop **diverged** from the
notebook (predicted −408 vs actual −383). When tuning, sweep parameters against the **notebook's own
code** (exec its cells), not a paraphrase.

## Converting a remaining script (checklist)

The un-converted scripts (`LJ-ELEC_MD*.py`, `LJ-ELEC-MCf*.py`, …) are the same family: a Tkinter GUI
around LJ+Coulomb particles with a different integrator/sampler in place of the minimizer. To port
one, reuse an existing notebook as the template:

1. **Read the `.py`** and identify what differs from the EM notebooks — usually only the core loop
   (the Tk `Go` callback) and a few parameters. Energy/force/`InitConf`/`SignR`/`charge_color` are
   nearly identical (watch for Python-2 syntax: `print` statements, `Tkinter`/`Canvas`, `xrange`).
2. **Clone the closest notebook** with `nbformat` (steepest for a stepping method; simplex for a
   class-based sampler) and mutate only: title/theory, the minimizer/integrator section, and the
   Run cell. Keep all shared cells byte-identical.
3. Keep the **section order** and the **shared conventions** above (merged imports cell, `deltaE`
   name, split Parameter tables, Parameters-before-Init/Run, same system params + `Seed=100`).
4. Replace the Tk event loop with a headless `for step in range(max_iter)` loop that records `traj`
   and `E_hist`; reuse the §9–11 plotting/animation cells unchanged.
5. Port Python-2 → 3 faithfully; only add safety (e.g. a `distsquare` floor) where the dynamics can
   actually hit it, and note it.
6. **Verify** with `nbconvert --execute` (0 error outputs) and, if tuning, sweep against the
   notebook's own code.
7. Append the shared **"Comparison of the … minimizers"** note if it fits, and update this file.

## PyMolTris (Tetris in PyMOL) — Python 2 → 3 port

`pymoltris2.py` is an old **Python-2** script (game "version 2", not Python 2 in the name) that plays
Tetris inside PyMOL. `pymoltris3.py` is its **Python-3 / numpy / PyMOL-3.x port** (original preserved).
Run: `pymol pymoltris3.py` (opens a Tkinter control panel; click Start; Left/Right move, Up/Down rotate,
Ctrl drop, p pause). Tested against **PyMOL 3.1 / Python 3.14 / numpy 2.4.3** (`/opt/homebrew/bin/pymol`,
runs headless via `pymol -cq`).

Port changes (done via a transformer script with count-checked string replacements):
- **Py2→3:** `Tkinter`/`tkFont`/`tkFileDialog`/`tkMessageBox`/`cPickle` → `tkinter`/`tkinter.font`/`…filedialog`/`…messagebox`/`pickle`; removed debug `print` statements; menu builders using `exec("…lambda…")` → real closures/`getattr`; integer division `/2` → `//2` for grid indices; `Event.isSet()` → `is_set()`; `askquestion` result `"YES"` → `"yes"`.
- **Numeric → numpy:** `from Numeric import *` removed, `import numpy as np` added; `zeros(area)` → `np.zeros(area, dtype=int)` (int is essential — grid holds resi numbers/indices); `sum(grid,axis)` → `np.sum`; `take(...)` → `np.take`; **`grid[list]` element access → `grid[tuple(list)]`** (numpy treats a list index as fancy row-indexing, Numeric treated it as element access); row truthiness `if grid[k]:` → `if grid[k].any():`; `if max(max(grid)):` → `if grid.any():`; bare `cos(` (was Numeric's) → `math.cos(`.
- **PyMOL 3.x API:** `from pymol.cgo import cyl_text` (the game's text CGO helper); dropped bare `import _cmd` and its `_cmd.get_atom_coords`/`_cmd.load_png` (use `cmd.get_atom_coords`); **gutted `Checkup`** (removed the obsolete `pymol_argv`-based "run under PyMOL" detection — `pymol_argv` no longer exists in PyMOL 3.1 and would have misfired — plus the dead nmr.chem.uu.nl version-check/auto-update network code); removed the dead splash-PNG download in `Fall.run()`, keeping the CGO "Click on Start" fallback.
- **GUI: Tkinter removed, replaced with native PyMOL controls.** The original built a `tkinter.Tk()`
  control panel. On macOS a Tk window **cannot coexist with PyMOL's default Qt GUI** — Tk's Cocoa
  `TkpGetColor` sends `-[QNSApplication macOSVersion]`, an unrecognized selector → the whole app
  **hard-crashes** (`NSInvalidArgumentException`, `Abort trap: 6`). So a second pass (`port2.py`) **removed
  Tkinter entirely** (dropped the `import tkinter/…font/…filedialog/…messagebox`, the control-panel window,
  menus, buttons, checkbuttons, speed scale, `focout`, and the `bind_all` key bindings) and drives the game
  through **`cmd.set_key`** (Left/Right = move, Up/Down = rotate, PgDn/Insert/CTRL-D = drop) + **`cmd.extend`
  commands** typed in the PyMOL command line (`start`, `pause`, `drop`, `savegame`/`loadgame [file]`,
  `quitgame`, `periodic`, `shownext`, `moreblocks`, `blockcolor <name>`, `blockstyle <n>`, jokers
  `atwa`/`collapse`/`turnaround`/`destroyc`/`destroyr`, `tetris_help`). IntVar checkbuttons → plain int
  attrs (`self.per`/`shn`/`mbl`); button text → `self._status()` prints; the Start button's rebindable
  command → `self._start_action` (PoP normally, `restart` after game over); `_tkM.*` dialogs → prints;
  file dialogs → a fixed `pymoltris.pmt` (pickle, binary mode) or a filename argument. **Now runs with the
  plain `pymol pymoltris3.py`** (Qt) — verified: launches with **no crash/traceback**, prints its help,
  registers the `start` command + key handlers, and the core logic (build_field/newblock/cyl_text/grid
  ops) still passes headlessly. Interactive gameplay (actual keypresses/rendering) still needs a display
  and wasn't click-tested.
- **`_self` gotcha:** PyMOL invokes `cmd.extend` commands with a `_self=` keyword, so every command
  callback must swallow it — the lambdas are `lambda *a, **k: …` (and `lambda arg, **k: …`). Without
  `**k` every typed command (`start`, `pause`, …) raised `TypeError: unexpected keyword argument '_self'`
  and did nothing. (`cmd.set_key` callbacks are called with no args, so the movement handlers `def
  left(self, event=None)` are fine as-is.) Verified with `cmd.do("start")` → `_AS=1`, `_PAUSE=0`, game
  thread launched. There is **no Start button** — the player types `start` in the PyMOL command line; the
  on-screen splash text was changed from "Click on Start" to "Type  start".
- **Keyboard focus gotcha:** the movement keys are correctly bound (`set_key` invokes `fn()` with no args;
  verified via `internal._invoke_key("left"/"right"/"up"/"down"/"pgdn")` that the handlers move/rotate/drop
  the block — the code path is fine). But PyMOL's Qt 3-D viewer uses `Qt.ClickFocus`
  (`pmg_qt/pymol_gl_widget.py`), so it only receives keystrokes **after you click on the 3-D display**.
  After typing `start` the focus is on the command line, so keys seem dead until you click the viewer.
  The help text and header docstring now say "first CLICK on the 3-D display". (Up/Down do nothing on the
  square piece — it is in `norot` — that is correct, not a bug.)
- **Up/Down arrows unreliable on macOS:** a user reported Left/Right worked but Up/Down (rotation) did not,
  despite the code being verified correct (every non-square piece rotates via `internal._invoke_key("up")`,
  and 100→left/101→up/102→right/103→down all map to the handlers). The vertical arrows appear to be
  swallowed by macOS/Qt before reaching PyMOL — not fixable in the script. **Mitigation:** rotation is now
  bound to **several keys** (`up`, `pgup`, `home`, `CTRL-R`; other way `down`, `end`, `CTRL-E`) and exposed
  as typed commands **`rotate`** / `rotback` — verified PgUp/Home/`rotate` all rotate. So use **PgUp** (or
  type `rotate`) if the Up arrow does nothing.

**Verification (headless, `pymol -cq`):** compiles under Py3.14; full script loads with no code errors; a
logic harness (`exec` the class defs, inject `cmd`/`cgo`/`plain`/`numpy`) confirmed numpy grid ops
(tuple-index/`.any()`/`take`/`sum`), `Fall()`→`build_field` (6 sphere atoms), `newblock` (4 `_block`
atoms), `cyl_text`/`load_cgo`, and `checkunderneath` all work. Interactive gameplay (Tk responsiveness +
rendering) needs a display and was not verified.

## HADDOCK3 antibody-antigen tutorial — Colab notebook → JupyterLite (2026-09-23)

`Notebooks/HADDOCK3-antibody-antigen-lite.ipynb` is an **in-browser (JupyterLite/Pyodide)** HADDOCK3
tutorial — not part of the LJ family. It is **the only version in this repo**: it was derived from
the upstream **Colab** tutorial (`HADDOCK3-antibody-antigen.ipynb`, which `pip install`s haddock3,
`wget`s a 54 MB data bundle, runs `!haddock3 …` docking / scoring / alanine-scanning workflows and
visualises with py3Dmol + plotly), and that full notebook was **removed from this repo on
2026-09-23**, the same day both were added. It lives upstream with the HADDOCK3 tutorials.

**Consequences for future edits:**

- Do **not** reintroduce the Colab notebook or a banner cell pointing between the two — the "pair"
  no longer exists here, and the lite notebook must stand alone (its intro should not assume the
  reader has the full version next to it).
- The build script that produced the lite notebook reused the Colab notebook's markdown **by index**
  (skipping the Colab-only banner cell). That script is therefore **not re-runnable from this repo**;
  the lite `.ipynb` is now the source of truth — edit it directly (via `nbformat`) rather than
  rebuilding from an upstream copy.

### What can and cannot run in Pyodide (all checked against the v314.0.6 `pyodide-lock.json`)

- **haddock3 itself: never.** It drives the **CNS** Fortran engine; there is no wasm build. So
  `haddock3 <cfg>`, `haddock3-score` and the `alascan` module are *described* in the lite notebook
  and their **pre-calculated results** analysed. There is no shell either, so every `!command` had
  to go.
- **But `freesasa` AND `biopython` AND `pandas` ARE in the Pyodide distribution** (2.2.1 / 1.87 /
  3.0.2) — that is what makes the restraint sections genuinely runnable. `py3Dmol` and `pdb-tools`
  are pure-Python wheels installed from PyPI. `plotly` is **not** needed at all (see the iframe
  trick). The install cell is a plain `%pip install -q py3Dmol pdb-tools biopython freesasa pandas`,
  the same idiom as the other notebooks in this repo.
- **The tutorial data cannot be downloaded from the browser**: `surfdrive.surf.nl` serves the zip
  with **no `Access-Control-Allow-Origin` header**, so the fetch is CORS-blocked. Fixed by shipping
  a **1.6 MB trimmed subset** as `Notebooks/HADDOCK3-antibody-antigen-data/` (35 files: the
  pre-processed PDBs, restraints, the 3 workflow cfgs, the shell scripts, and only the run outputs
  the notebook reads). `jupyter lite build --contents ../Notebooks` copies it into the site, and the
  Pyodide kernel reads it straight off the JupyterLite drive — **verified in a real browser**: the
  http server log shows `GET /files/HADDOCK3-antibody-antigen-data/…` for all 22 files the run touches.
- **`files.rcsb.org` DOES send `access-control-allow-origin: *`**, so `pdb_fetch` works in the
  browser — via `pyodide.http.open_url` (urllib does not work in Pyodide), falling back to `urllib`
  off-browser. Use the plain `.pdb` URL, not `.pdb.gz`.

### The iframe trick (the key rendering insight)

**`py3Dmol`'s `view.show()` renders nothing in JupyterLab 4 / JupyterLite.** It publishes raw
`<script>` in `text/html` plus an `application/3dmoljs_load.v0` mime type that only **Colab** has a
handler for, and JupyterLab does not execute scripts that come out of a cell.

Fix, used for *everything* interactive in this notebook: wrap a complete HTML document in
`<iframe srcdoc="…">` — the browser builds the iframe's own document and **its** scripts run
normally. `show_html()` in the toolbox cell does this (`html.escape(doc, quote=True)`), and
`show3d(view)` feeds it `view.write_html()` (the public API; `_make_html()` is private).

The same trick removes the plotly dependency: the haddock3 html reports embed their figures as
`<script id="data1|data2|datatable2" type="application/json">{data, layout}</script>` and load
plotly from a CDN. `show_plotly()` pulls that JSON out and calls `Plotly.newPlot` in a tiny
generated page that loads `cdn.plot.ly/plotly-3.0.1.min.js` — no `plotly` Python package, no
`iplot`, no jupyterlab extension. `datatable2` is not a figure: it is read with `pd.json_normalize`.

Two gotchas: **wrap the iframe in a `<div>`**, or IPython prints `UserWarning: Consider using
IPython.display.IFrame instead` on every viewer (it sniffs `content.startswith('<iframe')`); and
escape `</` to `<\/` in the JSON payload so it cannot close the `<script>` early.

### Pure-Python `haddock3-restraints` (verified byte-identical)

The restraint sub-commands are plain Python, so they were reproduced from
`haddock3/libs/librestraints.py` (Apache-2.0, same licence as this repo) into one folded cell:
`passive_from_active` (+ `get_surface_resids`, freesasa + the NACCESS `REL_ASA` bb/sc table),
`parse_actpass_file`, `active_passive_to_ambig`, `read_structure`/`get_bodies`/`build_restraints`/
`generate_tbl` (= `restrain_bodies`), plus a light `validate_tbl` (the `--quick` parenthesis/quote
check + "every assign ends in 3 numbers"; the full 234-line selection parser was **not** vendored).
`align_full` is vendored from `haddock3/libs/libnotebooks.py` (Bio.PDB `Superimposer`).

Each one is checked **in the notebook itself** against the shipped reference file and all print
`True`:

| step | check |
|---|---|
| `pdb-tools` pipelines | regenerated `4G6K_clean.pdb` and `4I1B_clean.pdb` byte-identical to `pdbs/` |
| `passive_from_active` | `3 24 46 47 48 50 66 76 77 79 80 82 86 87 88 91 93 95 118 119 120` — exact |
| `active_passive_to_ambig` | 43 restraints, byte-identical to `restraints/ambig-paratope-NMR-epitope.tbl` |
| `restrain_bodies` | byte-identical (`random.seed(917)` inside haddock3 makes it reproducible) |
| traceback lookup | AF2 ranks 86/90/91/92/93, AF3 40/81/87/88/89 — matches the tutorial text |

**`freesasa` vs Biopython's `ShrakeRupley`**: the Biopython fallback is *not* equivalent — it flags
three extra borderline residues (85, 102, 134 at 16–20 % relative ASA against the 15 % cutoff).
Since freesasa is in Pyodide, use it; do not "simplify" this to Biopython.

**pdb-tools API notes** (the modules expose `run()` generators, far cleaner than driving `main()`):
`pdb_selres.run(fh, residue_range)` wants a **set of ints**, not `(lo, hi)` tuples; `pdb_merge.run`
wants **file-like objects** (it calls `.close()`), so pass `StringIO`, not iterators.

### haddock3-score numbers quoted in the notebook

Generated locally with **haddock3 2026.7.0** (the tutorial pins 2025.9.0):
`4G6M_matched.pdb` → **-149.84**, `4G6M_matched_S150W.pdb` → **-170.26** (S150W is the enriching
mutation PROT-ON proposes). The cell performs the mutation (a `str.replace` of `SER A 150`, 6 atom
lines) and prints these as the reference output.

### Verification (do both for any future change)

1. `nbconvert --execute` locally (needs py3Dmol, pdb-tools, biopython, freesasa, pandas) → **0
   errors**, and every "identical to the provided …" line prints `True`.
2. **In a real browser**, because the rendering path is the whole point: build the site with
   `jupyter lite build --contents ../Notebooks --lite-dir jupyterlite` (note `--contents` is
   resolved **relative to `--lite-dir`**), serve it with `http.server`, then drive it with
   Playwright (`channel="chrome"` avoids downloading a browser) — Run All Cells, then scrape
   `[data-mime-type="application/vnd.jupyter.error"]` and check the iframes for `typeof Plotly` /
   `typeof $3Dmol`.

   Result of that run (2026-09-23): **0 error outputs**, all 33 code cells executed, `PROJECT_DIR`
   resolved to `/drive/HADDOCK3-antibody-antigen-data`, the RCSB download succeeded and all four
   "identical to the provided …" checks printed `True`; 12 iframes alive — 8 with `$3Dmol` and a
   WebGL `<canvas>`, 4 with `Plotly` and 3 `svg.main-svg` each (CAPRI plots, distributions, chord
   chart, alanine scan) — plus 3 pandas tables, and the zip written to
   `/drive/HADDOCK3-antibody-antigen-output.zip`.

   Two gotchas when writing such a driver: JupyterLab gives **folded cells no `.jp-InputPrompt`
   text**, so "everything finished" is *no prompt contains `*`*, not `done == total`; and pipe the
   script's stdout straight to a file (`python -u … > log`), since a `| tail` swallows everything
   if the run is killed.
