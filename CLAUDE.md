# MolMod-education

Teaching material for the Utrecht Molecular Modelling courses. Most programs simulate a small 2-D
system of Lennard-Jones + Coulomb particles. The originals are **Tkinter GUI scripts** (`src/`);
this project ports them to **Jupyter notebooks** (`Notebooks/`) with identical physics, a headless
run loop and matplotlib figures, and adds standalone teaching notebooks.

## Layout

- `src/` — the 13 standalone `.py` GUI scripts (incl. `pymoltris3.py`).
- `Notebooks/` — all `.ipynb` notebooks (a notebook's paired script is `src/<same name>.py`), plus
  three data folders: `HADDOCK3-antibody-antigen-data/`, `GDNP-data/` and `scoring-data/`.
- `Solutions/` — instructor versions (`MD-analysis-GDNP-solutions.ipynb`, `Docking-scoring-solutions.ipynb`,
  `Docking-scoring-optimisation-solutions.ipynb`). Kept **outside
  `Notebooks/`** so the JupyterLite site does not publish them.
- `README.md` — user-facing overview **organised by topic**, with a Colab launch table. Update it
  whenever a notebook is added or its scope changes.
- `requirements.txt` — notebook deps only (`matplotlib`, `numpy`, `scikit-learn` for Docking-scoring-optimisation, jupyter). GUI scripts need only
  stdlib + Tk; `pymoltris3.py` needs PyMOL.
- `jupyterlite/` + `.github/workflows/` — the in-browser site, built with
  `jupyter lite build --contents ../Notebooks` on every push to `main`.

## Tooling and verification

- Python with `nbformat`/`jupyter`/`matplotlib`: `/Users/abonvin/haddock_git/haddock3/.venv-3.14/bin/python3`.
  The default shell is **tcsh**; system python lacks these packages.
- Build or modify notebooks **programmatically with `nbformat`** (or minimal raw-text patches), not
  by hand.
- Verify every change by executing end-to-end and scanning for errors:
  ```
  <venv>/bin/jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=900 \
      --output <scratch>/out.ipynb <notebook>.ipynb
  ```
  then check no output has `output_type == 'error'`. Look at the figures (extract the PNGs) — "no
  errors" is not enough; labels and annotations regularly need repositioning.
- Notebooks with exercises: verify **both paths** — un-filled (must run clean) and filled-in.
- When tuning or quoting numbers, run the **notebook's own code**, never a re-implementation:
  `x += a + b` vs `x = x + a + b` differ in the last ulp and chaotic runs then diverge.
- For numerical refactors: execute before/after, diff stdout and md5 the output PNGs.

## Conventions (all notebooks)

- **No mention of the original Tk GUI** in notebook prose (sliders, `Go` callback, "live view"…).
  The `.py` pairing is documented in the README only.
- **LJ is written `4ε[(σ/r)¹² − (σ/r)⁶]`** with `Sigma` (zero crossing) and `Epsilon` (well depth).
  The `src/*.py` scripts still use the old `ε·Z(Z−1)`, `Rmin`, `Epsilon=25` — same curve
  (`Epsilon_new = 6.25 = 25/4`). The notebook code keeps the multiplication order
  `4*epsilon * u * (u-1)` so floats are **bit-identical** to the old form; the expanded
  `4*epsilon*(u*u-u)` is not, and would change every chaotic trajectory.
- Shared LJ system parameters: `nAtoms=20`, `Radius=25`, `Sigma=2.24*Radius`, `BoxDim=[500,500]`,
  `Epsilon=6.25`, `Dielec=1`, `qat=Radius`, `frac_neg=0.5`, `CutOff=250`, `Seed=100`.
- LJ notebook section order: imports → helpers → energy → forces → method → **Parameters** (split
  tables) → initialisation → run → plots → animation → comparison/take-homes.
- `distsquare = max(distsquare, 1e-6)` floor only in the simplex notebook; the soft core only in
  the velocity-Verlet/SoftCore notebooks. Do not add floors elsewhere.
- **Incidental code cells ship collapsed** via `metadata.jupyter.source_hidden = true` (install/imports,
  pure plotting). Cells students read or edit stay open. Do not hide inputs with CSS — it makes cells
  unrunnable.
- **Pyodide warning filter**: every import cell has, right after `import matplotlib.pyplot as plt`, an
  `if sys.platform == "emscripten":` filter for `"The \w+ parameter as float was deprecated"`.
  JupyterLite's matplotlib 3.10.8 wasm build warns on every text draw; desktop 3.10/3.11 do not.
  Keep it in new notebooks.
- **Exercise mechanism**: stubs return NaN/`None`; a guard (`works()`/`solved()`) detects the
  placeholder; every downstream cell degrades gracefully ("Not yet: …"); worked solutions in folded
  `<details>` blocks. The notebook must execute cleanly un-filled.
- **Markdown: never write "~" for "about"** — two tildes in a paragraph render as strikethrough in
  JupyterLab. Use "≈".
- **Encoding**: follow the file. `LJ-ELEC_EM-steepest.ipynb` is plain UTF-8 (re-saved from
  JupyterLab); most others use `\uXXXX` escapes. `src/LJ-ELEC_EM-conjugate.py` has **CRLF** line
  endings (patch with `newline=''`).
- **Before regenerating a notebook, diff it against HEAD** — hand edits made in Jupyter would be lost.

## Notebooks

### Interactions
- **`LJ-ELEC_Potentials`** (no `.py`): plots the LJ (split into r⁻¹² / r⁻⁶) and Coulomb potentials
  in real units (Å, kcal/mol); no simulation. §16 has two exercises (8-4 vs 12-6 potential; the
  dielectric that balances Coulomb against vdW). **§14 (animated approaching pair) must not use the
  words kinetic energy, velocity-Verlet, thermostat or energy conservation** — students have not met
  them yet; it is told in separation and speed only.

### Energy minimisation (steepest / conjugate / simplex)
- Results at `Seed=100`: SD ≈ −409 in ≈1110 steps, CG ≈ −409 in ≈590, simplex ≈ −187 in ≈620.
  These are local minima: over seeds 100–129 SD gives −312 ± 57, so −409 is the luckiest seed.
- **Conjugate gradient normalises the whole 2N search vector**, like steepest descent (the original
  per-atom normalisation was removed from notebook and script). Do not reintroduce the per-atom
  variant or a comparison with it.
- Simplex keeps the original's quirks (2·nAtoms-vertex simplex, cumulative `n2conv` counter); tuned
  `Simplex_step=300`, `FracShrimp1=0.5`, `FracShrimp2=0.4`, `FracExpend=2.0`, `n2conv=600`.
- `LJ-ELEC_EM-steepest` §13 exercises: tune `drmax/alpha/beta`; steepest *ascent* (blows up — the
  surface is unbounded above); "how local is local" over student-chosen `SEEDS` (`minimize_with_seed`
  is provided; it reads the global `Seed`, which the cell restores). `LJ-ELEC_EM-conjugate` §13:
  steepest-descent warm-up (`numsteep`) — it never helps.
- The steepest notebook's exercise section was built by a script that drops `## **13. Exercises**`
  onward and re-applies hand edits via a `BODY_EDITS` list.

### Molecular dynamics
- **`LJ-ELEC_MD-Verlet`**: position Verlet, shared system params, `timestep=5e-3` (the `.py` uses 1e-3).
- **`LJ-ELEC_MD-Verlet-noPBC`**: with standard params PBC makes **no visible difference** (no wall
  crossings); the notebook shows the effect with controlled demos A/B and the full run C. The `.py`'s
  apparent energy drift was a `calc_temp` indentation bug, now fixed — not the missing PBC.
- **`LJ-ELEC_MD-VelocityVerlet`**: velocity Verlet + Berendsen thermostat (`tauT=0.1`), 50 000-step
  run, soft core `SoftCore = 0.8*Sigma` (a colliding ± pair otherwise explodes at t≈34.4 at any
  timestep). Explains why Etot is not conserved under a thermostat (NVE vs NVT overlay) and that
  temperature is an average over particles.
- **`LJ-ELEC_MD-SoftCore`** (no `.py`): soft core on/off; without it the run **explodes at step
  34442**. This only reproduces because `InitVel` **continues** InitConf's RNG stream — do not reseed it.

### Monte Carlo
- **`LJ-ELEC_MMC`**: displacement + charge-swap moves; rejected moves restore the saved coordinate.
- **`LJ-ELEC_MMC-dipole`**: adds oriented dipoles and rotation moves; the original's ΔE bugs were
  fixed by recomputing the full energy for every move type. Notebook and `.py` default to **neutral**
  particles (`qat=0`) so ordering is purely dipolar; order is measured with a **local** order
  parameter (net polarisation is misleading — dipoles form head-to-tail rings). `mm=1000`.

### Thermodynamics
- **`Potential-well`**: 1-D double-well Metropolis MC — enthalpy = depth, entropy = width. The
  original loop was not Metropolis (retried until acceptance); fixed to one trial per step, counting
  rejected configurations. Defaults hide the entropy effect at 300 K; the notebook scans T and wells.
- **`MD-analysis-GDNP`** — see its own section below.

### Algorithms and randomness
- **`MC-PI-random`**: π from darts + BBP series. **`GA-optimisation`**: genetic algorithm on a string
  (`MUTATION_RATE=0.05`, optional `SWAP_RATE`). **`Random-number`**: Mersenne Twister vs a bad LCG
  (`a=1229, c=1, m=2048`) — passes a 1-D uniformity test, fails lag-1 autocorrelation.

## `MD-analysis-GDNP` — analysis of a real MD simulation

Built from the course exercise `GDNP-analysis.pdf` (untracked, repo root). Analyses the
Cα(1)–Cα(4) distance of the tetrapeptide GDNP (type VIII β-turn) at 280–360 K; cites Fuchs, Bonvin
*et al.*, *Biophys. J.* 90, 2745–2759 (2006). No `.py`.

- **Data**: only `Notebooks/GDNP-data/GDNP-CA-CA-dist-130ns.csv` (6501 frames, 20 ps). European
  export: `;` separator, decimal comma, BOM, CRLF, empty trailing columns — read by a provided numpy
  loader (no pandas). Falls back to downloading from `raw.githubusercontent.com/.../Notebooks/GDNP-data/`
  for Colab. The 70 ns file (same run, first 70 ns at 10 ps) was **deliberately dropped**.
- **Structure**: theory → "Python you will need" primer → data → plots → Ex1 choose `CUTOFF` →
  Ex2 `turn_fraction`, `delta_G` (only the `return` line; hint points to §1 for K) → Ex3 `fit_dH_dS`
  → §8 cutoff scan → Ex4 `running_dG` → §10 block averaging (provided; questions (a)–(g)) → §11 is
  320 K an outlier? → §12 kinetics → §13 take-homes. Four exercises; stubs carry step-by-step
  `# Hints:`; self-checks diagnose common mistakes; `ex2_ok`/`ex3_ok` are strict.
- **Decisions to keep**:
  - Results are reported honestly: at 0.70 nm, ΔG = 6.6–7.6 kJ/mol (> 0 everywhere); weighted fit
    ΔH = +4.7 ± 1.9 kJ/mol, ΔS = −8 ± 6 J/(mol K) — the ΔH/ΔS split is barely resolved even at 130 ns.
    Don't "fix" this into a cleaner story.
  - ΔG **rises** with T but P_turn also **rises** (5.5 → 7.6 % at 0.70 nm): the population follows
    ΔG/RT = ΔH/RT − ΔS/R, and ΔH > 0. Ex2 (d) and Ex3 (d) both say so; don't "correct" it.
  - The cutoff choice is almost purely **entropic** (ΔS −28 → −4 J/(mol K) from 0.60 to 0.80 nm,
    ΔH ≈ 3 kJ/mol).
  - Block averaging is done on **P_turn**, not per-block ΔG (empty blocks give ΔG = ∞). §10 ties the
    minimum block length to Ex4's settling times (113/49/53/74/10 ns at 280…360 K).
  - 320 K is **not** an outlier (1.6σ); dropping it only makes the unweighted fit look better.
  - A two-cutoff *population* count was tried and dropped (no improvement). Two thresholds are used
    only in §12 **kinetics**, where a single cutoff fails (55–70 % of turn visits last one frame).
- **Solutions**: `Solutions/MD-analysis-GDNP-solutions.ipynb` — same cells with solutions filled in,
  an instructor banner, `<details open>`, outputs kept. Both notebooks came from one scratchpad
  `nbformat` script (`--solved` flag, not in the repo); the `.ipynb` files are now the source of
  truth — keep them cell-for-cell aligned.

## `Docking-scoring` + `Docking-scoring-optimisation` — scoring docking models

How well does the HADDOCK score rank docking models, and can it be beaten? Two notebooks, no `.py`.

- **Part 1 `Docking-scoring`** (§1–11): theory, data, term distributions; Ex1 weights + `haddock_score`,
  Ex2 `top_n`, §7 refined, Ex3 `MY_WEIGHTS_*`, §9 random search + cross-validation, §10 Ex4 (questions; `success_rate` provided) on the
  10 BM5 complexes + transfer of the Ex3/§9 weights. Ex4 is questions only: `success_rate` is provided (no
  self-check, no `ex4_ok`). numpy/matplotlib only.
- **Part 2 `Docking-scoring-optimisation`** (§1–8): recap; §3 data + part-1 tools in two collapsed cells (solved
  code from part 1 — keep in sync if part 1's functions change); §4 recomputes `BEST_W`/`CV`; §5 Ex1 `zscore`;
  §6 machine learning (scikit-learn); §7 BM5: transfer revisited (z-scored + ML) and leave-one-complex-out.
- Split on 2026-10-27 from one notebook by a scratchpad script (cells copied, guards and cross-references
  renumbered); the four `.ipynb` files are the source of truth.

- **Data** (`Notebooks/scoring-data/`, downloaded from raw.githubusercontent.com on Colab):
  - case study (part 1 §4–9, reused in part 2): `protein-protein-data-{rigidbody,finalref}.csv` from the user: **1GGR**, E. coli
    IIA(Glc)–HPr (NMR; Wang et al., EMBO J. 2000). Not one of the BM5 `uu_ti5` runs.
    European export: `;`, BOM, CRLF, and the dot is both thousands separator and decimal point — every
    number has 3 decimals, so the loader strips dots and divides by 1000. First column is the rank (header
    mislabelled `HADDOCK-score` in some files).
  - part 1 §10 + part 2 §7: `bm5-{rigidbody,refined}.tsv`, 10 BM5 complexes (1T6B 3WD5 2GAF 2AJF 3LVK 4LW4 5Y9J 3R9A 1HE8 1BGX)
    extracted from tintin `/trinity/login/abonvin/benchmarking/docking_benchmarks/protein_protein/output/uu_ti5/uu_ti5/<PDB>/run1_analysis.tgz`
    (`02_caprieval_analysis` = rigid body, `08_caprieval_analysis` = after flexref + emref). Scenario `uu_ti5`
    = unbound docking with **true-interface (5 Å) restraints** — this is why every optimisation learns to
    up-weight AIR; the notebook says so. Chosen for a mix of successes, scoring failures (5Y9J, 3R9A, 1HE8)
    and one sampling failure (1BGX).
  - A protein–DNA set (1GCC, AtERF1–GCC-box DNA) was tried and **deleted**: its rigid-body set is
    degenerate (914/1000 models one pose). The BM5 complexes replaced it.
- **Weights** (verified by refitting): rigid `vdw 0.01, elec 1, desolv 1, air 0.01, bsa −0.01`; refined
  `vdw 1, elec 0.2, desolv 1, air 0.1, bsa 0`.
- **Main criterion is ranking, not correlation**: `RQ10` = mean over N=1..10 of the mean DockQ of the top N;
  top-N curves; success rates. The random searches optimise RQ10 directly.
- **Decisions / results to keep** (report honestly):
  - case study: in-sample search beats HADDOCK easily; cross-validation: real gain in rigid body (9/10
    splits), overfitting after refinement (3/10).
  - z-scoring within one set cannot reach new rankings (linear); it makes weights readable (w·σ) and
    equal z-weights are a parameter-free baseline (bad in rigid body because of vdW clashes).
  - ML within one complex looks spectacular in rigid body only because near-identical poses put
    siblings in train and test (nearest-neighbour DockQ r = 0.78) — not a docking rule.
  - BM5 leave-one-complex-out: modest gain in rigid body (mostly via AIR), none reliable after
    refinement; in-sample on all ten overestimates. Take-home: optimise on a large benchmark, validate on
    held-out complexes, under realistic restraints.
- **Wide tables**: JupyterLab wraps long printed lines. `print_table`/`print_benchmark` collect their lines
  and call `show_wide(text)` (defined in the imports cell of both parts), which displays a
  `<pre style="white-space: pre; overflow-x: auto">` block (plus a text/plain fallback): a horizontal
  scroll bar instead of wrapping. Keep other printed lines ≲ 95 characters.
- **Top-N curves plot DockQ**: legend entries read "ranked by <score>: mean/best DockQ of top N" — the label
  names the score used for *sorting*, never the plotted quantity. Searched weights are shown with
  `print_weights` (HADDOCK vs search vs ratio, via `show_wide`), defined in part 1 §9 and part 2 §4.
- Runtime: part 1 ≈7 s locally (3 s in Pyodide); part 2 ≈38 s (52 s in Pyodide, the leave-one-complex-out
  cell ≈39 s). scikit-learn is `%pip install`ed in part 2 §6 (Pyodide has it).
- **Solutions**: `Solutions/Docking-scoring-solutions.ipynb` and `Solutions/Docking-scoring-optimisation-solutions.ipynb`
  (outputs kept, banner, `<details open>`) — keep each cell-for-cell aligned with its student version.

## `HADDOCK3-antibody-antigen-lite` — HADDOCK3 tutorial in JupyterLite

- The only version in this repo; the full Colab tutorial lives upstream. Do not reintroduce it or a
  banner pointing to it. Edit the `.ipynb` directly (its build script is not re-runnable here).
- haddock3 itself cannot run in Pyodide (CNS): workflows are described and **pre-computed results**
  analysed from the trimmed `Notebooks/HADDOCK3-antibody-antigen-data/` (the original data host
  blocks CORS). `freesasa`, `biopython`, `pandas` are in Pyodide; `py3Dmol` and `pdb-tools` install
  from PyPI. RCSB downloads work via `pyodide.http.open_url`.
- **Rendering**: `py3Dmol`'s `view.show()` renders nothing in JupyterLab — wrap a full HTML document
  in `<iframe srcdoc>` (inside a `<div>`). Plotly figures are rebuilt from the reports' embedded JSON
  in the same way, loading plotly from a CDN (no plotly package).
- Restraint tools are vendored from haddock3 and checked **byte-identical** in the notebook. Use
  freesasa, not Biopython's ShrakeRupley (it flags extra residues).
- Verify locally with nbconvert **and** in a real browser (build the lite site, serve it, drive with
  Playwright, scrape for error outputs and live iframes).

## PyMolTris (`src/pymoltris3.py`)

Python-3 / numpy / PyMOL-3 port of `pymoltris2.py`. Tkinter was **removed** (a Tk window crashes
PyMOL's Qt GUI on macOS); the game is driven by `cmd.set_key` and typed `cmd.extend` commands
(`start`, `pause`, `rotate`, …). Command callbacks must accept `**kwargs` (PyMOL passes `_self`).
Keys only work after clicking the 3-D view; Up/Down arrows can be swallowed on macOS, so rotation is
also bound to PgUp/Home/CTRL-R and the `rotate` command. Numpy port gotcha: index with
`grid[tuple(idx)]`, not `grid[list]`.

## Porting a remaining script

Un-converted: other `LJ-ELEC_MD*.py` / `LJ-ELEC-MCf*.py` variants. Clone the closest notebook with
`nbformat`, change only theory, the method section and the run cell, keep shared cells identical and
the conventions above, replace the Tk loop by a headless loop recording `traj`/`E_hist`, port Py2→3
faithfully, verify, then update the README and this file.
