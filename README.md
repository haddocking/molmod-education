# MolMod-education

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![Deploy JupyterLite](https://github.com/haddocking/molmod-education/actions/workflows/jupyterlite-deploy.yml/badge.svg)](https://github.com/haddocking/molmod-education/actions/workflows/jupyterlite-deploy.yml)

[![JupyterLite](https://jupyterlite.rtfd.io/en/latest/_static/badge.svg)](https://bonvinlab.org/molmod-education/)
**→ Run the notebooks in your browser: <https://bonvinlab.org/molmod-education/>**

Teaching material for the Utrecht University  **Molecular Modelling (MolMod)** courses. Most of these programs
simulate a small 2D system of **Lennard-Jones + Coulomb** charged particles and let you *watch* how
different algorithms explore or minimise its energy — energy minimisation, molecular dynamics and
Monte Carlo — plus a few standalone demos on algorithms, randomness and thermodynamics.

Each topic comes in **two flavours**:

- **`Notebooks/*.ipynb`** — self-contained **Jupyter notebooks**: theory, a headless run loop and
  matplotlib plots/animations. Best for reading, experimenting and re-running. *(This is the
  recommended way to explore the material.)*
- **`src/*.py`** — the original **interactive Tkinter GUI** programs the notebooks are based on. The
  physics is identical; you drive them with buttons and sliders.

> 🚀 **No installation needed:** all notebooks also run **directly in your browser** via JupyterLite
> at **<https://bonvinlab.org/molmod-education/>** — see [Run in your browser](#run-in-your-browser-no-installation-needed).

---

## Repository layout

```
MolMod-education/
├── Notebooks/         Jupyter notebooks (.ipynb)  — read/run these
│   └── HADDOCK3-antibody-antigen-data/   data for the HADDOCK3 tutorial notebook
├── src/               standalone GUI scripts (.py) — the source programs
├── jupyterlite/       build config for the in-browser JupyterLite site
├── .github/workflows/ GitHub Action that builds & deploys that site
├── requirements.txt   Python dependencies for the notebooks
└── README.md          this file
```

Each notebook is paired with the source script of the same name
(`Notebooks/Potential-well.ipynb` ↔ `src/Potential-well.py`), with three exceptions:
`LJ-ELEC_Potentials` and `LJ-ELEC_MD-SoftCore` are **notebook only** (derived, no script) and
`pymoltris3.py` is a **script only** (no notebook).

---

## Getting started

### Run in your browser (no installation needed)

All notebooks are published as a **[JupyterLite](https://jupyterlite.readthedocs.io/) site** at

**<https://bonvinlab.org/molmod-education/>**

Open the link, pick a notebook from the file browser and run it — everything (Python, matplotlib,
numpy) executes **inside your browser** via [Pyodide](https://pyodide.org); nothing is installed and
nothing is sent to a server. This is the quickest way to get going, and the recommended route for
students.

A few things to know:

- **Your changes are saved in the browser's local storage**, not in the repository. They survive a
  reload on the same browser/machine, but clearing site data wipes them — use *File → Download* to
  keep anything you want to hand in or reuse elsewhere.
- **The first load is slow** (the Python runtime is downloaded once, tens of MB) and each notebook
  pulls its packages on first run; afterwards it is cached.
- The longer runs (the 50 000-step velocity-Verlet MD, the animations) are **noticeably slower** than
  native Python. If you want full speed, install locally as below.

The site is rebuilt and redeployed automatically by a GitHub Action on every push to `main`, so it
always matches the notebooks in `Notebooks/`.

### Run the notebooks locally

```bash
python3 -m venv .venv && source .venv/bin/activate     # optional but recommended
python3 -m pip install -r requirements.txt
jupyter lab            # then open anything in Notebooks/
```

The only third-party packages are **matplotlib**, **numpy** and a **Jupyter** runtime.

### Run the GUI scripts

The scripts need **no pip packages** — just Python's standard library plus **Tkinter** for the GUI:

```bash
python3 src/LJ-ELEC_EM-steepest.py
```

Tkinter ships with CPython but on some systems needs an OS package (not pip):
`brew install python-tk` (macOS) · `sudo apt install python3-tk` (Debian/Ubuntu) ·
`sudo dnf install python3-tkinter` (Fedora).

### Run the PyMOL game

`src/pymoltris3.py` needs a **PyMOL** installation (which brings its own Python + numpy):

```bash
pymol src/pymoltris3.py     # then type 'start' in the PyMOL command line
```

Install PyMOL via `conda install -c conda-forge pymol-open-source` or `brew install pymol`.

---

## Contents by topic

### 1. The interactions — what drives everything else

Before any algorithm: the two pair potentials all the `LJ-ELEC_*` programs are built on, plotted
on their own.

| Notebook | Script | What you learn |
|---|---|---|
| `LJ-ELEC_Potentials` | *(notebook only)* | The **Lennard-Jones** potential split into its **repulsive** (r⁻¹²) and **attractive** (r⁻⁶) components, the **Coulomb** potential, what ε, σ and the dielectric constant do, energy vs force, and why a cutoff is safe for van der Waals but not for electrostatics — in **real force-field units** (Å, kcal/mol, AMBER carbon parameters) |

The **plotting cells start collapsed** (§3), so you see each figure rather than the matplotlib
that drew it. Collapsing is per cell and changes nothing else — click the blue collapser bar to
open or close any one of them, and they execute either way. The physics-bearing cells (parameters,
energy functions, and both exercise stubs) are left open.

It ends with **two graded exercises** (§16) for students to code themselves, each self-checking and
with a folded worked solution:

1. **Build a softer potential** — derive the general *n-m* prefactor, code an **8-4** potential, and
   compare it with the 12-6 at the same ε and σ (wall height, well width, tail range). Why do
   docking and coarse-grained models use softened potentials?
2. **Balance electrostatics against van der Waals** — find the dielectric constant that makes the
   two equal at contact, then identify which real solvent that is. (Spoiler: for a full ±1 e ion
   pair, *none* — it would take ε_r ≈ 800, ten times water.)

*Take-away: read these curves first — the clustering, avoidance and collisions seen in all the
other notebooks are already written into them.*

### 2. Energy minimisation — finding the lowest-energy arrangement

Minimise the LJ + Coulomb energy of 20 charged particles in a box; compare how different
optimisers reach (different) local minima.

| Notebook | Script | What you learn |
|---|---|---|
| `LJ-ELEC_EM-steepest`  | `LJ-ELEC_EM-steepest.py`  | **Steepest descent** — follow the force downhill |
| `LJ-ELEC_EM-conjugate` | `LJ-ELEC_EM-conjugate.py` | **Conjugate gradient** (Fletcher–Reeves), optional steepest-descent warm-up |
| `LJ-ELEC_EM-simplex`   | `LJ-ELEC_EM-simplex.py`   | **Downhill simplex** (Nelder–Mead) — a derivative-free minimiser |

The steepest-descent notebook ends with **three graded exercises** (§13) for students to code
themselves, each self-checking and with a folded worked solution:

1. **Tune the step-size controller** — sweep `drmax`, `alpha` and `beta` for the lowest energy in
   the fewest steps, and learn to tell a *converged* run from a merely *stalled* one.
2. **Run the minimiser backwards** — write a steepest *ascent* and find that there is nothing to
   find: the r⁻¹² wall makes the energy unbounded above, so the run explodes in 21 steps.
3. **How local is "local"?** — minimise ten different random starting configurations and look at
   the spread of final energies (here: 107 kcal/mol, a third of the mean — and the seed the
   notebook uses turns out to be the *best* of the ten). The deepest run ends as one compact
   cluster, the shallowest as three fragments a downhill move can never merge.

The conjugate-gradient notebook ends with **one exercise** (§13) of its own: does a
**steepest-descent warm-up** (`numsteep`) give a deeper minimum, or reach it in fewer steps? The
sweep answers cleanly — the warm-up never changes *where* the run ends up, and past ~20 steps it
only adds cost, sliding from 589 steps back towards the 1109 of pure steepest descent — and ten
starting configurations show the rest of the differences sitting inside their own error bars.
Setting `numsteep` beyond `max_iter` turns the notebook into the steepest-descent one and
reproduces it exactly, which is the cross-check the exercise builds on.

*Take-away: gradient methods find deeper minima than the simplex in this high-dimensional search.*

### 3. Molecular dynamics — letting the system move in time

Integrate Newton's equations for the same particle system; look at energy conservation,
temperature, periodic boundaries and thermostats.

| Notebook | Script | What you learn |
|---|---|---|
| `LJ-ELEC_MD-Verlet`         | `LJ-ELEC_MD-Verlet.py`         | **Verlet integrator**, kinetic energy & temperature, Maxwell–Boltzmann velocities |
| `LJ-ELEC_MD-Verlet-noPBC`   | `LJ-ELEC_MD-Verlet-noPBC.py`   | The effect of **dropping periodic boundary conditions** (minimum-image convention) |
| `LJ-ELEC_MD-VelocityVerlet` | `LJ-ELEC_MD-VelocityVerlet.py` | **Velocity-Verlet** + **Berendsen thermostat** (NVE vs NVT; temperature is an average) |
| `LJ-ELEC_MD-SoftCore`       | *(notebook only)*              | The **soft core** — capping the r⁻¹² wall so a collision can't blow up the simulation |

### 4. Monte Carlo — sampling configurations by chance

Metropolis Monte Carlo of the particle system: accept/reject random moves to sample the
Boltzmann distribution.

| Notebook | Script | What you learn |
|---|---|---|
| `LJ-ELEC_MMC`        | `LJ-ELEC_MMC.py`        | **Metropolis Monte Carlo** (displacement + charge-swap moves) |
| `LJ-ELEC_MMC-dipole` | `LJ-ELEC_MMC-dipole.py` | Monte Carlo with **oriented dipoles** (rotation moves, dipole–dipole energy, head-to-tail ordering) |

### 5. Thermodynamics — enthalpy vs entropy

| Notebook | Script | What you learn |
|---|---|---|
| `Potential-well` | `Potential-well.py` | A 1D **double-well** Monte Carlo: **well depth = enthalpy**, **well width = entropy**, and temperature as the referee of ΔG = ΔH − TΔS |

### 6. Algorithms & randomness — standalone demos

Small self-contained programs, not the LJ particle system.

| Notebook | Script | What you learn |
|---|---|---|
| `MC-PI-random`  | `MC-PI-random.py`  | **Monte Carlo estimation of π** (darts in a circle, 1/√N convergence) |
| `GA-optimisation` | `GA-optimisation.py` | A **genetic algorithm** evolving a target string (selection, crossover, mutation) |
| `Random-number` | `Random-number.py` | **How random are random numbers?** — a good RNG (Mersenne Twister) vs a deliberately bad one (*uniform ≠ random*) |

### 7. Docking — modelling an antibody-antigen complex with HADDOCK3

A full [HADDOCK3](https://www.bonvinlab.org/haddock3/) tutorial: preparing the structures,
turning a predicted paratope and an NMR-mapped epitope into distance restraints, running the
docking workflow and analysing the models. Not an LJ particle system — this one uses the real
docking software.

| Notebook | Runs where | What you learn |
|---|---|---|
| `HADDOCK3-antibody-antigen-lite` | **in your browser** (JupyterLite) | The tutorial runs everything except the haddock3 runs: structure preparation, restraint generation, all analysis, plots and 3D views run in the browser; the docking, scoring and alanine-scanning results are pre-calculated |

`haddock3` drives the [CNS](https://cns-online.org) Fortran engine, which has no browser build —
that is the one thing the `-lite` version cannot do. Everything else is real: the PDB files are
fetched from the RCSB and prepared with `pdb-tools`, and the restraint files are generated by
pure-Python reimplementations of the `haddock3-restraints` sub-commands that give byte-identical
output. Its data (1.6 MB, a trimmed copy of the tutorial bundle) ships in
`Notebooks/HADDOCK3-antibody-antigen-data/`.

### 8. Just for fun — Tetris in PyMOL

| Script | What it is |
|---|---|
| `src/pymoltris3.py` | **Tetris played inside PyMOL** (Python 3 / numpy / PyMOL 3.x). Run `pymol src/pymoltris3.py`, click the 3-D view, then type `start`. Keys: Left/Right = move, **PgUp** = rotate, PgDn = drop. |

---

## Notes

- The notebooks keep the **physics identical** to their GUI scripts; they only replace the Tk event
  loop with a headless run loop and matplotlib visualisation. All notebooks execute end-to-end with
  no errors.
- **Lennard-Jones convention.** All notebooks write the van der Waals energy in the standard form
  **U = 4ε[(σ/r)¹² − (σ/r)⁶]**, with `Sigma` the separation where U = 0 and `Epsilon` the well
  depth (reached at R_min = 2^(1/6)·σ). The `src/*.py` scripts still use the equivalent compact
  form the originals were written in, where the distance parameter is called `Rmin` and the energy
  parameter is 4ε — so `Epsilon = 6.25` in a notebook and `Epsilon = 25` in its script describe the
  **same** curve. Numbers and trajectories are unchanged; only the notation differs.
- The **browser (JupyterLite) version covers the notebooks only**. The `src/*.py` programs need a
  real desktop Python — their Tkinter windows (and PyMOL for the Tetris game) cannot run in the
  browser. The HADDOCK3 tutorial ships **only** in its browser form
  (`HADDOCK3-antibody-antigen-lite`).
- These are **teaching toys**: small systems, mixed/loose units and modest step counts, chosen to
  make the concepts visible rather than to be production simulation code.

---

## Contributors

Some of the original python scripts were adapted from scripts written by Dr. Patrick Fuchs, University Paris Diderot.
The scripts and their notebooks are originating from the [Bonvin group](https://bonvinlab.org) at Utrecht University.
