# MolMod-education

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

---

## Repository layout

```
MolMod-education/
├── Notebooks/         Jupyter notebooks (.ipynb)  — read/run these
├── src/               standalone GUI scripts (.py) — the source programs
├── requirements.txt   Python dependencies for the notebooks
└── README.md          this file
```

Each notebook is paired with the source script of the same name
(`Notebooks/Potential-well.ipynb` ↔ `src/Potential-well.py`), with two exceptions:
`LJ-ELEC_MD-SoftCore` is a **notebook only** (derived, no script) and `pymoltris3.py` is a
**script only** (no notebook).

---

## Getting started

### Run the notebooks (recommended)

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

### 1. Energy minimisation — finding the lowest-energy arrangement

Minimise the LJ + Coulomb energy of 20 charged particles in a box; compare how different
optimisers reach (different) local minima.

| Notebook | Script | What you learn |
|---|---|---|
| `LJ-ELEC_EM-steepest`  | `LJ-ELEC_EM-steepest.py`  | **Steepest descent** — follow the force downhill |
| `LJ-ELEC_EM-conjugate` | `LJ-ELEC_EM-conjugate.py` | **Conjugate gradient** (Fletcher–Reeves), optional steepest-descent warm-up |
| `LJ-ELEC_EM-simplex`   | `LJ-ELEC_EM-simplex.py`   | **Downhill simplex** (Nelder–Mead) — a derivative-free minimiser |

*Take-away: gradient methods find deeper minima than the simplex in this high-dimensional search.*

### 2. Molecular dynamics — letting the system move in time

Integrate Newton's equations for the same particle system; look at energy conservation,
temperature, periodic boundaries and thermostats.

| Notebook | Script | What you learn |
|---|---|---|
| `LJ-ELEC_MD-Verlet`         | `LJ-ELEC_MD-Verlet.py`         | **Verlet integrator**, kinetic energy & temperature, Maxwell–Boltzmann velocities |
| `LJ-ELEC_MD-Verlet-noPBC`   | `LJ-ELEC_MD-Verlet-noPBC.py`   | The effect of **dropping periodic boundary conditions** (minimum-image convention) |
| `LJ-ELEC_MD-VelocityVerlet` | `LJ-ELEC_MD-VelocityVerlet.py` | **Velocity-Verlet** + **Berendsen thermostat** (NVE vs NVT; temperature is an average) |
| `LJ-ELEC_MD-SoftCore`       | *(notebook only)*              | The **soft core** — capping the r⁻¹² wall so a collision can't blow up the simulation |

### 3. Monte Carlo — sampling configurations by chance

Metropolis Monte Carlo of the particle system: accept/reject random moves to sample the
Boltzmann distribution.

| Notebook | Script | What you learn |
|---|---|---|
| `LJ-ELEC_MMC`        | `LJ-ELEC_MMC.py`        | **Metropolis Monte Carlo** (displacement + charge-swap moves) |
| `LJ-ELEC_MMC-dipole` | `LJ-ELEC_MMC-dipole.py` | Monte Carlo with **oriented dipoles** (rotation moves, dipole–dipole energy, head-to-tail ordering) |

### 4. Thermodynamics — enthalpy vs entropy

| Notebook | Script | What you learn |
|---|---|---|
| `Potential-well` | `Potential-well.py` | A 1D **double-well** Monte Carlo: **well depth = enthalpy**, **well width = entropy**, and temperature as the referee of ΔG = ΔH − TΔS |

### 5. Algorithms & randomness — standalone demos

Small self-contained programs, not the LJ particle system.

| Notebook | Script | What you learn |
|---|---|---|
| `MC-PI-random`  | `MC-PI-random.py`  | **Monte Carlo estimation of π** (darts in a circle, 1/√N convergence) |
| `GA-optimisation` | `GA-optimisation.py` | A **genetic algorithm** evolving a target string (selection, crossover, mutation) |
| `Random-number` | `Random-number.py` | **How random are random numbers?** — a good RNG (Mersenne Twister) vs a deliberately bad one (*uniform ≠ random*) |

### 6. Just for fun — Tetris in PyMOL

| Script | What it is |
|---|---|
| `src/pymoltris3.py` | **Tetris played inside PyMOL** (Python 3 / numpy / PyMOL 3.x). Run `pymol src/pymoltris3.py`, click the 3-D view, then type `start`. Keys: Left/Right = move, **PgUp** = rotate, PgDn = drop. |

---

## Notes

- The notebooks keep the **physics identical** to their GUI scripts; they only replace the Tk event
  loop with a headless run loop and matplotlib visualisation. All notebooks execute end-to-end with
  no errors.
- These are **teaching toys**: small systems, mixed/loose units and modest step counts, chosen to
  make the concepts visible rather than to be production simulation code.

---

## Contributors

Some of the original python scripts were adapted from scripts written by Dr. Patrick Fuchs, University Paris Diderot.
The scripts and their notebooks are originating from the [Bonvin group](https://bonvinlab.org) at Utrecht University.
