# Burglence 2.0 — one equation, two methods, three plots

**DIY project, PHY341/745.** Aryan Bandyopadhyay (2411014), NISER Bhubaneswar.

A one-dimensional fluid that steepens, breaks, and then settles into a shock —
solved two ways, checked against exact answers, and Fourier-analysed.

---

## The two keywords, and only two

Sir's list was: self-avoiding walk, adaptive Runge–Kutta, SVD, chaos, Fourier
analysis, 1-D Navier–Stokes.

| Role | Keyword | Where it is used |
|---|---|---|
| **Primary** | **1-D Navier–Stokes** | the Burgers equation `u_t + u u_x = ν u_xx`, solved with finite differences |
| **Secondary** | **Fourier analysis** | the FFT of the shock profile gives `E(k) ∝ k⁻²` |

That is the whole project. Everything else in the folder exists to support
these two.

**Deliberately not used** (this is v1 of the project's problem, not a bug):
spectral methods, Strang splitting beyond one half-step, adaptive Runge–Kutta,
Cole–Hopf integrals, structure functions, fractal exponents. If a sentence in
this project cannot be explained in one line, it was cut.

---

## The physics, in three formulas

```
u_t + u u_x = ν u_xx        0 ≤ x ≤ 2π        u(0,t) = u(2π,t) = 0
```

Each bit of fluid moves at its own speed `u`. A part where `u` decreases
(`u_x < 0`) gets steeper as it moves, and the slope blows up at

```
t_s = -1 / min u₀'(x)   =   1        for u₀ = sin x
```

After that the classical solution has a discontinuity. Viscosity refuses to
allow one, and replaces it with a smooth front of thickness

```
δ = 4ν/U          U = half the jump across the shock
```

whose shape is exactly `u = -U tanh[U(x - x₀)/2ν]`. Both formulas are verified
in Experiment 2 — they are the "replications" this project is built on.

**Why `u = 0` at both ends is exact, not a fudge:** `u₀ = sin x` is odd about
`x = 0`; the equation preserves that symmetry, so the true solution is zero
there for all time.

---

## The three experiments

### 1 — Is the code right?  (`figures/fig01_sanity.png`)

| check | result |
|---|---|
| diffusion only, `u_t = νu_xx` (uses `mylib.pde_heat_1d_crank_nicolson`) | max error **8.0 × 10⁻⁶** |
| convergence before breaking, vs a fine reference | orders **2.02, 2.00, 2.02** |
| explicit scheme at `0.9` and `1.1 × dx²/2ν` | survives / explodes |

### 2 — The shock  (`figures/fig02_shock.png`)

| quantity | theory | measured |
|---|---|---|
| breaking time | 1 | **1.00** (strain follows `1/(1-t)` until then) |
| thickness `δ` at `ν = 0.01` | `4ν/U = 0.0544` | **0.0518** (5 % low) |
| scaling | `δ ∝ ν¹` | **`ν^1.00`** |
| shape | `-U tanh` | max deviation **2 % of U** |
| peak strain | `U²/2ν = 27.1` | **28.4** |

### 3 — Fourier analysis  (`figures/fig03_spectrum.png`)

A sawtooth has Fourier coefficients `c_k ∝ 1/k`, so `E(k) = |c_k|² ∝ k⁻²`.
Measured over `4 ≤ k ≤ 12`, below the dissipative cut-off:

| ν | fitted slope | cut-off `k_cut` |
|---|---|---|
| 0.01 | **−2.06** | 31 |
| 0.005 | **−2.03** | 63 |

Halving `ν` doubles `k_cut` — viscosity decides *where the spectrum ends*,
not what its slope is.

---

## The code: four small files

| file | lines | what it does |
|---|---|---|
| `burgers.py` | 257 | grid, initial condition, the two solvers, the measurements |
| `experiments/e1_sanity.py` | ~150 | accuracy and stability |
| `experiments/e2_shock.py` | ~180 | the four shock panels |
| `experiments/e3_spectrum.py` | ~120 | the two spectral panels |

`burgers.py` has **two solvers**:

* **`explicit`** — forward Euler on `-u u_x + ν u_xx`. Stable only if
  `dt < min(dx/|u|, dx²/2ν)`. The second limit is the killer: halving `dx`
  costs four times the steps.
* **`cn`** — half a step of Crank–Nicolson diffusion, a full step of advection
  (predictor–corrector), half a step of diffusion. Diffusion is implicit, so
  the `dx²/2ν` limit disappears; a run that the explicit scheme cannot do at
  all (dt = 50 × the explicit limit) is routine here.

### What comes from `mylib.py`

| routine | used for |
|---|---|
| `thomas_algorithm` | solving the tridiagonal system inside Crank–Nicolson, every step |
| `poly_fit` | the log–log fits (`δ` vs `ν`, and `E(k)` vs `k`) |
| `pde_heat_1d_crank_nicolson` | the pure-diffusion check in Experiment 1 |

`numpy` is used for array arithmetic and for `np.fft` in Experiment 3.
`matplotlib` for the figures. Nothing else.

---

## Quick start

```bash
pip install -r requirements.txt

python tests/test_burgers.py     # 11 checks, 0.3 s
python run_all.py                # all three experiments, ~45 s
python run_all.py e2             # just one
```

Outputs: `figures/fig01–03.{pdf,png}` and `results/e1–e3.json`.

---

## Documents

| file | contents |
|---|---|
| [`report/report.pdf`](report/report.pdf) | the four-page report |
| [`slides/slides.pdf`](slides/slides.pdf) | the eight-slide deck |
| [`EXPLAINER.md`](EXPLAINER.md) | **everything, in depth** — the story, the physics and numerics derived from scratch, the code line by line, the error analysis, thirty viva questions |
| [`TALK_SCRIPT.md`](TALK_SCRIPT.md) | the 13-minute talk, word for word |
| [`STUDY_GUIDE.md`](STUDY_GUIDE.md) | a six-hour plan to understand all of it, plus viva questions |
| [`REFERENCES.md`](REFERENCES.md) | the three papers and exactly what was taken from each |

---

## Honest limitations

* One equation, one initial condition, one shock — this is a test case, not a
  turbulence simulation.
* The thickness is 5 % below `4ν/U`; the shock spans about 7 grid cells, and a
  discrete maximum slope is a slightly blunt instrument at that resolution.
* The spectrum is measured over one decade of `k` (4–12). Enough to tell
  `k⁻²` from `k⁻³`, not enough for anything subtle.
* Dirichlet boundaries let the wave decay, so the amplitude `U` is not constant
  in time; the `U²/2ν` curve in the figure uses the instantaneous `U(t)`.
