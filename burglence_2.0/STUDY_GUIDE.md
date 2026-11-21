# Six hours to understand this project

You do not need to know any of this before you start. Everything below is
either first-year calculus or something you have already coded in this course.

---

## Hour 1 — the equation (no computer)

Read `burgers.py` lines 1–60 only. Then, on paper:

1. Write `u_t + u u_x = ν u_xx`.
2. Drop `ν` and explain each term out loud: "the fluid moves at speed u, so
   a point that is faster than the point behind it pulls away, and the profile
   steepens."
3. Derive `t_s = -1/min u₀'` from characteristics. It takes four lines:
   `dx/dt = u`, so along a characteristic `u` is constant and the slope obeys
   `du_x/dt = -u_x²`, which solves to `u_x = u₀'/(1 + t u₀')`. It diverges at
   `t = -1/u₀'`. For `u₀ = sin x` the worst point is `x = π`, where
   `u₀' = -1`, so `t_s = 1`.
4. Be able to say what `δ = 4ν/U` means: thin fluid ⇒ sharp shock, thick fluid
   ⇒ smeared shock.

**Check yourself:** can you explain, without notes, why the solution breaks at
`t = 1`, and why it does not break when `ν > 0`?

---

## Hour 2 — finite differences and stability

Read `advection`, `diffusion`, `step_explicit` in `burgers.py`.

1. Write down the three difference formulas (first derivative centred, second
   derivative centred). You have done these in class.
2. Work out `u[i] + dt*(-u[i]*(u[i+1]-u[i-1])/(2dx) + ν*(u[i+1]-2u[i]+u[i-1])/dx²)`.
   That is the whole explicit solver — one line.
3. Stability: the diffusion term adds `1 - 2ν dt/dx²` to the amplification
   factor of the highest mode; it must stay above `-1`, hence `dt < dx²/2ν`.
4. Run `python experiments/e1_sanity.py` and look at panel (c): at `0.9` of the
   limit the solution decays, at `1.1` it explodes within 200 steps.

**Check yourself:** why does halving `dx` make the run 4× longer, not 2×?

---

## Hour 3 — Crank–Nicolson and the Thomas algorithm

Read `_implicit_diffusion` and `step_cn`.

1. Implicit means: evaluate the diffusion at the *new* time level. That turns
   `u_new = u + dt ν D2 u_new` into a linear system `(I - r/2 D2) u_new = (I + r/2 D2) u`.
2. The matrix is tridiagonal, so `mylib.thomas_algorithm` solves it in `O(N)`
   operations — same cost per step as the explicit scheme, but with no `dt`
   restriction at all.
3. Why the halves? Splitting "diffusion then advection" once per step is only
   first order; sandwiching the advection between two *half* steps of diffusion
   restores second order. That is why panel (b) of fig. 1 reads `2.02, 2.00, 2.02`
   instead of `1.8, 1.6, 1.4`. (If sir asks, this is the one "trick" in the
   project, and you can prove it by running with `0.5*r` replaced by `r` in the
   first call.)

**Check yourself:** run `python -c "import burgers, numpy as np; x,dx=burgers.make_grid(256); print(burgers.solve('cn', N=256, nu=0.01, T=0.2, dt=50*dx**2/0.02)['nsteps'])"`.
The explicit scheme could never take that step. Say why.

---

## Hour 4 — the shock

Read `shock_width`, `max_slope`, then `experiments/e2_shock.py`.

1. Panel (a): the wave steepens, breaks at `t = 1`, then stops moving —
   the profile is frozen and only decays slowly.
2. Panel (b): before `t = 1` the strain follows `1/(t_s - t)`; after, viscosity
   caps it at `U²/2ν`. Measured 28.4 against 27.1.
3. Panel (c): the profile is a tanh. Verify by hand: substitute
   `u = -U tanh[U(x-x₀)/2ν]` into `u u_x = ν u_xx`; both sides give `U²/2`.
4. Panel (d): `δ = 4ν/U` — the jump divided by the maximum slope — measured at
   four viscosities, slope `1.00`.

**Check yourself:** why is `δ` measured 5 % low, and why is that the fault of
the grid rather than of the physics? (A discrete slope between two points
underestimates the true slope when the shock spans only a few cells.)

---

## Hour 5 — the Fourier part

Read `energy_spectrum` and `fit_power_law`, then `experiments/e3_spectrum.py`.

1. `c_k = (1/N) Σ u_j exp(-i k x_j)`, `E(k) = |c_k|²`. One line with numpy.
2. Why `k⁻²`: a sawtooth is a straight line with a jump. Its coefficients fall
   like `1/k` (integrate by parts once — the jump survives), so the energy
   falls like `1/k²`.
3. Panel (b) multiplies by `k²` to flatten the law; the knee is where viscosity
   eats the harmonics, and it sits at `k ≈ 1/δ`, which doubles when `ν` halves.
4. Be able to say: "viscosity does not change the slope, only where the line
   ends."

**Check yourself:** if the shock were smooth (Gaussian-like) instead of sharp,
would the slope be steeper or shallower? (Steeper — smoother functions have
faster-decaying coefficients.)

---

## Hour 6 — rehearse

Open `slides/slides.pdf` and talk through it twice with a timer: 8 slides,
13 minutes, leaving 2 minutes spare. Then answer the questions below out loud.

---

## Twelve questions you will be asked

1. **"Where is the Navier–Stokes equation here?"**
   In 1-D, incompressibility forces `u_x = 0`, which kills advection and the
   pressure gradient. Dropping the pressure gradient from the 1-D momentum
   equation leaves Burgers — the standard one-dimensional model of a
   compressible flow.

2. **"Why not a 2-D or 3-D Navier–Stokes?"**
   Time. In 1-D there is an exact steady shock (tanh) and an exact breaking
   time, so every number can be checked. In 3-D nothing is exact.

3. **"Is this chaos?"** No. There is no sensitive dependence: the shock
   position is fixed, and the equation is a stable relaxation to a steady
   state.

4. **"Why the boundary conditions?"** `u₀ = sin x` is odd about `x = 0`, and the
   equation preserves oddness, so the exact solution vanishes at both ends.
   Dirichlet zero is exact here, not an approximation.

5. **"What is the order of your scheme and how do you know?"** Second order.
   Fig. 1(b): error against a fine reference, slopes 2.02, 2.00, 2.02.

6. **"Why is Crank–Nicolson better than the explicit scheme?"** Diffusion is
   implicit, so the `dt < dx²/2ν` limit disappears. Same cost per step,
   far fewer steps. Fig. 1(c) shows the explicit scheme exploding just past
   the limit.

7. **"Why is there a Thomas algorithm here?"** The implicit step is a
   tridiagonal linear system; the Thomas algorithm solves it in `O(N)`.

8. **"You split the operators. Doesn't that ruin the accuracy?"** A plain split
   does (first order). Half-diffusion, full advection, half-diffusion is
   symmetric and stays second order — visible in the measured orders.

9. **"Why `k⁻²`?"** The profile is a sawtooth: straight lines joined by jumps.
   One integration by parts leaves the jump behind and gives `c_k ∝ 1/k`.

10. **"Where does the spectrum end?"** At `k ≈ 1/δ`, where `δ` is the shock
    thickness. Halving `ν` doubles the cut-off (31 → 63 in the figure).

11. **"Your thickness is 5 % low. Is your code wrong?"** No — the shock spans
    about seven grid cells and `max|Δu|/Δx` underestimates the true slope.
    The *scaling* `δ ∝ ν^1.000` is unaffected.

12. **"What would you do with more time?"** Smaller `ν` at higher resolution
    (a wider `k⁻²` range), and a forced version so the amplitude does not
    decay — but the physics needed for that is beyond this course.

---

## If you have only 20 minutes before the viva

Read: `t_s = 1`, `δ = 4ν/U`, the table in `README.md`, and questions 1, 5, 6, 9.
That covers 80 % of what will be asked.
