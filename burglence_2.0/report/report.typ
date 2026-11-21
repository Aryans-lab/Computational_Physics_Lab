// ---------------------------------------------------------------------------
//  DIY project report -- PHY341/745, NISER Bhubaneswar
//  Aryan Bandyopadhyay (2411014)
//
//  Typeset with Typst:   typst compile report.typ
//  Figures: ../figures/*.png   (produced by ../run_all.py)
//
//  This is scaffolding with the correct numbers -- rewrite the prose.
// ---------------------------------------------------------------------------

#set document(title: "A shock in a one-dimensional fluid", author: "Aryan Bandyopadhyay")
#set page(paper: "a4", margin: (x: 1.8cm, top: 1.5cm, bottom: 1.5cm),
  numbering: "1")
#set text(font: "STIX Two Text", size: 9.8pt, lang: "en")
#show math.equation: set text(font: "STIX Two Math", size: 10pt)
#set math.equation(numbering: "(1)")
#set par(justify: true, leading: 0.63em, spacing: 0.66em)
#set heading(numbering: "1.")
#show heading: it => block(above: 1.0em, below: 0.5em, {
  set text(size: 10.8pt, weight: 700, fill: rgb("#1f4e79"))
  it
})
#show heading.where(level: 2): it => block(above: 0.7em, below: 0.3em, {
  set text(size: 9.9pt, weight: 700, fill: rgb("#1f4e79"))
  it
})
#set figure.caption(separator: [ ])
#show figure.caption: set text(size: 8.4pt)
#show figure: set block(breakable: false, above: 0.6em, below: 0.6em)
#show ref: it => text(fill: rgb("#1f4e79"), it)

#align(center)[
  #text(size: 15pt, weight: 700, fill: rgb("#1f4e79"))[A shock in a one-dimensional fluid]
  #v(-0.45em)
  #text(size: 10.5pt)[Steepening, breaking and the $k^(-2)$ spectrum of the viscous Burgers equation]
  #v(-0.25em)
  #text(size: 8.8pt)[Aryan Bandyopadhyay (2411014) \ — PHY341/745 Physics Computer Lab, NISER Bhubaneswar]
]
#v(-0.5em)
#line(length: 100%, stroke: 0.6pt + rgb("#1f4e79"))

#pad(left: 1.5em, right: 1.5em)[
  #text(size: 8.9pt)[
    *Abstract.* In one spatial dimension the incompressibility condition removes
    both the advection term and the pressure gradient from the Navier--Stokes
    equation, leaving the Burgers equation $u_t + u u_x = nu u_(x x)$ as the only
    non-trivial one-dimensional model of a compressible flow. It contains two
    competing effects: the nonlinear term steepens a smooth wave until it breaks,
    and viscosity replaces the resulting discontinuity by a front of finite
    thickness. We solve the equation two ways --- an explicit finite-difference
    scheme and a Crank--Nicolson scheme built on the Thomas algorithm --- check
    both against exact results, and Fourier-analyse the final shock. The
    breaking time, the $tanh$ shock profile and its $4 nu \/ U$ thickness are
    reproduced to within a few per cent, and the energy spectrum of the shock
    train is measured to be $E(k) prop k^(-2)$.
  ]]
#v(-0.2em)

= Introduction

The Burgers equation
$ u_t + u u_x = nu u_(x x), quad 0 <= x <= 2 pi, quad u(0,t) = u(2 pi,t) = 0 $ <eq:burgers>
describes a fluid in which every point moves at its own speed $u$. Where the
velocity decreases with $x$, the point ahead moves slower than the point
behind, so the profile steepens: the slope $u_x$ obeys $"d"u_x\/"d"t = -u_x^2$
along a characteristic and therefore
$ u_x (t) = u_0'(x) / (1 + t thin u_0'(x)) , $ <eq:slope>
which diverges at
$ t_s = -1 / op("min") u_0'(x) = 1 quad "for" quad u_0 = sin x . $ <eq:ts>
After $t_s$ the inviscid solution is multivalued. Viscosity does not allow a
discontinuity, and replaces it by a steady front
$ u = -U tanh[U (x - x_0) \/ 2 nu] $ <eq:tanh>
of thickness $delta = 4 nu \/ U$, where $U$ is half the jump across it. Both the breaking time and the shock profile are checked below.

The boundary conditions are exact rather than convenient: $u_0 = sin x$ is odd
about $x = 0$, and @eq:burgers preserves that symmetry, so the true solution
vanishes at $x = 0$ and $x = 2 pi$ for all time.

A shock train is a sawtooth, and the Fourier coefficients of a sawtooth fall
off as $1\/k$, so its energy spectrum should obey $E(k) prop k^(-2)$ over the
range between the box and the shock thickness [3]. Measuring that exponent is
the second half of the project.

= Numerical method

== Two schemes

The interval is divided into $N$ cells of width $Delta x$, and the equation is
marched in time. Writing the right-hand side as
$F(u) = -u u_x + nu u_(x x)$ with centred differences for both derivatives, the
explicit scheme is simply
$ u^(n+1) = u^n + Delta t thin F(u^n) . $
It is stable only if $Delta t < min(Delta x \/ |u|, thin Delta x^2 \/ 2 nu)$;
the second condition is severe, because halving $Delta x$ costs four times as
many steps.

The second scheme treats the diffusion implicitly. One step is split
symmetrically,
$ "diffusion(" Delta t \/ 2 ")" quad -> quad "advection(" Delta t ")" quad -> quad "diffusion(" Delta t \/ 2 ")" , $
and each diffusion half-step solves the tridiagonal system
$ (I - r/2 thin D_2) thin u^(n+1) = (I + r/2 thin D_2) thin u^n, quad r = nu Delta t \/ Delta x^2 , $ <eq:cn>
with `mylib.thomas_algorithm`. Because the diffusion is implicit, the
$Delta x^2 \/ 2 nu$ restriction disappears. The symmetric splitting matters:
stepping diffusion and advection once each would be only first order in time,
whereas the halved version keeps second order.

== Tools

Only three routines are borrowed from the course library:
`thomas_algorithm` (the tridiagonal solve in @eq:cn, once per half-step),
`poly_fit` (the log--log slope fits) and `pde_heat_1d_crank_nicolson` (used in
the pure-diffusion check). NumPy supplies array arithmetic and the FFT used for
the spectrum; matplotlib draws the figures. Everything else --- the schemes,
the stability limit, the diagnostics --- is written in `burgers.py`.

= Results

== Is the code right?

@fig:sanity collects three independent checks. Panel (a) switches the
nonlinear term off, so the equation becomes $u_t = nu u_(x x)$ with the exact
solution $exp(-nu t) sin x$; the numerical profile matches it to
$8.0 times 10^(-6)$. Panel (b) compares runs at $N = 64, 128, 256, 512$ with a
reference run at $N = 4096$, and the errors fall with measured orders $2.02$,
$2.00$ and $2.02$. Panel (c) is the stability boundary: at $0.9$ of the
explicit limit the solution decays smoothly, at $1.1$ it blows up within a few
hundred steps, exactly where von Neumann theory places the boundary.

#figure(
  image("../figures/fig01_sanity.png", width: 100%),
  caption: [(a) The diffusion step against the exact solution.
  (b) Convergence of the full equation before breaking, against a fine
  reference. (c) The explicit scheme either side of $Delta t = Delta x^2 \/ 2 nu$.],
) <fig:sanity>

== Steepening and the shock

@fig:shock (a) shows profiles at five times for $nu = 10^(-2)$. The wave
steepens, breaks at $t = 1$ and then stops: the asymptotic state is a
sawtooth, straight lines joined by a sharp front. Panel (b) measures the
steepest slope. Before breaking it follows $1\/(t_s - t)$ and after breaking it
saturates at $U^2 \/ 2 nu$: at the final time the measured value is $28.4$
against a prediction of $27.1$ from the instantaneous amplitude. Panel (c)
compares the front with @eq:tanh and the two are indistinguishable to within
$2%$ of $U$. Panel (d) repeats the measurement at four viscosities: the
thickness obeys $delta prop nu^(1.00)$ and sits a few per cent below the
prediction $4 nu \/ U$, the deficit shrinking as the shock becomes better
resolved. At $nu = 10^(-2)$ the measured thickness is $0.0518$ against
$0.0544$.

#figure(
  image("../figures/fig02_shock.png", width: 100%),
  caption: [(a) Steepening, breaking and the asymptotic sawtooth.
  (b) Maximum strain against time, with the inviscid divergence and the viscous
  saturation. (c) The front at $t = 3$ against $-U tanh[U(x-pi)\/2 nu]$.
  (d) Shock thickness against viscosity.],
) <fig:shock>

== The spectrum of a shock

@fig:spectrum (a) shows $E(k)$ for two viscosities on log--log axes. Both are
straight lines over more than a decade, with fitted slopes $-2.06$ and $-2.03$
in the window $4 <= k <= 12$. Viscosity does not change the slope; it changes
where the line ends. Panel (b) makes this visible by plotting $E(k) k^2$,
which is flat while the $k^(-2)$ law holds and falls once viscosity removes the
harmonics. The knee moves from $k approx 31$ at $nu = 10^(-2)$ to
$k approx 63$ at $nu = 5 times 10^(-3)$, i.e. as $1\/delta$.

#figure(
  image("../figures/fig03_spectrum.png", width: 96%),
  caption: [(a) Energy spectra for two viscosities, with the $k^(-2)$ line.
  (b) The same data multiplied by $k^2$: flat while the power law holds, and a
  knee at the shock thickness.],
) <fig:spectrum>

== Summary of the measurements

#align(center)[
  #table(
    columns: (auto, auto, auto),
    align: (left, center, center),
    stroke: none,
    inset: (x: 8pt, y: 2.4pt),
    table.hline(stroke: 0.7pt),
    [*quantity*], [*prediction*], [*measured*],
    table.hline(stroke: 0.4pt),
    [breaking time $t_s$], [$1$], [$1.00$],
    [peak strain, $nu = 10^(-2)$], [$U^2 \/ 2 nu = 27.1$], [$28.4$],
    [thickness, $nu = 10^(-2)$], [$4 nu \/ U = 0.0544$], [$0.0518$],
    [scaling of $delta$], [$nu^1$], [$nu^(1.00)$],
    [deviation from $tanh$], [$0$], [$2%$ of $U$],
    [spectral slope], [$-2$], [$-2.06$, $-2.03$],
    [cut-off $k_"cut"$], [$prop 1 \/ delta$], [$31 -> 63$ as $nu$ halves],
    table.hline(stroke: 0.7pt),
  )
]

= Discussion

Three things are worth keeping from this. First, every quantity in the project
has an exact value, so the code can be wrong only in ways that show up
immediately --- the breaking time, the thickness and the spectral slope are not
matters of opinion. Second, the two schemes differ only in cost: both are
second order, but the implicit one is not limited by $Delta t < Delta x^2\/2nu$
and reaches the same answer in far fewer steps. Third, the Fourier view
separates two things that look merged in real space: the slope of the spectrum
is a property of the shock shape, while the position of the cut-off is a
property of the viscosity.

The limitations are equally clear. The thickness comes out a few per cent low
because the shock spans only about seven grid cells and a discrete maximum
slope underestimates the true one. The spectral range covers one decade, which
is enough to tell $k^(-2)$ from $k^(-3)$ and nothing more. And the boundaries
let the wave decay, so the amplitude is not constant and the saturation value
$U^2\/2nu$ has to be evaluated with the instantaneous $U$.

= Conclusion

The one-dimensional Navier--Stokes equation is small enough to be solved in
twenty lines and rich enough to contain a shock, an exact profile, a finite
breaking time and a clean power law. Two finite-difference schemes, an exact
solution for each regime, and one FFT are enough to reproduce all of it.

= Code

Four files, about 700 lines in total. `burgers.py` holds the grid, the initial
condition, both solvers and the three measurements (`max_slope`,
`shock_width`, `energy_spectrum`); the three scripts in `experiments/` each
produce one figure and write their numbers to `results/*.json`.

#align(center)[
  #table(
    columns: (auto, auto),
    align: (left, left),
    stroke: none,
    inset: (x: 8pt, y: 2.4pt),
    table.hline(stroke: 0.7pt),
    [*file*], [*contents*],
    table.hline(stroke: 0.4pt),
    [`burgers.py`], [grid, solvers, diagnostics],
    [`experiments/e1_sanity.py`], [figure 1: accuracy and stability],
    [`experiments/e2_shock.py`], [figure 2: the shock],
    [`experiments/e3_spectrum.py`], [figure 3: the spectrum],
    [`tests/test_burgers.py`], [eleven checks, 0.3 s],
    table.hline(stroke: 0.7pt),
  )
]

`python run_all.py` regenerates every figure and every number quoted above in
about forty-five seconds; `python tests/test_burgers.py` verifies the solver,
the Thomas algorithm, the $tanh$ formulas and the spectral slope.

#v(0.3em)
#line(length: 100%, stroke: 0.4pt)
#v(-0.3em)
#text(size: 8.6pt)[
  *References.*
  #par(hanging-indent: 1.1em)[
  [1] J. M. Burgers, "A mathematical model illustrating the theory of
  turbulence", _Adv. Appl. Mech._ *1*, 171 (1948). \
  [2] G. B. Whitham, _Linear and Nonlinear Waves_ (Wiley, 1974), ch. 2--4. \
  [3] J. Bec and K. Khanin, "Burgers turbulence", _Phys. Rep._ *447*, 1 (2007). \
  ]
]
