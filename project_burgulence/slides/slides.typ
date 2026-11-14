// ---------------------------------------------------------------------------
//  DIY presentation -- 15 + 5 minutes, 16:9
//  Aryan Bandyopadhyay (2411014), PHY341/745
//
//  Typeset with Typst.  Same content as deck.html / slides.tex.
// ---------------------------------------------------------------------------

#set document(title: "Shock waves in a one-dimensional fluid — slides", author: "Aryan Bandyopadhyay")
#set page(width: 33.87cm, height: 19.05cm, margin: (x: 1.5cm, y: 1.1cm),
  fill: white, numbering: "1", number-align: right)
#set text(font: "STIX Two Text", size: 15pt, fill: rgb("#1b1b1b"))
#show math.equation: set text(font: "STIX Two Math", size: 15pt)
#set par(justify: false, leading: 0.65em)

#let accent = rgb("#1f4e79")
#let accent2 = rgb("#a8322d")
#let faint = rgb("#6b6b6b")
#let hilite(body) = text(fill: accent2, weight: 700, body)

#let slide(title, body, halign: center) = {
  grid(
    rows: (auto, 1fr),
    gutter: 0.4cm,
    {
      set text(size: 22pt, weight: 700, fill: accent)
      title
    },
    align(halign + horizon, body),
  )
}

// -------------------------------------------------------------- 1. title ---
#block(above: 2.6cm)[
  #text(size: 11pt, weight: 700, fill: accent2, tracking: 1.5pt)[DIY PROJECT#h(0.5em)·#h(0.5em)PHY341/745]
  #v(0.2cm)
  #text(size: 40pt, weight: 700, fill: accent)[Shock waves in a\ one-dimensional fluid]
  #v(0.1cm)
  #text(size: 20pt, fill: rgb("#3a3a3a"))[How much of what you see is the grid?]
  #v(0.85cm)
  #text(size: 14pt)[Aryan Bandyopadhyay#h(0.5em)·#h(0.5em)2411014]
  #v(0.12cm)
  #text(size: 11pt, fill: faint)[Physics Computer Lab#h(0.5em)·#h(0.5em)NISER Bhubaneswar]
]

#pagebreak()

// ------------------------------------------------------------ 2. equation ---
#slide([One equation], halign: left)[
  #grid(
    columns: (1fr, 1.05fr),
    gutter: 1.0cm,
    [
      #v(0.6cm)
      #align(center)[#text(size: 27pt)[$u_t + u u_x = nu u_(x x)$]]
      #v(0.8cm)
      #set text(size: 16pt)
      - nonlinearity *steepens* the wave \
      - viscosity *smooths* it \
      - the balance is a #hilite[shock]
      #v(0.7cm)
      #text(size: 12pt, fill: faint)[
        In one dimension, incompressibility forces $partial_x u = 0$: advection and
        the pressure gradient both vanish. Dropping only the pressure gradient
        leaves the *1-D Navier--Stokes* equation.]
    ],
    [
      #align(center)[#image("../figures/fig02_shock_formation.png", width: 100%)]
      #v(0.2cm)
      #align(center)[#text(size: 10.5pt, fill: faint)[The whole talk in one
      picture: steepening, breaking at $t = 1$, and the asymptotic profile.]]
    ],
  )
]

#pagebreak()

// --------------------------------------------------------------- 3. theory ---
#slide([The only theory needed], halign: left)[
  #grid(
    columns: (1fr, 1fr),
    gutter: 1.2cm,
    [
      #text(size: 15pt, weight: 700, fill: accent2)[When does it break?]
      #v(0.15cm)
      Characteristics $x = xi + u_0(xi) thin t$ cross at
      #align(center)[#text(size: 18pt)[$t_s = -1 / op("min") u_0'(x) = 1$]]
      #v(0.6cm)
      #text(size: 15pt, weight: 700, fill: accent2)[How thick is what remains?]
      #v(0.15cm)
      #align(center)[#text(size: 18pt)[$delta = 4 nu \/ U$]]
      #v(0.2cm)
      #align(center)[#text(size: 16pt)[$u = -U tanh[U (x - x_0) \/ 2 nu]$]]
    ],
    [
      #text(size: 15pt, weight: 700, fill: accent2)[Four predictions, no free parameters]
      #v(0.45cm)
      #align(center)[#text(size: 21pt)[$E(k) prop k^(-2)$]]
      #v(0.35cm)
      #align(center)[#text(size: 21pt)[$zeta_p = op("min")(p, 1)$]]
      #v(0.45cm)
      #text(size: 12.5pt, fill: faint)[
        Everything in this talk is measured against one of these. That is what
        makes the _numerics_ measurable rather than decorative.]
    ],
  )
]

#pagebreak()

// -------------------------------------------------------------- 4. methods ---
#slide([Five solvers, three diagnostics], halign: left)[
  #grid(
    columns: (1.05fr, 1fr),
    gutter: 1.2cm,
    [
      #text(size: 13pt)[
        #table(
          columns: (auto, auto),
          align: (left, left),
          stroke: none,
          inset: (x: 6pt, y: 3.5pt),
          table.hline(stroke: 0.7pt),
          [*scheme*], [*order*],
          table.hline(stroke: 0.4pt),
          [upwind], [1],
          [Lax--Wendroff], [2],
          [FTCS], [2 (unstable)],
          [Crank--Nicolson + Heun], [2],
          [spectral, Strang + RK4], [exponential],
          table.hline(stroke: 0.7pt),
        )
      ]
      #v(0.3cm)
      #text(size: 11.5pt, fill: faint)[
        All written for this project. NumPy only for arrays --- no `np.fft`,
        no `scipy`. The FFT is radix-2 and checked against a naive DFT to
        $10^(-15)$.]
    ],
    [
      #text(size: 15pt)[*Diagnostics*]
      #v(0.15cm)
      #text(size: 14.5pt)[
        $delta$ = jump / max slope \
        $E(k)$ from the from-scratch FFT \
        $S_p (r)$, moments of increments]
      #v(0.6cm)
      #text(size: 15pt)[*Reference*]
      #v(0.15cm)
      #text(size: 14.5pt)[exact Cole--Hopf solution, grid-independent to $10^(-12)$]
    ],
  )
]

#pagebreak()

// ----------------------------------------------------------- 5. validation ---
#slide([Result 1 --- are the solvers right?])[
  #image("../figures/fig01_validation.png", width: 97%)
  #v(0.1cm)
  #text(size: 13pt)[
    measured orders #hilite[0.94, 2.02, 2.00], spectral
    $7.5 times 10^(-4) -> 10^(-9)$; growth rates match von Neumann theory to
    #hilite[0.3%].]
]

#pagebreak()

// ------------------------------------------------------------- 6. breaking ---
#slide([Result 2 --- breaking, and what is left behind])[
  #image("../figures/fig03_shock_scaling.png", width: 97%)
  #v(0.1cm)
  #text(size: 13pt)[
    strain $= 1\/(t_s - t)$ to $0.5%$, saturating at #hilite[$247.5$ vs $250$];
    $delta prop nu^(0.985)$, prefactor $5$--$10%$ above $4 nu \/ U$.]
]

#pagebreak()

// -------------------------------------------------------------- 7. spectra ---
#slide([Result 3 --- the spectrum has a movable end])[
  #image("../figures/fig04_spectra.png", width: 97%)
  #v(0.1cm)
  #text(size: 13pt)[
    $E(k) prop k^(-2.030 plus.minus 0.001)$ for every $nu$; the cut-off
    $k_"cut" prop nu^(-0.971)$. Viscosity changes where the spectrum _ends_, not
    its shape.]
]

#pagebreak()

// ------------------------------------------------------------ 8. burgulence ---
#slide([Result 4 --- many shocks, one exponent])[
  #image("../figures/fig05_structure_functions.png", width: 93%)
  #v(0.1cm)
  #text(size: 13pt)[
    $zeta_p = 0.52, thin 0.92, thin 1.00, thin 0.94, thin 0.92$ versus
    $op("min")(p,1)$ and $p\/3$: a factor of three apart at $p = 4$, with a $5%$
    error bar.]
]

#pagebreak()

// -------------------------------------------------------- 9. grid viscosity ---
#slide([Result 5 --- how much of the viscosity is the grid?])[
  #image("../figures/fig06_grid_viscosity.png", width: 97%)
  #v(0.1cm)
  #text(size: 13pt)[
    $nu_"num" = U thin Delta x thin (1 - C) \/ 2$ verified to #hilite[1%] over
    five decades; at $N = 256$ the grid is #hilite[$3.7 times$] more viscous than
    the fluid.]
]

#pagebreak()

// ------------------------------------------------------------ 10. takeaways ---
#slide([What I take away], halign: left)[
  #set text(size: 16pt)
  - Convergence _order_ does not predict shock _quality_: all four schemes have
    known orders and disagree qualitatively --- smeared vs ringed.
  - Integral quantities are robust, local ones are not: the $k^(-2)$ slope
    survives a $10^3$ change in $nu$; the shock thickness does not.
  - Numerical dissipation is measurable, has a formula, and should be computed
    rather than hoped away.
  #v(0.55cm)
  #block(width: 100%, fill: rgb("#faf7f7"), stroke: (left: 3pt + accent2),
    inset: (x: 1em, y: 0.8em))[
    #set text(size: 15pt)
    *Practical criterion.* \ Resolve the shock, not the wavelength:
    $Delta x lt.eq 0.2 thin delta$ --- at least five points per shock.
  ]
]

#pagebreak()

// --------------------------------------------------------- 11. limitations ---
#slide([Limitations, and what I would do next], halign: left)[
  #grid(
    columns: (1fr, 1fr),
    gutter: 1.4cm,
    [
      #text(weight: 700, fill: accent2, size: 15pt)[Limitations]
      #v(0.3cm)
      #set text(size: 15pt)
      - inertial range only $approx 1.5$ decades
      - statistics at one time, not stationary
      - one equation; no vortex stretching (absent in 1-D)
    ],
    [
      #text(weight: 700, fill: accent2, size: 15pt)[Next]
      #v(0.3cm)
      #set text(size: 15pt)
      - forced, statistically stationary run
      - WENO/TVD scheme, same $zeta_p$ test
      - extension to 2-D
    ],
  )
]

#pagebreak()

// ---------------------------------------------------------- 12. references ---
#slide([References & code], halign: left)[
  #set text(size: 13.5pt)
  + J. M. Burgers, _Adv. Appl. Mech._ *1*, 171 (1948).
  + J. Bec and K. Khanin, _Phys. Rep._ *447*, 1 (2007).
  + G. B. Whitham, _Linear and Nonlinear Waves_ (Wiley, 1974), ch. 2--4.
  + R. J. LeVeque, _Finite Volume Methods for Hyperbolic Problems_ (CUP, 2002), ch. 4, 6.
  + J. P. Boyd, _Chebyshev and Fourier Spectral Methods_ (Dover, 2001), ch. 2.
  + J. D. Cole, _Quart. Appl. Math._ *9*, 225 (1951).
  #v(0.7cm)
  #block(width: 100%, fill: rgb("#faf7f7"), stroke: (left: 3pt + accent2),
    inset: (x: 1em, y: 0.8em))[
    #set text(size: 14pt)
    Code, figures and every number: \
    #text(weight: 700, fill: accent)[github.com/Aryans-lab/Computational_Physics_Lab]
    $->$ `project_burgulence/`
  ]
]
