// ---------------------------------------------------------------------------
//  DIY presentation -- 15 + 5 minutes, 16:9, eight slides
//  Aryan Bandyopadhyay (2411014), PHY341/745
//
//  typst compile slides.typ      (figures in ../figures/)
// ---------------------------------------------------------------------------

#set document(title: "A shock in a one-dimensional fluid — slides", author: "Aryan Bandyopadhyay")
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
#block(above: 2.4cm)[
  #text(size: 11pt, weight: 700, fill: accent2, tracking: 1.5pt)[DIY PROJECT#h(0.5em)·#h(0.5em)PHY341/745]
  #v(0.2cm)
  #text(size: 38pt, weight: 700, fill: accent)[A shock in a\ one-dimensional fluid]
  #v(0.1cm)
  #text(size: 19pt, fill: rgb("#3a3a3a"))[Steepening, breaking, and a $k^(-2)$ spectrum]
  #v(0.8cm)
  #text(size: 14pt)[Aryan Bandyopadhyay#h(0.5em)·#h(0.5em)2411014]
  #v(0.12cm)
  #text(size: 11pt, fill: faint)[Physics Computer Lab#h(0.5em)·#h(0.5em)NISER Bhubaneswar]
]

#pagebreak()

// ----------------------------------------------------------- 2. equation ---
#slide([One equation, two keywords], halign: left)[
  #grid(
    columns: (1fr, 1.05fr),
    gutter: 1.0cm,
    [
      #align(center)[#text(size: 26pt)[$u_t + u u_x = nu u_(x x)$]]
      #v(0.5cm)
      #text(size: 15pt)[
        *1-D Navier–Stokes* (primary) #h(0.4em)·#h(0.4em) *Fourier analysis* (secondary)
      ]
      #v(0.4cm)
      #text(size: 13.5pt)[
        - nonlinearity *steepens* the wave
        - viscosity stops it becoming a discontinuity
        - the result is a *shock* of thickness $delta = 4 nu \/ U$
      ]
      #v(0.4cm)
      #text(size: 12pt, fill: faint)[
        In one dimension, incompressibility kills advection and pressure; what
        is left is Burgers.]
    ],
    [
      #align(center)[#image("../figures/fig02a_profiles.png", width: 100%)]
      #v(0.15cm)
      #align(center)[#text(size: 10pt, fill: faint)[The whole project on one page:
      steepening, breaking at $t = 1$, then a steady shock.]]
    ],
  )
]

#pagebreak()

// ------------------------------------------------------------ 3. methods ---
#slide([Two solvers, one exact answer each], halign: left)[
  #grid(
    columns: (1.15fr, 1fr),
    gutter: 1.2cm,
    [
      #text(size: 14pt)[
        #table(
          columns: (auto, auto, auto),
          align: (left, left, left),
          stroke: none,
          inset: (x: 6pt, y: 3.5pt),
          table.hline(stroke: 0.7pt),
          [*scheme*], [*order*], [*limit*],
          table.hline(stroke: 0.4pt),
          [explicit], [2], [$Delta t < Delta x^2 \/ 2 nu$],
          [Crank–Nicolson], [2], [none],
          table.hline(stroke: 0.7pt),
        )
      ]
      #v(0.35cm)
      #text(size: 13.5pt)[
        CN = half diffusion step, full advection step, half diffusion step.
        The diffusion is implicit, so each step costs one tridiagonal solve ---
        `mylib.thomas_algorithm` --- and the $Delta x^2$ limit disappears.
      ]
    ],
    [
      #text(size: 14.5pt)[*Checked against*]
      #v(0.2cm)
      #text(size: 14.5pt)[
        - $u_t = nu u_(x x)$: exact $exp(-nu t) sin x$
        - before breaking: second-order convergence
        - after breaking: $tanh$ shock, $delta = 4 nu \/ U$
        - spectrum: $E(k) prop k^(-2)$
      ]
      #v(0.3cm)
      #text(size: 11.5pt, fill: faint)[
        `mylib`: `thomas_algorithm`, `poly_fit`, `pde_heat_1d_crank_nicolson`.]
    ],
  )
]

#pagebreak()

// --------------------------------------------------------- 4. is it right ---
#slide([Is the code right?])[
  #image("../figures/fig01_sanity.png", width: 97%)
  #v(0.1cm)
  #text(size: 13pt)[
    diffusion: error #hilite[$8 times 10^(-6)$]  ·  convergence orders
    #hilite[2.02, 2.00, 2.02]  ·  explicit scheme dies just past
    $Delta x^2 \/ 2 nu$.]
]

#pagebreak()

// ------------------------------------------------------------- 5. shock ---
#slide([The shock])[
  #image("../figures/fig02_shock.png", height: 11.4cm)
  #v(0.1cm)
  #text(size: 12.5pt)[
    breaks at $t_s = 1$  ·  peak strain #hilite[28.4] against $U^2\/2nu = 27.1$
    ·  $delta = 0.0518$ against $4nu\/U = 0.0544$  ·  $delta prop nu^1$ measured
    as #hilite[$nu^(1.00)$].]
]

#pagebreak()

// ---------------------------------------------------------- 6. spectrum ---
#slide([The spectrum of a shock])[
  #image("../figures/fig03_spectrum.png", width: 92%)
  #v(0.1cm)
  #text(size: 13pt)[
    slope #hilite[$-2.06$ and $-2.03$] for $nu = 10^(-2)$ and $5 times 10^(-3)$;
    the cut-off moves $31 -> 63$. Viscosity changes where the line *ends*,
    not its slope.]
]

#pagebreak()

// --------------------------------------------------------- 7. takeaway ---
#slide([What it means], halign: left)[
  #set text(size: 15.5pt)
  - Two second-order schemes, but only one of them survives a large time step.
  - Every number in the project is predicted in advance: $t_s$, $delta$,
    the $tanh$ shape, the $k^(-2)$ slope.
  - Fourier separates two things that look merged in real space: the *slope*
    is the shock shape, the *cut-off* is the viscosity.
  #v(0.45cm)
  #block(width: 100%, fill: rgb("#faf7f7"), stroke: (left: 3pt + accent2),
    inset: (x: 1em, y: 0.7em))[
    #set text(size: 13.5pt)
    *Limitations.* One equation, one shock; thickness 5 % low because the shock
    spans seven cells; the spectrum covers one decade.
  ]
]

#pagebreak()

// ------------------------------------------------------- 8. references ---
#slide([References], halign: left)[
  #set text(size: 13.5pt)
  + J. M. Burgers, "A mathematical model illustrating the theory of turbulence",
    _Adv. Appl. Mech._ *1*, 171 (1948).
  + G. B. Whitham, _Linear and Nonlinear Waves_ (Wiley, 1974), ch. 2–4.
  + J. Bec and K. Khanin, "Burgers turbulence", _Phys. Rep._ *447*, 1 (2007).
  #v(0.6cm)
  #text(size: 12.5pt, fill: faint)[
    Code, figures and every number:
    #text(fill: accent, weight: 700)[github.com/Aryans-lab/Computational\_Physics\_Lab]
    → `burglence_2.0/`]
  #v(0.6cm)
  #align(center)[#text(size: 26pt, weight: 700, fill: accent)[Thank you — questions?]]
]
