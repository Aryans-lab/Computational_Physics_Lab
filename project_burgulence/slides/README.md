# Presentation (15 + 5 minutes)

**`slides.pdf`** — 12 slides, 16:9, ready to show. This is the file to use.

| File | What it is |
|---|---|
| `slides.pdf` | the deck, compiled |
| `slides.typ` | its source (Typst) — edit and recompile with `typst compile slides.typ` |
| `deck.html` | same content as a self-contained browser deck: arrow keys / space to navigate, `p` to print |
| `slides.tex` | same content as beamer, if you prefer LaTeX |

The `.typ` sources use the STIX Two fonts by name; without them Typst falls back
to its default font and the layout shifts slightly. To get an identical PDF,
pass the fonts explicitly: `typst compile --font-path /path/to/stix slides.typ`.

Twelve frames for fifteen minutes is deliberate: about one minute per frame,
with three minutes on the grid-viscosity result.

| Slides | Minutes | Content |
|---|---|---|
| 1–3 | 2 | the equation, and the four formulas everything is tested against |
| 4 | 1 | methods in one slide |
| 5 | 3 | verification: orders, stability, cost |
| 6–8 | 5 | breaking, spectra, Burgulence |
| 9 | 2.5 | grid viscosity — the punchline |
| 10–12 | 1.5 | conclusions, limitations, references |

Leave the last five minutes for questions. The one that always comes is
*“how do you know your reference solution is right?”* — answer with the PDE
residual ($4\times10^{-9}$, `tests/test_core.py`) and the band-limited
interpolation of the initial data.
