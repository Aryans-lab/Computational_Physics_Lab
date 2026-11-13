# Presentation (15 + 5 minutes)

Two decks with identical content — use whichever suits you.

| File | How to use |
|---|---|
| `deck.html` | Open in any browser. Arrow keys / space to navigate, `p` to print (then “Save as PDF”, landscape, no margins). No installation, no network. |
| `slides.tex` | `pdflatex slides.tex` (twice). Beamer, 16:9, figures pulled from `../figures/`. |

Twelve frames for fifteen minutes is deliberate: roughly one minute per frame
plus three minutes on the grid-viscosity result.

Suggested timing:

| Frames | Minutes | Content |
|---|---|---|
| 1–3 | 2 | what the equation is, and the two formulas everything is tested against |
| 4 | 1 | methods in one slide |
| 5 | 3 | verification: orders, stability, cost |
| 6–8 | 5 | breaking, spectra, Burgulence |
| 9 | 2.5 | grid viscosity — the punchline |
| 10–12 | 1.5 | conclusions, limitations, references |

Leave the last five minutes for questions. The one that always comes is
*“how do you know your reference solution is right?”* — answer with the PDE
residual ($4\times10^{-9}$, `tests/test_core.py`) and the band-limited
interpolation of the initial data.
